from __future__ import division
from itertools import islice, imap
from operator import itemgetter

import argparse
import casacore.tables as ct
import numpy as np
import math
import subprocess
import os
import sys

def form_baselines(antennas):
    ''' Calculate all possible baselines.
    Args:
        antennas (list) : a list containing all antennas as an integer.
    '''
    baselines = []
    for i in antennas:
        for j in antennas:
            if j <= i:
                # It's either the same antenna or the baseline is already covered.
                continue
            else:
                baselines.append((i,j))
    return baselines

desc = '''Process the visibilities in an MS file in order to determine better estimates of the sigmas and weights for them.\nCASA assigns the same value to each baseline as 1  / sqrt(BW*T) whereas this script will estimate the sigmas for each baseline by looking at the scatter in the visibilities through time. The weights are then calculated as 1 / sigma**2.'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('--backup', action='store_true', dest='backup', help='create a backup of the MS file before operating on it')
parser.add_argument('--subtract', action='store_true', dest='subtract', help='calculate with the subtracted visibilities instead')
parser.add_argument('--hdf5', action='store_true', dest='use_hdf5', help='write out a Pandas dataframe in HDF5 format instead of a txt file')
parser.add_argument('filename', action='store')
args = parser.parse_args()
SUBTRACT = args.subtract
filename = args.filename
BACKUP = args.backup
HDF5 = args.use_hdf5
if SUBTRACT:
    print '[MaSER] Using the subtracted visibilties to calculate sigmas and weights.'
else:
    print '[MaSER] Using the regular visibilties to calculate sigmas and weights.'

if BACKUP:
    print '[MaSER] Backing up MS file...'
    subprocess.call('cp -r' + ' ' + filename + ' ' + filename + '.BACKUP', shell=True)
if HDF5:
    import pandas as pd

'''
Read in MS file and create table object.
'''
print '[MaSER] Reading in file...'
msfile = ct.table(filename, readonly=False)

'''
Calculate the number of baselines present in the measurement set.
'''
print '[MaSER] Forming baselines...'
ANTENNA1 = msfile.getcol('ANTENNA1')
ANTENNA2 = msfile.getcol('ANTENNA2')

ANTENNAS = set(ANTENNA1).union(set(ANTENNA2))
baselines = form_baselines(ANTENNAS)

print len(ANTENNA1)
print len(ANTENNA2)
print len(baselines)

print 'Formed %d baselines.' % len(baselines)

# Determine the number of correlations present.
polarizations = ct.taql('SELECT FROM %s::POLARIZATION'%(filename))
correlations = polarizations.getcol('NUM_CORR')
print 'Found %d correlations.' % correlations[0]

'''
Extract the visibilities per baseline. The measurement set is accessed using TaQL extracting
the columns into numpy arrays.
The final data is stored in regular text files, with each baseline and correlation having its
own file.
The columns will have the uvw coordinates, channel frequency, real and imaginary parts of the
visibilities and their corresponding errors.
'''
print '[MaSER] Extracting visibilities per baseline...'
try:
    os.mkdir('visibilities')
except:
    pass

for corr in range(correlations):
    if HDF5:
        df = None
    print 'Processing correlation %d/%d:' % (corr+1, correlations[0])
    # Select the spectral window keyword from the main table.
    spws = ct.taql('SELECT FROM %s::SPECTRAL_WINDOW'%(filename))
    Nspw = len(spws)
    # Frequencies corresponding to each channel.
    frequencies = spws.getcol('CHAN_FREQ')
    frequencies_flat_ghz = frequencies.flatten() * 1e-9
    # Calculate channel frequencies.
    nu = frequencies_flat_ghz
    # Determine channels in a spectral window
    Nchan = len(nu) // Nspw

    progress = 0; end = len(baselines); printed = False
    for k,(ant1,ant2) in enumerate(baselines):
        p = int(progress / end * 100)
        if ((p % 10) == 0) and not printed:
            print '%d%%' % (p),
            printed = True
        elif printed and (p%10 == 9):
            printed = False
        sys.stdout.flush()
        # Select wanted columns (possibly slightly redundant step).
        baseline = ct.taql('SELECT UVW,DATA FROM $msfile WHERE ANTENNA1=$ant1 AND ANTENNA2=$ant2')
        data = baseline.getcol('DATA')
        # uvw is an array with three values: the u, v and w coordinates.
        uvw = baseline.getcol('UVW')
        # Split the u,v,w coordinates for each baseline in meters.
        uu, vv, ww = uvw[...,0], uvw[...,1], uvw[...,2]
        u = np.asarray([i for i in uu for x in xrange(Nchan)])
        v = np.asarray([i for i in vv for x in xrange(Nchan)])
        w = np.asarray([i for i in ww for x in xrange(Nchan)])

        # Select the correlation, specified by the last index.
        # data has the shape (timestamps, channels, correlations) where timestamps repeats after each spectral window.
        subdata = data[...,corr]
        # Flatten data into subdata changing the shape from (timestamps, channels) to (timestamps*channels)
        subdata = subdata.flatten()
        # Take the real and imaginary parts of the data and handle them separately.
        data_real = subdata.real
        data_imag = subdata.imag

        if not SUBTRACT:
            # Calculate standard deviations.
            stdr = data_real.std()
            stdi = data_imag.std()
            # Because we have multiple channels, but only one standard deviation per spectral window we need to pad the array to match the length of the frequencies.
            std_real = np.zeros(len(data_real)); std_real.fill(stdr)
            std_imag = np.zeros(len(data_real)); std_imag.fill(stdi)
            # The MS file only has one sigma per correlation, so take the largest.
            sigma = max(stdr, stdi)
            # We need correct frequencies written out for the uv coordinates. Here we determine when we should switchs to a new spw.
            Nswitch = (len(uu) // Nspw)
            freq = []
            for i in range(Nspw):
                f = list(nu[Nchan * i:Nchan * (i+1)]) * Nswitch
                freq.extend(f)
            freq = np.asarray(freq)
            # Save data to file.
            #print len(u), len(v), len(w), len(data_real), len(data_imag), len(std_real), len(std_imag), len(freq)
            if not HDF5:
                FILEHEADER = 'Baseline: %d-%d\nEntries: %d\nu [m], v [m], w [m], frequency [GHz], real, imag, std(real), std(imag)' % (ant1, ant2, len(data_real))
                with open('visibilities/visibilities_corr_%.2d.txt'%(corr,), 'ab') as f:
                    np.savetxt(f, zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), header=FILEHEADER)
            else:
                try:
                    tdf = pd.DataFrame(zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), columns=['u', 'v', 'w', 'freq', 'data_real', 'data_imag', 'sigma_real', 'sigma_imag'])
                    df = df.append(tdf, ignore_index=True)
                except:
                    df = pd.DataFrame(zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), columns=['u', 'v', 'w', 'freq', 'data_real', 'data_imag', 'sigma_real', 'sigma_imag'])
            weights = sigma ** -2
            ct.taql('UPDATE $msfile SET SIGMA[$corr]=$sigma WHERE ANTENNA1=$ant1 AND ANTENNA2=$ant2')
            ct.taql('UPDATE $msfile SET WEIGHT[$corr]=$weights WHERE (ANTENNA1=$ant1 AND ANTENNA2=$ant2)')
            #ct.taql('UPDATE $msfile SET WEIGHT_SPECTRUM[$corr]=$weights WHERE (ANTENNA1=$ant1 AND ANTENNA2=$ant2)')
        elif SUBTRACT:
            real_shifted = np.roll(data_real, -1)
            imag_shifted = np.roll(data_imag, -1)
            sub_real = real_shifted - data_real
            sub_imag = imag_shifted - data_imag
            # Calculate new standard deviations with sqrt(2) correction, as the number of visibilities goes down by a factor 2.
            stdr = sub_real.std() / np.sqrt(2)
            stdi = sub_imag.std() / np.sqrt(2)
            # Because we have multiple channels, but only one standard deviation per spectral window we need to pad the array to match the length of the frequencies.
            std_real = np.zeros(len(data_real)); std_real.fill(stdr)
            std_imag = np.zeros(len(data_real)); std_imag.fill(stdi)
            # The MS file only has one sigma per correlation, so take the largest.
            sigma_sub = max(stdr, stdi)
            # We need correct frequencies written out for the uv coordinates. Here we determine when we should switchs to a new spw.
            Nswitch = (len(uu) // Nspw)
            freq = []
            for i in range(Nspw):
                f = list(nu[Nchan * i:Nchan * (i+1)]) * Nswitch
                freq.extend(f)
            freq = np.asarray(freq)
            # The subtracted visibilities.
            if not HDF5:
                FILEHEADER = 'Baseline: %d-%d\nEntries: %d\nu [m], v [m], w [m], frequency [GHz], real, imag, std(real), std(imag)' % (ant1, ant2, len(data_real))
                with open('visibilities/visibilities_subtracted_corr_%.2d.txt'%(corr,), 'ab') as f:
                    np.savetxt(f, zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), header=FILEHEADER)
            else:
                try:
                    tdf = pd.DataFrame(zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), columns=['u', 'v', 'w', 'freq', 'data_real', 'data_imag', 'sigma_real', 'sigma_imag'])
                    df = df.append(tdf, ignore_index=True)
                except:
                    df = pd.DataFrame(zip(u, v, w, freq, data_real, data_imag, std_real, std_imag), columns=['u', 'v', 'w', 'freq', 'data_real', 'data_imag', 'sigma_real', 'sigma_imag'])
            # Write back errors and weights to the SIGMA and WEIGHT columns of the MS file.
            weights = sigma_sub ** -2
            test =  np.zeros(2)
            ct.taql('UPDATE $msfile SET SIGMA[$corr]=$sigma_sub WHERE ANTENNA1=$ant1 AND ANTENNA2=$ant2')
            ct.taql('UPDATE $msfile SET WEIGHT[$corr]=$weights WHERE (ANTENNA1=$ant1 AND ANTENNA2=$ant2)')
            #ct.taql('UPDATE $msfile SET WEIGHT_SPECTRUM[$corr]=$weights WHERE (ANTENNA1=$ant1 AND ANTENNA2=$ant2)')
        progress += 1
    if HDF5 and not SUBTRACT:
        df.to_hdf('./visibilities/vis_corr_%.2d.hdf5'%corr, 'df', mode='w', format='table', data_columns=True)
    elif HDF5 and SUBTRACT:
        df.to_hdf('./visibilities/vis_subtracted_corr_%.2d.hdf5'%corr, 'df', mode='w', format='table', data_columns=True)
    print '100%\n'

print '[MaSER] Closing MS file...'
msfile.close()
print '[MaSER] Finished'
