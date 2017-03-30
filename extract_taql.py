from __future__ import division
from itertools import islice, imap
from operator import itemgetter

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

'''
Read in MS file and create table object.
'''
print '[CVE] Reading in file...'
filename = sys.argv[1]#'TESTVIS.ms'
msfile = ct.table(filename, readonly=False)

'''
Calculate the number of baselines present in the measurement set.
'''
print '[CVE] Forming baselines...'
ANTENNA1 = msfile.getcol('ANTENNA1')
ANTENNA2 = msfile.getcol('ANTENNA2')

ANTENNAS = set(ANTENNA1).union(set(ANTENNA2))
baselines = form_baselines(ANTENNAS)
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
print '[CVE] Extracting visibilities per baseline...'
try:
    os.mkdir('visibilities')
except:
    pass

for corr in range(correlations):
    print 'Processing correlation %d/%d:' % (corr+1, correlations[0])
    # Select the spectral window keyword from the main table.
    spws = ct.taql('SELECT FROM %s::SPECTRAL_WINDOW'%(filename))
    # Frequencies corresponding to each channel.
    frequencies = spws.getcol('CHAN_FREQ')
    frequencies_flat_ghz = frequencies.flatten() * 1e-9
    # Calculate channel frequencies.
    nu = frequencies_flat_ghz

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
        u, v, w = uvw[...,0], uvw[...,1], uvw[...,2]
        # Select the correlation, specified by the last index.
        # data has the shape (timestamps, channels, correlations).
        subdata = data[...,corr]
        # subdata now has the shape (timestamps, channels)
        subdata = subdata.flatten()
        # Take the real and imaginary parts of the data and handle them separately.
        data_real = subdata.real
        data_imag = subdata.imag
        # Calculate standard deviations.
        stdr = data_real.std()
        stdi = data_imag.std()
        # Because we have multiple channels, but only one standard deviation per spectral window we need to pad the array to match the length of the frequencies.
        std_real = np.zeros(len(nu)); std_real.fill(stdr)
        std_imag = np.zeros(len(nu)); std_imag.fill(stdi)
        # The MS file only has one sigma per correlation, so take the largest.
        sigma = max(stdr, stdi) 
        '''
        # Subtract every first visibility from the second, i.e. 2-1, 4-3, 6-5 etc.
        i = 0
        sub_real = []
        sub_imag = []
        while i < (len(data_real)//2 - 1):
            subr = data_real[2*i+1] - data_real[2*i]
            subi = data_imag[2*i+1] - data_imag[2*i]
            sub_real.append(subr)
            sub_imag.append(subi)
            i += 1
        sub_real = np.asarray(sub_real)
        sub_imag = np.asarray(sub_imag)
        # Calculate new standard deviations.
        stds_real = sub_real.std()
        stds_imag = sub_imag.std()
        # Padd the arrays and write out to file.
        std_sub_real = np.zeros(len(nu)); std_sub_real.fill(stds_real)
        std_sub_imag = np.zeros(len(nu)); std_sub_imag.fill(stds_imag)
        # The subtracted visibilities.
        FILEHEADER = 'Baseline: %d-%d\nEntries: %d\nu [m], v [m], w [m], frequency [GHz], real, imag, std(real), std(imag)' % (ant1, ant2, nu.shape[0])
        with open('visibilities/visibilities_subtracted_corr_%.2d.txt'%(corr,), 'ab') as f:
            np.savetxt(f, zip(u, v, w, nu, sub_real, sub_imag, std_sub_real, std_sub_imag), header=FILEHEADER)
        del sub_real
        del sub_imag
        '''
        # Save data to files.
        # The regular data.
        FILEHEADER = 'Baseline: %d-%d\nEntries: %d\nu [m], v [m], w [m], frequency [GHz], real, imag, std(real), std(imag)' % (ant1, ant2, nu.shape[0])
        with open('visibilities/visibilities_corr_%.2d.txt'%(corr,), 'ab') as f:
            np.savetxt(f, zip(u, v, w, nu, data_real, data_imag, std_real, std_imag), header=FILEHEADER)
        # Write back errors and weights to the SIGMA and WEIGHT columns of the MS file.
        weights = sigma ** -2
        ct.taql('UPDATE $msfile SET SIGMA[$corr]=$sigma WHERE ANTENNA1=$ant1 AND ANTENNA2=$ant2')
        ct.taql('UPDATE $msfile SET WEIGHT[$corr]=$weights WHERE (ANTENNA1=$ant1 AND ANTENNA2=$ant2)')
        progress += 1
    print '100%\n'

print '[CVE] Closing MS file...'
msfile.close()
print '[CVE] Finished'
