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
filename = 'TESTVIS.ms'
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

progress = 0; end = len(baselines); printed = False
for i,j in baselines:
    p = int(progress / end * 100)
    if ((p % 10) == 0) and not printed:
        print '%d%%' % (p),
        printed = True
    elif printed and (p%10 == 9):
        printed = False
    sys.stdout.flush()
    # Select wanted columns (possibly slightly redundant step).
    baseline = ct.taql('SELECT UVW,DATA FROM $msfile WHERE ANTENNA1=$i AND ANTENNA2=$j')
    data = baseline.getcol('DATA')
    uvw = baseline.getcol('UVW')
    # Select the spectral window keyword from the main table.
    spws = ct.taql('SELECT CHAN_FREQ FROM %s::SPECTRAL_WINDOW'%(filename))
    # Frequencies corresponding to each channel.
    frequencies = spws.getcol('CHAN_FREQ')
    nu = frequencies.flatten() * 1e-9
    # u,v,w coordinates for each baseline in meters.
    u, v, w = uvw[...,0], uvw[...,1], uvw[...,2]
    # Determine the number of correlations present.
    correlations = data.shape[-1]
    # Store the sigmas to write back to the MS file here.
    sigmas = []
    # Loop over every correlation.
    for corr in range(correlations):
        # Select the correlation, specified by the last index.
        subdata = data[...,corr]
        subdata = subdata.flatten()
        data_real = subdata.real
        data_imag = subdata.imag
        # Calculate standard deviation.
        stdr = data_real.std()
        stdi = data_imag.std()
        std_real = np.zeros(len(data_real)); std_real.fill(stdr)
        std_imag = np.zeros(len(data_imag)); std_imag.fill(stdi)
        # The MS file only has one sigma per correlation, so take the largest.
        std = np.asarray([r if r > i else i for r,i in zip(std_real, std_imag)])
        sigmas.append(std)

        # Save data to file.
        FILEHEADER = 'Baseline: %d-%d\nEntries: %d\nu [m], v [m], w [m], frequency [GHz], real, imag, std(real), std(imag)' % (i, j, nu.shape[0])
        #with open('visibilities/baseline%.2d-%.2d_corr%.2d.txt'%(i,j,corr), 'ab') as f:
        with open('visibilities/visibilities.txt', 'ab') as f:
            np.savetxt(f, zip(u, v, w, nu, data_real, data_imag, std_real, std_imag), header=FILEHEADER)
    # Write back errors and weights to the SIGMA and WEIGHT columns of the MS file.
    sigmas = np.asarray(sigmas)
    weights = np.asarray(sigmas) ** -2
    ct.taql('UPDATE $msfile SET SIGMA=$sigmas')
    ct.taql('UPDATE $msfile SET WEIGHTS=$weights')
    progress += 1
print '100%'

print '[CVE] Closing MS file...'
msfile.close()
print '[CVE] Finished'
