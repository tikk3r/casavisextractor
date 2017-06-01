from __future__ import division
import argparse
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

desc = '''Process the visibilities in an MS file in order to determine better estimates of the sigmas and weights for them.\nCASA assigns the same value to each baseline as 1  / sqrt(BW*T) whereas this script will estimate the sigmas for each baseline by looking at the scatter in the visibilities through time. The weights are then calculated as 1 / sigma**2.'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('--subtract', action='store_true', dest='subtract', help='calculate with the subtracted visibilities instead')
parser.add_argument('filename', action='store')
args = parser.parse_args()
SUBTRACT = args.subtract
filename = args.filename
if SUBTRACT:
    print '[MSERR] Using the subtracted visibilties to calculate sigmas and weights.'
else:
    print '[MSERR] Using the regular visibilties to calculate sigmas and weights.'

'''
Read in MS file and create table object.
'''
print '[MSERR] Reading in file...'
msfile = ct.table(filename, readonly=False)

'''
Calculate the number of baselines present in the measurement set.
'''
print '[MSERR] Forming baselines...'
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
print '[MSERR] Extracting visibilities per baseline...'
try:
    os.mkdir('visibilities')
except:
    pass

for corr in range(correlations):
    print 'Processing correlation %d/%d:' % (corr+1, correlations[0])
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
        baseline = ct.taql('SELECT DATA,WEIGHT FROM $msfile WHERE ANTENNA1=$ant1 AND ANTENNA2=$ant2')
        data = baseline.getcol('DATA')
        weight = baseline.getcol('WEIGHT')
        # Select the correlation, specified by the last index.
        # data has the shape (timestamps, channels, correlations).
        subdata = data[...,corr]
        subweight = weight[...,corr]
        # subdata now has the shape (timestamps, channels)
        subdata = subdata.flatten()
        # Take the real and imaginary parts of the data and handle them separately.
        data_real = subdata.real
        data_imag = subdata.imag
        # Save data to file.
        std = 1 / np.sqrt(subweight)
        FILEHEADER = 'Baseline: %d-%d\nreal, imag, std' % (ant1, ant2)
        with open('visibilities/sigmas_corr_%.2d.txt'%(corr,), 'ab') as f:
            np.savetxt(f, zip(data_real, data_imag, std), header=FILEHEADER)
        progress += 1
    print '100%\n'

print '[MSERR] Closing MS file...'
msfile.close()
print '[MSERR] Finished'
