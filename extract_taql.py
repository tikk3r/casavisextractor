from __future__ import division
import casacore.tables as ct
import numpy as np

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
Read in MS file and extract the required columns into numpy arrays.
'''
print '[CVE] Reading in file...'
msfile = ct.table('TESTVIS.ms', readonly=False)

print '[CVE] Forming baselines...'
ANTENNA1 = msfile.getcol('ANTENNA1')
ANTENNA2 = msfile.getcol('ANTENNA2')

ANTENNAS = set(ANTENNA1).union(set(ANTENNA2))
baselines = form_baselines(ANTENNAS)

# Extract the data for each baseline.
print '[CVE] Extracting visibilities per baseline...'
progress = 0; end = len(baselines); printed = False
f = open('visibilities_test.txt', 'wb'); f.close()
for i,j in baselines:
    p = int(progress / end * 100)
    if ((progress % 10) == 0):
        print 'Progress: %d%%' % (p)
    bl = ct.taql('SELECT DATA FROM $msfile WHERE ANTENNA1=$i AND ANTENNA2=$j')
    data = bl.getcol('DATA')
    # Average polarizations together.
    data_avg_pol = np.average(data, axis=2)
    data_avg_pol_flat = data_avg_pol.flatten()
    data_real = data_avg_pol_flat.real
    data_imag = data_avg_pol_flat.imag

    n = len(data_avg_pol_flat)
    antenna1 = np.zeros(n); antenna1.fill(i)
    antenna2 = np.zeros(n); antenna2.fill(j)
    # Save data to file.
    with open('visibilities_test.txt', 'ab') as f:
        np.savetxt(f, zip(antenna1, antenna2, data_real, data_imag))
    progress += 1

print '[CVE] Calculating statistics...'
real, imag = np.loadtxt('visibilities_test.txt', usecols=(2,3), unpack=True)
mean_re = real.mean(); mean_im = imag.mean()
std_re = real.std(); std_im = imag.std()

sig = [std_re, std_im]

print 'Re mean: ', real.mean()
print 'Re std: ', real.std()
print 'Im mean: ', imag.mean()
print 'Im std: ', imag.std()

print '[CVE] Writing back to MS file...'
ct.taql('UPDATE $msfile SET SIGMA=$sig')

print '[CVE] Finished.'
