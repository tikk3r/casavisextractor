# Measurement Set Error Estimator
The Measurement Set Error Estimator (MaSER for short) will extract the real and imaginary parts of the visibilities from a Measurement Set and calculates their corresponding errors to be used for further analysis. Errors are estimated for the real and imaginary parts separately, but as there is only one value stored for each correlation the largest of the two is written back to the MS file.

Note that this script does not work on concatenated datasets, because of the antenna numbers increasing beyond the actual number of antennas after using concat.

This script modifies three columns in the MS file: the _SIGMA_, _WEIGHT_ and _WEIGHT_SPECTRUM_ columns. The sigmas are estimated from the visibilities and the weights are determined as 1 / sigma\*\*2.

Usage
=====
To estimate errors and weights in the visibilities run

    python extract_taql.py observation.ms

Optional flags that can be passed along are:

    --subtract: operate on pairwise-subtracted visibilities instead of the raw visibilities. That is instead of using visibilities 1, 2, 3, 4, 5, 6, ... to determine the errors, 1-2, 2-3, 3-4, ... etc. are used. This also writes out the subtracted visibilities.
    --backup: make a backup of the input measurement set before operating on it. This will be stored as observation.ms.BACKUP.
    --hdf5: experimental feature to store the data in the HDF5 file format instead of plain text. This results in much faster I/O.

Functionality
=============
- The following data is extracted: u, v, w, channel frequency, real part, imaginary part, real error, imaginary error in the *visibilities* file. The *visibilities_subtracted* files contain the same data, except the visibilities were subtracted from each other as 1-2, 2-3, 3-4 etc. before calculating the uncertainties.
- Visibilities are stored in txt files with separate files per correlation. The data themselves are sorted by frequency on a per baseline basis.

Requirements
============
The following software is required to run the scripts.

extract_taql.py
---------------
- [casacore](https://github.com/casacore/casacore)
- [python-casacore](https://github.com/casacore/python-casacore)
- numpy
