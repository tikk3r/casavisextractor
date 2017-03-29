# casavisextractor
The CASA Visibility Extractor (CVE for short) will extract the real and imaginary parts of the visibilities with their corresponding errors from a Measurement Set to be used for further analysis.

Functionality
=============
- The following data is extracted: u, v, w, channel frequency, real part, imaginary part, real error, imaginary error in the *visibilities* file. The *visibilities_subtracted* files contain the same data, except the visibilities were subtracted from each other as 2-1, 4-3, 6-5 etc. before calculating the uncertainties.
- Visibilities are stored in txt files with separate files per correlation. The data themselves are sorted by frequency on a per baseline basis.
- Histograms of the visibilities and the errors can be plotted.

Requirements
============
The following software is required to run the scripts.

extract_taql.py
---------------
- casacore
- python-casacore
- numpy

plot_sigmas.py
--------------
- matplotlib
