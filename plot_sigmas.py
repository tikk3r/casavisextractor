from matplotlib.pyplot import figure, show
import argparse
import numpy as np

correlations = 2
for i in range(correlations):
    # Load and plot the normal visibilities and their statistics.
    real, imag, sigma_real, sigma_imag = np.loadtxt('./visibilities/visibilities_corr_%.2d.txt'%i, usecols=(4,5,6,7), unpack=True)
    fig = figure('Visibilities Correlation %d/%d'%(i+1,correlations), figsize=(12,8), dpi=100)
    fig.suptitle('Visibilities Correlation %d/%d'%(i+1,correlations), fontweight='bold')
    ax = fig.add_subplot(221)
    ax.hist(real, bins=100)
    ax.set_title('Real', fontweight='bold')

    ax2 = fig.add_subplot(222)
    ax2.hist(sigma_real, bins=100)
    ax2.set_title('Real $\\sigma$', fontweight='bold')

    ax3 = fig.add_subplot(223)
    ax3.hist(imag, bins=100)
    ax3.set_title('Imag', fontweight='bold')

    ax4 = fig.add_subplot(224)
    ax4.hist(sigma_imag, bins=100)
    ax4.set_title('Imag $\\sigma$', fontweight='bold')
    fig.savefig('./visibilities/visibilities_corr_%d.pdf'%i, dpi=400)

    # Load and plot the subtracted visibilities and their statistics.
    try:
        real, imag, sigma_real, sigma_imag = np.loadtxt('./visibilities/visibilities_subtracted_corr_%.2d.txt'%i, usecols=(4,5,6,7), unpack=True)
        fig = figure('Subtracted Visibilities Correlation %d/%d'%(i+1,correlations), figsize=(12,8), dpi=100)
        fig.suptitle('Subtracted Visibilities Correlation %d/%d'%(i+1,correlations), fontweight='bold')
        ax = fig.add_subplot(221)
        ax.hist(real, bins=100)
        ax.set_title('Real', fontweight='bold')

        ax2 = fig.add_subplot(222)
        ax2.hist(sigma_real, bins=100)
        ax2.set_title('Real $\\sigma$', fontweight='bold')

        ax3 = fig.add_subplot(223)
        ax3.hist(imag, bins=100)
        ax3.set_title('Imag', fontweight='bold')

        ax4 = fig.add_subplot(224)
        ax4.hist(sigma_imag, bins=100)
        ax4.set_title('Imag $\\sigma$', fontweight='bold')
        fig.savefig('./visibilities/visibilities_subtracted_corr_%d.pdf'%i, dpi=400)

    except IOError:
        print 'Subtracted visibilities not found, skipping...'
#show()
try:
    print 'Attempting to create pdf files with pdfunite...'
    import subprocess
    subprocess.call('pdfunite ./visibilities/visibilities_corr_*.pdf ./visibilities/visibilities.pdf', shell=True)
except Exception as e:
    print e
    print 'PDF creation failed.'
else:
    print 'PDF created.'

try:
    print 'Attempting to create pdf files with pdfunite...'
    import subprocess
    subprocess.call('pdfunite ./visibilities/visibilities_subtracted_corr_*.pdf ./visibilities/visibilities_subtracted.pdf', shell=True)
except Exception as e:
    print e
    print 'PDF creation failed.'
else:
    print 'PDF created.'
