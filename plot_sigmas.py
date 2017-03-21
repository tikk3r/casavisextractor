from matplotlib.pyplot import figure, show
import numpy as np

real = np.loadtxt('./visibilities/visibilities.txt', usecols=(4,), unpack=True)
imag = np.loadtxt('./visibilities/visibilities.txt', usecols=(5,), unpack=True)
sigma_real = np.loadtxt('./visibilities/visibilities.txt', usecols=(6,), unpack=True)
sigma_imag = np.loadtxt('./visibilities/visibilities.txt', usecols=(7,), unpack=True)
fig = figure(figsize=(12,8), dpi=100)
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
fig.savefig('visibilities.png', dpi=400)
show()
