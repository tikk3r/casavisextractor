from matplotlib.pyplot import figure, show
import numpy as np

sigma_real = np.loadtxt('./visibilities/visibilities.txt', usecols=(6,), unpack=True)
sigma_imag = np.loadtxt('./visibilities/visibilities.txt', usecols=(7,), unpack=True)
fig = figure()
ax = fig.add_subplot(211)
ax.hist(sigma_real)
ax.set_title('Real', fontweight='bold')

ax2 = fig.add_subplot(212)
ax2.hist(sigma_imag)
ax2.set_title('Imag', fontweight='bold')
show()
