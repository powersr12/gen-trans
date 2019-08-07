#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="../res/CLEAN_RIBBON_1e-06/ZZ80_clean_l2_200/POTDIS_c_0.0010_d_2.000_xi_2.0/.Eloop_-0.40_to_+0.40_Bfixed_+0.000_run.conf00.disprof"
outputname=filename+"_plot.png"


#print dosname

##load coordinates and DOS
x, y, pm = np.loadtxt(filename, unpack=True)
#ldos = np.loadtxt(dosname, delimiter=' ', usecols = (2,), unpack=True)

absmax = max ( pm.max(), abs(pm.min() ))

##regular grid
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 800)
xi, yi = np.meshgrid(xi, yi)

## Interpolate DOS and map, output
plot3 = plt.figure()
pmi = scipy.interpolate.griddata((x, y), pm , (xi, yi), method='linear')

plt.imshow(pmi, vmin=pm.min(), vmax=pm.max(), origin='lower', cmap=cm.RdBu, 
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plot3.savefig(outputname, dpi=300)
plt.show(plot3)


