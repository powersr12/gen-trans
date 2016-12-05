#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/ICN2/spower/projects/gen-trans/res/CLEAN_rnltest_4_LEADS_1_to_0__1e-06/ZZ200_clean_l2_1000/POTDIS_c_0.01_d_0.200_xi_3.0/Bloop_-20.00_to_+20.00_Efixed_+0.000_NEW.conf00.disprof"
outputname=filename+"_plot.png"


#print dosname

##load coordinates and DOS
x, y, pm = np.loadtxt(filename, unpack=True)
#ldos = np.loadtxt(dosname, delimiter=' ', usecols = (2,), unpack=True)

absmax = max ( pm.max(), abs(pm.min() ))

##regular grid
xi, yi = np.linspace(x.min(), x.max(), 500), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)

## Interpolate DOS and map, output
plot3 = plt.figure()
pmi = scipy.interpolate.griddata((x, y), pm , (xi, yi), method='linear')

plt.imshow(pmi, vmin=-absmax, vmax=absmax, origin='lower', cmap=cm.PuOr, 
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.show(plot3)
plot3.savefig(outputname, dpi=300)


