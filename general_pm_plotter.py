#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/ICN2/spower/projects/hueso_sys/res/res/GB558_NLINF4RIB_4_LEADS_1_to_0__1e-06/AC320_GB558_l2_200_ANDL_2.0_ANDW_3.00_yc_20/POTDIS_c_0.0001_d_0.020_xi_15.0/xymap.dat"
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

#plt.imshow(pmi, vmin=pm.min(), vmax=pm.max(), origin='lower', cmap=cm.jet, 
#           extent=[x.min(), x.max(), y.min(), y.max()])
plt.imshow(pmi, vmin=-0.04, vmax=0.04, origin='lower', cmap=cm.jet, 
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plt.show(plot3)
plot3.savefig(outputname, dpi=300)


