#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/powersr/projects/gen-trans/res/BUBBLES_RIBBON_1e-06/AC240_rect_lat_L_72_120_membrane_bub_R_40.6_H_12.2_1x1_xyf_0.0_rf_0.0_zf_0.0/POTDIS_c_0.0200_d_0.020_xi_10.0/abs_pot/Eloop_+0.01_to_+0.10_Bfixed_+0.000_abs5.conf00.disprof"
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

plt.imshow(pmi, vmin=pm.min(), vmax=pm.max(), origin='lower', cmap=cm.jet, 
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
plot3.savefig(outputname, dpi=300)
plt.show(plot3)


