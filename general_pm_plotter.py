#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/stpower/projects/gated_triangles/dtures/SUBLDOTS_CUSTOM_4_LEADS_1_to_0__1e-06/ZZ324_rect_lat_L_54_96_SUBA_1.00x0.100_SUBB_1.00x0.100_triZZ_dot_up_R_0.0_3x3_xyf_0.0_rf_0.0/POTDIS_c_0.0010_d_0.010_xi_8.0/map.dat"
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


