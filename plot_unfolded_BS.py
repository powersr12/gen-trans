#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/spow/projects/gen_trans/res/ANTIDOTS_RIBBON_1e-06/ZZ56_rect_lat_L_7_12_circ_dot_R_2.0_4x1_xyf_0.0_rf_0.0/Eloop_+0.00_to_+0.50_Bfixed_+200.000_run.conf00.wbands"
cleanbands="/home/spow/projects/gen_trans/res/CLEAN_RIBBON_1e-06/ZZ56_clean_l2_1/Eloop_+0.00_to_+0.50_Bfixed_+200.000_run.conf00.bands"
outputname=filename+"_plot.png"

x, band, weight = np.loadtxt(filename, unpack=True)

cleandata = np.loadtxt(cleanbands)
cleanys = cleandata[:,1:]
cleanx = cleandata[:,0]

print cleanx.shape
print cleandata.shape
print cleanys.shape

for bandc in cleanys.T:
  plt.plot(cleanx, bandc, color="red", linewidth=1)
  


rgba_colors = np.zeros((x.shape[0],4))
## for red the first column needs to be one
rgba_colors[:,0] = 0.0
rgba_colors[:,1] = 0.0
rgba_colors[:,2] = 0.0


## the fourth column needs to be your alphas
rgba_colors[:, 3] = weight.clip(0, 10)

plt.scatter(x, band, color=rgba_colors, facecolors=rgba_colors, marker='o', s=30, lw=0)
plt.ylim(-0.1, 1.0)

plt.show()



