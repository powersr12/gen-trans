#!/usr/bin/python


import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate

filename="/home/ICN2/spower/projects/vqhe/res/SUBLMOIRE_RIBBON_1e-06/ZZ118_L2_59_rescaled_off_0.0_0.0/Eloop_-0.30_to_+0.30_Bfixed_+0.000_cont.conf00.wbands"
cleanbands="/home/ICN2/spower/projects/vqhe/res/SUBLMOIRE_RIBBON_1e-06/ZZ118_L2_59_rescaled_off_0.0_0.0/Eloop_-0.30_to_+0.30_Bfixed_+0.000_cont.conf00.bands"
outputname=filename+"_plot.png"




x, band, weight = np.loadtxt(filename, unpack=True)

#cleandata = np.loadtxt(cleanbands)
#cleanys = cleandata[:,1:]
#cleanx = cleandata[:,0]

#print cleanx.shape
#print cleandata.shape
#print cleanys.shape

#for bandc in cleanys.T:
#  plt.plot(cleanx, bandc, color="red", linewidth=1)
  


rgba_colors = np.zeros((x.shape[0],4))
## for red the first column needs to be one
rgba_colors[:,0] = 0.0
rgba_colors[:,1] = 0.0
rgba_colors[:,2] = 0.0


## the fourth column needs to be your alphas
rgba_colors[:, 3] = weight.clip(0, 10)

plt.scatter(x, band, color=rgba_colors, facecolors=rgba_colors, marker='o', s=30, lw=0)
plt.ylim(-0.2, 0.2)

plt.show()



