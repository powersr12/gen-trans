#!/usr/bin/python
 
import matplotlib
import matplotlib.pyplot as plt
print matplotlib.__version__
from matplotlib.pyplot import cm
import numpy as np
from matplotlib.colors import LogNorm
import scipy.interpolate
 
#basefolder="/home/ICN2/spower/projects/gen-trans/res/"
##basefolder="/home/spow//GEN_TRANS_RESULTS/"
#projfolder="ANTIDOTS_RIBBON_1e-06/ZZ288_rect_lat_L_36_64_circ_dot_R_5.0_4x8_xyf_0.0_rf_0.0/"
#projname="E_+0.40_B_+161.000_longbuf.conf00"
 
basefolder="/home/powersr/projects/"
projfolder="bubbles/"
projname="E_+0.0090_B_+0.000_abs5.conf00"
 
 
 
basefilename=basefolder+projfolder+projname
dosname=basefilename+".ldos_red.dat"
odosname=dosname+"__.png"
filename=basefilename+"_currentsum.dat"
outputname=filename+"_lr.png"
##leads = {'l0', 'l2', 'l3', 'l1', 'l4', 'l5', 'multi'}
leads = {'l0'}
#leads = {'multi'}
#leads = {'l0', 'l1', 'l2', 'l3', 'multi'}
runzoom = 0
 
#print dosname
 
##load coordinates and DOS
x, y, ldos = np.loadtxt(dosname, delimiter=' ', unpack=True)
#ldos = np.loadtxt(dosname, delimiter=' ', usecols = (2,), unpack=True)
x, y = x*0.246, y*0.246
#print jx, ldos
#print ldos.min()
 
##regular grid
xi, yi = np.linspace(x.min(), x.max(), 800), np.linspace(y.min(), y.max(), 800 )
xi, yi = np.meshgrid(xi, yi)
 
#arrow grid
xim, yim = np.linspace(x.min(), x.max(), 30), np.linspace(y.min(), y.max(), 30)
xim, yim = np.meshgrid(xim, yim)
 
## Interpolate DOS and map, output
plot3 = plt.figure()
ldosi = scipy.interpolate.griddata((x, y), ldos, (xi, yi), method='linear')
 
plt.imshow(ldosi, vmin=0.0, vmax=0.8*ldosi.max(), origin='lower', cmap=cm.inferno, aspect='equal', extent=[xi.min(), xi.max(), yi.min(), yi.max()]            )
##plt.colorbar()
plot3.savefig(odosname, dpi=300, bbox_inches='tight')
plt.show(plot3)
 
 
#for lead in leads:
  #filename=basefilename+".cmaps_"+lead+"_red.dat"
  #jx, jy =np.loadtxt(filename, delimiter=' ', usecols = (2,3), unpack=True)
  #curmag = np.sqrt(jx**2 + jy**2)
  #maxcur = max(curmag)
  ##plot2 = plt.figure()
  #zi = scipy.interpolate.griddata((x, y), curmag, (xi, yi), method='linear')
  ##plt.imshow(zi, vmin=curmag.min(), vmax=0.9*curmag.max(), origin='lower', cmap=cm.hot,
           ##extent=[x.min(), x.max(), y.min(), y.max()])
  ##plt.show(plot2)
  ##outputname=filename+".png"
  ##plot2.savefig(outputname, dpi=300)
   
   
   
  ##plot4 = plt.figure()
  #jxi = scipy.interpolate.griddata((x, y), jx, (xi, yi), method='linear')
  #jyi = scipy.interpolate.griddata((x, y), jy, (xi, yi), method='linear')
  #curmagi=np.sqrt(jxi**2 + jyi**2)
  ##plt.quiver(xi, yi, jxi, jyi,        # data
           ##curmagi,                   # colour the arrows based on this array
           ##cmap=cm.seismic,     # colour map
           ##headlength=3)        # length of the arrows
 
  ###plt.colorbar()                  # adds the colour bar
 
  ##plt.title('Quive Plot, Dynamic Colours')
  ##plt.show(plot4)                 # display the plot
  #lw = 2.0*curmagi/curmag.max()
  ##plot5 = plt.figure()
  ##plt.imshow(zi, vmin=curmag.min(), vmax=0.25*curmag.max(), origin='lower', cmap=cm.coolwarm,
           ##extent=[x.min(), x.max(), y.min(), y.max()])
  ##plt.streamplot(xi, yi, jxi, jyi,          # data
               ##color='k',
               ##density=7.0,
               ##linewidth=lw,         # line thickness
               ###arrowstyle='->',     # arrow style
               ##arrowsize=3.0)       # arrow size
 
  ###plt.colorbar()                      # add colour bar on the right
 
  ###plt.title(projfolder+projname+" "+lead)
  ##plt.xlim(x.min(), x.max())
  ##plt.ylim(y.min(), y.max())
  ##plt.show(plot5)   
  ##outputname=filename+".png"
  ##plot5.savefig(outputname, dpi=300, bbox_inches='tight')
   
   
  ##plot6=plt.figure()
  ##plt.imshow(zi, vmin=curmag.min(), vmax=0.25*curmag.max(), origin='lower', cmap=cm.coolwarm,
           ##extent=[x.min(), x.max(), y.min(), y.max()])
  ##jxim = scipy.interpolate.griddata((x, y), jx, (xim, yim), method='linear')
  ##jyim = scipy.interpolate.griddata((x, y), jy, (xim, yim), method='linear')
  ##curmagim=np.sqrt(jxi**2 + jyi**2)
  ##plt.quiver(xim, yim, jxim, jyim,        # data
           ##headlength=2, headaxislength=1.5, width=0.005, headwidth=2, scale=7.5, minshaft=2)        # length of the arrows
  ##plt.show(plot6)   
  ##plot6.savefig(filename+"alt.png", dpi=300, bbox_inches='tight')
   
   
   
   
### Interpolate current mag
##plot2 = plt.figure()
##zi = scipy.interpolate.griddata((x, y), curmag, (xi, yi), method='linear')
 
##plt.imshow(zi, vmin=curmag.min(), vmax=curmag.max(), origin='lower', cmap=cm.hot,
           ##extent=[x.min(), x.max(), y.min(), y.max()])
###plt.colorbar()
##plt.show(plot2)
 
##plot2.savefig(outputname, dpi=300)
 
 
 
 
 
##left vs right moving channel comparison
##xi, yi = np.linspace(x.min(), x.max(), 2000), np.linspace(y.min(), y.max(), 2000)
##xi, yi = np.meshgrid(xi, yi)
 
### Interpolate
##plot2 = plt.figure()
##zi = scipy.interpolate.griddata((x, y), jx1, (xi, yi), method='linear')
 
##plt.imshow(zi, vmin=jx1.min(), vmax=jx1.max(), origin='lower', cmap=cm.seismic,
           ##extent=[x.min(), x.max(), y.min(), y.max()])
##plt.colorbar()
##plt.show(plot2)
 
##plot2.savefig(outputname, dpi=300)