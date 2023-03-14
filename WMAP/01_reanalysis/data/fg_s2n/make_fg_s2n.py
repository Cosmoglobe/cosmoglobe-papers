# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 8.8

# Load data
cmbT   = np.loadtxt('cmb_s2n_T.dat')
synchT = np.loadtxt('synch_s2n_T.dat')
ffT    = np.loadtxt('ff_s2n_T.dat')
ameT   = np.loadtxt('ame_s2n_T.dat')
dustT  = np.loadtxt('dust_s2n_T.dat')
cmbP   = np.loadtxt('cmb_s2n_P.dat')
synchP = np.loadtxt('synch_s2n_P.dat')
dustP  = np.loadtxt('dust_s2n_P.dat')

x = np.arange(0,2)


vmin = -110
vmax =  160

# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 2.5*cm2inch(width)))


fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



# CMB T
ax1 = plt.subplot2grid((8, 1), (0, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax1.errorbar(cmbT[:,0]+0., cmbT[:,1], yerr=cmbT[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax1.errorbar([7], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red', label='Channel with highest S/N')
ax1.errorbar([1], [3.2], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red', label='Channel with highest S/N')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 4]);
ax1.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,2,3], [r"$1$", "2",r"$3$"])
plt.setp( ax1.get_xticklabels(), visible=True)
plt.text(8.3,3,r"CMB, $T$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.minorticks_off()
plt.text(1.2,3.,r"Channel with highest S/N", fontsize=8, ha='left')


# synch T
ax2 = plt.subplot2grid((8, 1), (1, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax2.errorbar(synchT[:,0]+0., synchT[:,1], yerr=synchT[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax2.errorbar([1], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 200]);
ax2.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,10,100], [r"$1$", r"$10$", r"100"])
plt.setp( ax2.get_xticklabels(), visible=True)
plt.text(8.3,1.3,r"Synchrotron, $T$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


# freefree T
ax3 = plt.subplot2grid((8, 1), (2, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax3.errorbar(ffT[:,0]+0., ffT[:,1], yerr=ffT[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax3.errorbar([2], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 30]);
ax3.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1, 10], [r"$1$", r"$10$"])
plt.setp( ax3.get_xticklabels(), visible=True)
plt.text(8.3,1.,r"Free-free, $T$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


# AME T
ax4 = plt.subplot2grid((8, 1), (3, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax4.errorbar(ameT[:,0]+0., ameT[:,1], yerr=ameT[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax4.errorbar([1], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 200]);
ax4.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,10,100], [r"$1$", r"$10$", r"100"])
plt.setp( ax4.get_xticklabels(), visible=True)
plt.text(8.3,1.2,r"AME, $T$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


# Dust T
ax5 = plt.subplot2grid((8, 1), (4, 0))
plt.locator_params(nbins=5)
plt.ylabel(r"Relative signal-to-noise ratio, $\frac{(S/N)_\mathrm{max}}{S/N}$", ha='center');
#plt.xlabel(r"Radiometer"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax5.errorbar(dustT[:,0]+0., dustT[:,1], yerr=dustT[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax5.errorbar([7], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 30]);
ax5.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,10], [r"$1$", r"$10$"])
plt.setp( ax5.get_xticklabels(), visible=True)
plt.text(8.3,17,r"Thermal dust, $T$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


# CMB P
ax6 = plt.subplot2grid((8, 1), (5, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax6.errorbar(cmbP[:,0]+0., cmbP[:,1], yerr=cmbP[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax6.errorbar([7], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 4]);
ax6.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,2,3], [r"$1$", "2", r"$3$"])
plt.setp( ax6.get_xticklabels(), visible=True)
plt.text(8.3,17,r"CMB, $P$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.minorticks_off()
ax6.set_facecolor((0.9,0.9,0.9))

# Synch P
ax7 = plt.subplot2grid((8, 1), (6, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
#plt.xlabel(r"Radiometer"); 
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax7.errorbar(synchP[:,0]+0., synchP[:,1], yerr=synchP[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax7.errorbar([1], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 200]);
ax7.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"30", r"44", r"70", r"$K$", r"$Ka$", r"$Q$", r"$V$", r"$W$"])
plt.yticks([1,10,100], [r"$1$", r"$10$", r"100"])
plt.setp( ax7.get_xticklabels(), visible=True)
plt.text(8.3,1.2,r"Synchrotron, $P$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax7.set_facecolor((0.9,0.9,0.9))


# Dust P
ax8 = plt.subplot2grid((8, 1), (7, 0))
plt.locator_params(nbins=5)
#plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
plt.xlabel(r"Frequency channel"); 
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
ax8.errorbar(dustP[:,0]+0., dustP[:,1], yerr=dustP[:,2], fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='black')
ax8.errorbar([7], [1], yerr=[0], fmt='.', ms=6, capsize=1, capthick=1, elinewidth=1, color='red')
plt.grid(False, which="major", axis="both")
plt.ylim([0.8, 30]);
ax8.set_yscale("log")
plt.xlim([0.5, 8.5]);
plt.xticks([1,2,3,4,5,6,7,8], [r"$K$", r"$30$", r"$Ka$", r"$Q$", r"$44$", r"$V$", r"$70$", r"$W$"])
plt.yticks([1,10], [r"$1$", r"$10$"])
plt.setp( ax8.get_xticklabels(), visible=True)
plt.text(8.3,17,r"Thermal dust, $P$", fontsize=8, ha='right')
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax8.set_facecolor((0.9,0.9,0.9))



# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# set vertical y axis ticklables
for ticklabel in ax8.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")    

#for ticklabel in ax1.xaxis.get_ticklabels():
#    ticklabel.set_rotation(45.)    








# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("fg_s2n_v2.pdf", bbox_inches='tight',
    bbox_extra_artists=[],pad_inches=0.03, dpi=100)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

