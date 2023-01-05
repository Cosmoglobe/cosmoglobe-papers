# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 8.8

# Load data
bp30 = np.loadtxt('bp_RIMO_R3_030.dat')
bp44 = np.loadtxt('bp_RIMO_R3_044.dat')
bp70 = np.loadtxt('bp_RIMO_R3_070.dat')

#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(1.4*cm2inch(width), cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



ax1 = plt.subplot2grid((1, 1), (0, 0))

#plt.plot(bp30[:,0], bp30[:,0],  linewidth=1, color='red', alpha=1, )

plt.locator_params(nbins=5)


# x axis
#plt.hlines(0, 0, 3300)

# labels
plt.xlabel(r"Frequency, $\mathrm{GHz}$", fontsize=12);
plt.ylabel(r"Normalized bandpass profile, $\tau$", fontsize=12); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

plt.plot(bp30[:,0], bp30[:,1], color='red', label='30 GHz', linewidth=0.5)
plt.plot(bp44[:,0], bp44[:,1], color='green', label='44 GHz', linewidth=0.5)
plt.plot(bp70[:,0], bp70[:,1], color='blue', label='70 GHz', linewidth=0.5)

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
plt.ylim([0, 0.4]);
#plt.xlim([0, 28]);

plt.xticks([30,50,70], [r"30", r"50", r"70"])
#plt.yticks([-0.75,-0.5,-0.25,0,0.25,0.5,0.75], [r"$-0.75$", r"$-0.50$", r"$-0.25$", "0.00", "0.25", "0.50", "0.75"])
#plt.setp( ax1.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

#for ticklabel in ax1.xaxis.get_ticklabels():
#    ticklabel.set_rotation(45.)    

# legend
leg = plt.legend(frameon=True, loc=1, fontsize=12)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)







# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("bp_LFI_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

