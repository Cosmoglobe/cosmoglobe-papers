# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp


width = 8.8

# Load data
wmap = np.loadtxt('dust_tempfit_WMAP.dat')
lfi = np.loadtxt('dust_tempfit_LFI.dat')
spdust = np.loadtxt('spdust2_cnm.dat')

nu = spdust[:,0]
spdust[:,0] = spdust[:,0] / 30. * 21.
spdust[:,1] = spdust[:,1] / (nu*nu)
spdust[:,1] = spdust[:,1] / N.max(spdust[:,1])*90.

expmod = 47 * N.exp(-3.57 * (nu-23.)/23.)
powlaw = 51 * (nu/23.)**(-4.9)

#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)

####################
#   SED
####################

# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 0.7*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)


ax1 = plt.subplot2grid((1, 1), (0, 0))

#plt.locator_params(nbins=5)


# x axis
#plt.hlines(0, 0, 3300)

# labels
plt.xlabel(r"Frequency, $\mathrm{GHz}$", fontsize=12);
plt.ylabel(r"Template amplitude", fontsize=12); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


plt.plot(nu, spdust[:,1], color='green', label=r"SPDust, $\nu_{\mathrm{p}} = 21$\,GHz", linewidth=0.5)
#plt.plot(nu, powlaw, color='blue', label=r'Power-law, $\beta = -4.9$', linewidth=0.5)
plt.plot(nu, expmod, color='red', label=r'Exponential, $\beta = -3.57$', linewidth=1)

ax1.errorbar(wmap[:,0], wmap[:,1], yerr=wmap[:,2], fmt='.', ms=5, capsize=1, capthick=1, elinewidth=1, color='black', label='WMAP')
ax1.errorbar(lfi[:,0], lfi[:,1], yerr=lfi[:,2], fmt='.', ms=5, capsize=1, capthick=1, elinewidth=1, color='darkgray', label='LFI')




# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
ax1.set_xscale("log")
ax1.set_yscale("log")
plt.ylim([0.011, 100]);
plt.xlim([20.1, 80]);

#plt.xticks([30,50,70], [r"30", r"50", r"70"])
#plt.yticks([-0.75,-0.5,-0.25,0,0.25,0.5,0.75], [r"$-0.75$", r"$-0.50$", r"$-0.25$", "0.00", "0.25", "0.50", "0.75"])
plt.setp( ax1.get_xticklabels(), visible=False)


# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

#for ticklabel in ax1.xaxis.get_ticklabels():
#    ticklabel.set_rotation(45.)    

# legend
leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)





# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
plt.savefig("tempfit_AME_single_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

