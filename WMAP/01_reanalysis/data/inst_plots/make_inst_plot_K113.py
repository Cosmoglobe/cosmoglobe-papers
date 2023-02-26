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
slopeK = np.loadtxt('baseslope_CG_023-WMAP_K_v1.dat')
gainK = np.loadtxt('gain_CG_023-WMAP_K_v1.dat')
maskK = np.loadtxt('mask_CG_023-WMAP_K_v1.dat')
baseK = np.loadtxt('baseline_CG_023-WMAP_K_v1.dat')
sigmaK = np.loadtxt('sigma0_CG_023-WMAP_K_v1.dat')
fkneeK = np.loadtxt('fknee_CG_023-WMAP_K_v1.dat')
alphaK = np.loadtxt('alpha_CG_023-WMAP_K_v1.dat')
chisqK = np.loadtxt('chisq_CG_023-WMAP_K_v1.dat')

sigmaK[:,1:5] = sigmaK[:,1:5] / gainK[:,1:5] * N.sqrt(1.536/12) / N.sqrt(2.)# K 
fkneeK[:,1:5] = 1000*fkneeK[:,1:5]


for i in range(4):
    inds = np.where(maskK[:,i+1] == 1)
    baseK[:,i+1] = baseK[:,i+1] - N.mean(baseK[inds,i+1])
    
inds = np.where(maskK == 0)
slopeK[inds] = np.nan
alphaK[inds] = np.nan
baseK[inds] = np.nan
chisqK[inds] = np.nan
sigmaK[inds] = np.nan
alphaK[inds] = np.nan
fkneeK[inds] = np.nan
gainK[inds] = np.nan



wmapgain = np.loadtxt('regressed_gains.txt')


vmin = -110
vmax =  160


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 2*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)

wmap = N.zeros((20,2))
wmap[0,:] = [0.66, 0.40]

gsfc = N.zeros((20,2))
gsfc[0,:] = [0.72, 6.13]

mjd_wmap = [52130, 52477]
mjd_gsfc = [52130, 55412]


###############
#   gain
###############

mjd_gain = gainK[0,0] + np.arange(1,3280+1)/3280. * (gainK[-1,0]-gainK[0,0])
wmap2 = wmapgain[0,:]
inds = (np.abs(wmap2) < 0.01)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax1 = plt.subplot2grid((7, 1), (0, 0))
#plt.plot(mjd_gain,1./np.abs(wmap2), linewidth=0.5, color='red')
plt.plot(gainK[:,0],gainK[:,1], linewidth=1, color='black')
data = np.loadtxt('K113_g0.txt')
plt.plot(data[::73,0], abs(data[::73,1]), linewidth=0.5, color='red')
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,1.25,r"K113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax1.yaxis.labelpad = 10*width/17.

###############
#   baseline
###############

ax2 = plt.subplot2grid((7, 1), (1, 0))
plt.plot(baseK[:,0],baseK[:,1], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,70,r"K113", fontsize=10)
plt.ylabel(r"$b - \left<b\right>$ [du]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-100, 100]);
plt.yticks([-60,0,60], [r"$-60$", r"$0$", r"$60$"])
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])


###############
#   slope
###############

ax3 = plt.subplot2grid((7, 1), (2, 0))
plt.plot(slopeK[:,0],slopeK[:,1], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,2.8,r"K113", fontsize=10)
plt.ylabel(r"$b_1$ [du]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-1, 1]);
plt.yticks([-0.5,0,0.5], [r"$-0.5$", r"$0$", r"$0.5$"])


###############
#   sigma0
###############

ax4 = plt.subplot2grid((7, 1), (3, 0))
plt.plot(sigmaK[:,0],sigmaK[:,1], linewidth=1, color='black', label='CG')
plt.plot(mjd_wmap,[wmap[0,0],wmap[0,0]], linewidth=1, color='red', linestyle=':', label='WMAP')
plt.plot(mjd_gsfc,[gsfc[0,0],gsfc[0,0]], linewidth=1, color='orange', linestyle=':', label='GSFC')
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.ylim([0.65, 0.75]);
#plt.text(52200,0.835,r"K113", fontsize=10)
plt.ylabel(r"$\sigma_{0}$ [mK\,$\mathrm{s}^{\frac{1}{2}}$]");
ax1.yaxis.labelpad = 10*width/17.
plt.yticks([0.67,0.70,0.73], [r"$0.67$", r"$0.70$", r"$0.73$"])


###############
#   fknee
###############

ax5 = plt.subplot2grid((7, 1), (4, 0))
plt.plot(fkneeK[:,0],fkneeK[:,1], linewidth=1, color='black', label='CG')
plt.plot(mjd_wmap,[wmap[0,1],wmap[0,1]], linewidth=1, color='red', linestyle=':', label='WMAP')
plt.plot(mjd_gsfc,[gsfc[0,1],gsfc[0,1]], linewidth=1, color='orange', linestyle=':', label='GSFC')
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,4.5,r"K113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([0, 7]);
plt.yticks([2,4,6], [r"$2$", r"$4$", r"$6$"])


###############
#   alpha
###############

ax6 = plt.subplot2grid((7, 1), (5, 0))
plt.plot(alphaK[:,0],alphaK[:,1], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,-0.6,r"K113", fontsize=10)
plt.ylabel(r"$\alpha$");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-1.4, -0.6]);
plt.yticks([-1.3,-1.0,-0.7], [r"$-1.3$", r"$-1.0$", r"$-0.7$"])

###############
#   chisq
###############

ax7 = plt.subplot2grid((7, 1), (6, 0))
plt.plot(chisqK[:,0],chisqK[:,1], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=True)
plt.setp( ax7.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,11.2,r"K113", fontsize=10)
plt.ylabel(r"$\chi^2$ [$\sigma$]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-9, 2]);
plt.yticks([-6,-3,0], [r"$-6$", r"$-3$", r"$0$"])
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])

plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
# labels
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax2.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# save to pdf with right bounding box
plt.savefig("../../figures/instpar_CG_K113_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

