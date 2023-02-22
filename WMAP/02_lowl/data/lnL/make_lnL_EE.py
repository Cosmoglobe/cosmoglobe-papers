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

width = 12.0

# Load data
wmap2 = np.loadtxt('lnL_EE_wmap-1.txt')
wmap3 = np.loadtxt('lnL_EE_wmap-2.txt')
wmap4 = np.loadtxt('lnL_EE_wmap-3.txt')
wmap5 = np.loadtxt('lnL_EE_wmap-4.txt')
wmap6 = np.loadtxt('lnL_EE_wmap-5.txt')
wmap7 = np.loadtxt('lnL_EE_wmap-6.txt')
wmap8 = np.loadtxt('lnL_EE_wmap-7.txt')
wmap9 = np.loadtxt('lnL_EE_wmap-8.txt')
wmap10 = np.loadtxt('lnL_EE_wmap-9.txt')

reprod2 = np.loadtxt('lnL_EE_reprod-1.txt')
reprod3 = np.loadtxt('lnL_EE_reprod-2.txt')
reprod4 = np.loadtxt('lnL_EE_reprod-3.txt')
reprod5 = np.loadtxt('lnL_EE_reprod-4.txt')
reprod6 = np.loadtxt('lnL_EE_reprod-5.txt')
reprod7 = np.loadtxt('lnL_EE_reprod-6.txt')
reprod8 = np.loadtxt('lnL_EE_reprod-7.txt')
reprod9 = np.loadtxt('lnL_EE_reprod-8.txt')
reprod10 = np.loadtxt('lnL_EE_reprod-9.txt')


tempcorr2 = np.loadtxt('lnL_EE_tempcorr-1.txt')
tempcorr3 = np.loadtxt('lnL_EE_tempcorr-2.txt')
tempcorr4 = np.loadtxt('lnL_EE_tempcorr-3.txt')
tempcorr5 = np.loadtxt('lnL_EE_tempcorr-4.txt')
tempcorr6 = np.loadtxt('lnL_EE_tempcorr-5.txt')
tempcorr7 = np.loadtxt('lnL_EE_tempcorr-6.txt')
tempcorr8 = np.loadtxt('lnL_EE_tempcorr-7.txt')
tempcorr9 = np.loadtxt('lnL_EE_tempcorr-8.txt')
tempcorr10 = np.loadtxt('lnL_EE_tempcorr-9.txt')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(1.3*cm2inch(width), cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



###############
#   ell = 2
###############

ax1 = plt.subplot2grid((3, 3), (0, 0))
plt.plot(wmap2[:,0],np.exp(wmap2[:,1]-np.max(wmap2[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod2[:,0],np.exp(reprod2[:,1]-np.max(reprod2[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr2[:,0],np.exp(tempcorr2[:,1]-np.max(tempcorr2[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax1.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=2$", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax1.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"MJD");
plt.ylabel(r"$\mathcal{L}/\mathcal{L}_\mathrm{max}$");
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


###############
#   ell = 3
###############

ax2 = plt.subplot2grid((3, 3), (0, 1))
plt.plot(wmap3[:,0],np.exp(wmap3[:,1]-np.max(wmap3[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod3[:,0],np.exp(reprod3[:,1]-np.max(reprod3[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr3[:,0],np.exp(tempcorr3[:,1]-np.max(tempcorr3[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax2.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=3$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax2.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


###############
#   ell = 4
###############

ax3 = plt.subplot2grid((3, 3), (0, 2))
plt.plot(wmap4[:,0],np.exp(wmap4[:,1]-np.max(wmap4[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='Official')
plt.plot(reprod4[:,0],np.exp(reprod4[:,1]-np.max(reprod4[:,1])), linewidth=1, color='red', label='Uncorrected')
plt.plot(tempcorr4[:,0],np.exp(tempcorr4[:,1]-np.max(tempcorr4[:,1])), linewidth=1, color='blue', label='Corrected')
plt.grid(False, which="major", axis="both")
plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax3.set_xscale("log")
plt.text(3,0.7,r"$\ell=4$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax3.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# legend
leg = plt.legend(frameon=True, loc=1, fontsize=6)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)


###############
#   ell = 5
###############

ax4 = plt.subplot2grid((3, 3), (1, 0))
plt.plot(wmap5[:,0],np.exp(wmap5[:,1]-np.max(wmap5[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod5[:,0],np.exp(reprod5[:,1]-np.max(reprod5[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr5[:,0],np.exp(tempcorr5[:,1]-np.max(tempcorr5[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax4.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=5$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax4.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
plt.ylabel(r"$\mathcal{L}/\mathcal{L}_\mathrm{max}$");
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

###############
#   ell = 6
###############

ax5 = plt.subplot2grid((3, 3), (1, 1))
plt.plot(wmap6[:,0],np.exp(wmap6[:,1]-np.max(wmap6[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod6[:,0],np.exp(reprod6[:,1]-np.max(reprod6[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr6[:,0],np.exp(tempcorr6[:,1]-np.max(tempcorr6[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax5.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=6$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax5.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


###############
#   ell = 7
###############

ax6 = plt.subplot2grid((3, 3), (1, 2))
plt.plot(wmap7[:,0],np.exp(wmap7[:,1]-np.max(wmap7[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod7[:,0],np.exp(reprod7[:,1]-np.max(reprod7[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr7[:,0],np.exp(tempcorr7[:,1]-np.max(tempcorr7[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax6.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=7$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax6.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


###############
#   ell = 8
###############

ax7 = plt.subplot2grid((3, 3), (2, 0))
plt.plot(wmap8[:,0],np.exp(wmap8[:,1]-np.max(wmap8[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod8[:,0],np.exp(reprod8[:,1]-np.max(reprod8[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr8[:,0],np.exp(tempcorr8[:,1]-np.max(tempcorr8[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=True)
plt.setp( ax7.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax7.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=8$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax7.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0,0.5,1], [r"$0.0$", r"$0.5$", r"$1.0$"])
plt.xticks([0.01,0.1,1], [r"$10^{-2}$", r"$10^{-1}$", r"$10^0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
plt.ylabel(r"$\mathcal{L}/\mathcal{L}_\mathrm{max}$");
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
plt.xlabel(r"Power spectrum, $D_{\ell}^{EE}$ [$\mu\mathrm{K}^2$]");
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

###############
#   ell = 9
###############

ax8 = plt.subplot2grid((3, 3), (2, 1))
plt.plot(wmap9[:,0],np.exp(wmap9[:,1]-np.max(wmap9[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod9[:,0],np.exp(reprod9[:,1]-np.max(reprod9[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr9[:,0],np.exp(tempcorr9[:,1]-np.max(tempcorr9[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax8.get_xticklabels(), visible=True)
plt.setp( ax8.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax8.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=9$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax8.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
plt.yticks([0,0.5,1], [r"$0.0$", r"$0.5$", r"$1.0$"])
plt.xticks([0.01,0.1,1], [r"$10^{-2}$", r"$10^{-1}$", r"$10^0$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
#plt.xlabel(r"MJD");
#ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
plt.xlabel(r"Power spectrum, $D_{\ell}^{EE}$ [$\mu\mathrm{K}^2$]");
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels


###############
#   ell = 10
###############

ax9 = plt.subplot2grid((3, 3), (2, 2))
plt.plot(wmap10[:,0],np.exp(wmap10[:,1]-np.max(wmap10[:,1])), linewidth=0.5, alpha=0.5, color='black', linestyle='--', label='CG')
plt.plot(reprod10[:,0],np.exp(reprod10[:,1]-np.max(reprod10[:,1])), linewidth=1, color='red', label='CG')
plt.plot(tempcorr10[:,0],np.exp(tempcorr10[:,1]-np.max(tempcorr10[:,1])), linewidth=1, color='blue', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax9.get_xticklabels(), visible=True)
plt.setp( ax9.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
ax9.set_xscale("log")
plt.text(2.5,0.95,r"$\ell=10$", fontsize=10)
#plt.text(52200,-0.6,r"K111", fontsize=10)
#plt.ylabel(r"$\alpha$");
ax9.yaxis.labelpad = 10*width/17.
plt.xlim([0.01, 10]); plt.ylim([0., 1.1]);
plt.yticks([0.5,1], [r"$0.5$", r"$1.0$"])
plt.yticks([0,0.5,1], [r"$0.0$", r"$0.5$", r"$1.0$"])
plt.xticks([0.01,0.1,1,10], [r"$10^{-2}$", r"$10^{-1}$", r"$10^0$", r"$10^1$"])

## legend
#leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"Power spectrum, $D_{\ell}^{EE}$ [$\mu\mathrm{K}^2$]");
ax9.yaxis.labelpad = 10*width/17.; ax9.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels




# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")



# save to pdf with right bounding box
plt.savefig("lnL_EE_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

