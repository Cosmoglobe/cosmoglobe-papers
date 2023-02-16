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

width = 14.5

# Load data
cgK30 = np.loadtxt('Clhat_K30_CG.txt')
wmapK30 = np.loadtxt('Clhat_K30_WM.txt')
cg30Ka = np.loadtxt('Clhat_30Ka_CG.txt')
wmap30Ka = np.loadtxt('Clhat_30Ka_WM.txt')
cgQ = np.loadtxt('Clhat_dQ_CG.txt')
wmapQ = np.loadtxt('Clhat_dQ_WM.txt')
cgV = np.loadtxt('Clhat_dV_CG.txt')
wmapV = np.loadtxt('Clhat_dV_WM.txt')
cgW = np.loadtxt('Clhat_dW_CG.txt')
wmapW = np.loadtxt('Clhat_dW_WM.txt')

#lcdm = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')

vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(0.8*cm2inch(width), cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)


###############
#   K-30
###############

# EE
ax1 = plt.subplot2grid((5, 2), (0, 0))

plt.plot(wmapK30[1,2:30],   label='WMAP',  linewidth=1, color='red')
plt.plot(cgK30[1,2:30],   label='CG',  linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
ax1.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax1.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$K$-30", fontsize=12, ha='right')
plt.title(r"EE", fontsize=12)

# BB
ax2 = plt.subplot2grid((5, 2), (0, 1))

plt.plot(wmapK30[2,2:30],   linewidth=1, color='red')
plt.plot(cgK30[2,2:30],     linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax2.yaxis.labelpad = 10*width/17.; ax2.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax2.set_xscale("log")
ax2.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$K$-30", fontsize=12, ha='right')
#plt.text(16,0.9,r"$K$-30, BB", fontsize=12)
plt.title(r"BB", fontsize=12)


###############
#   30-Ka
###############

# EE
ax3 = plt.subplot2grid((5, 2), (1, 0))

plt.plot(wmap30Ka[1,2:30],   label='WMAP',  linewidth=1, color='red')
plt.plot(cg30Ka[1,2:30],   label='CG',  linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax3.set_xscale("log")
ax3.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax3.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$Ka$-30", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, EE", fontsize=12)


# BB
ax4 = plt.subplot2grid((5, 2), (1, 1))

plt.plot(wmap30Ka[2,2:30],   linewidth=1, color='red')
plt.plot(cg30Ka[2,2:30],     linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax4.yaxis.labelpad = 10*width/17.; ax4.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax4.set_xscale("log")
ax4.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$Ka$-30", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, BB", fontsize=12)



###############
#   Delta Q
###############

# EE
ax5 = plt.subplot2grid((5, 2), (2, 0))

plt.plot(wmapQ[1,2:30],   label='WMAP',  linewidth=1, color='red')
plt.plot(cgQ[1,2:30],   label='CG',  linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
plt.ylabel(r"Power spectrum, $C_{\ell}\,[\mu\mathrm{K}^2]$"); 
ax5.yaxis.labelpad = 10*width/17.; ax5.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax5.set_xscale("log")
ax5.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax5.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$Q1-Q2$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, EE", fontsize=12)


# BB
ax6 = plt.subplot2grid((5, 2), (2, 1))

plt.plot(wmapQ[2,2:30],   linewidth=1, color='red')
plt.plot(cgQ[2,2:30],     linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax6.yaxis.labelpad = 10*width/17.; ax6.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax6.set_xscale("log")
ax6.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$Q1-Q2$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, BB", fontsize=12)



###############
#   Delta V
###############

# EE
ax7 = plt.subplot2grid((5, 2), (3, 0))

plt.plot(wmapV[1,2:30],   label='WMAP',  linewidth=1, color='red')
plt.plot(cgV[1,2:30],   label='CG',  linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax7.yaxis.labelpad = 10*width/17.; ax7.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax7.set_xscale("log")
ax7.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax7.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$V1-V2$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, EE", fontsize=12)


# BB
ax8 = plt.subplot2grid((5, 2), (3, 1))

plt.plot(wmapV[2,2:30],   linewidth=1, color='red')
plt.plot(cgV[2,2:30],     linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
#plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax8.yaxis.labelpad = 10*width/17.; ax8.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax8.set_xscale("log")
ax8.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
plt.setp( ax8.get_xticklabels(), visible=False)
plt.setp( ax8.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$V1-V2$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, BB", fontsize=12)



###############
#   Delta W
###############

# EE
ax9 = plt.subplot2grid((5, 2), (4, 0))

plt.plot(wmapW[1,2:30],   label='WMAP',  linewidth=1, color='red')
plt.plot(cgW[1,2:30],   label='CG',  linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax9.yaxis.labelpad = 10*width/17.; ax9.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax9.set_xscale("log")
ax9.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

plt.xticks([0,4,8,12], [r"$2$", r"$6$", r"$10$", r"$14$"])
#plt.xticks([0.01,0,010,100,1000,10000,100000], [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$"])
#plt.yticks([0.01,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax9.get_xticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$(W1+W2)-(W3+W4)$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, EE", fontsize=12)


# BB
ax10 = plt.subplot2grid((5, 2), (4, 1))

plt.plot(wmapW[2,2:30],   linewidth=1, color='red')
plt.plot(cgW[2,2:30],     linewidth=1, color='black')

# x axis
#plt.hlines(0, 0, 3300)

# legend
#leg = plt.legend(frameon=True, loc=3, fontsize=10)
# remove box around legend
#leg.get_frame().set_edgecolor("white")
#leg.get_frame().set_alpha(0)

# labels
plt.xlabel(r"Multipole moment, $\ell$");
#plt.ylabel(r"$C_{\ell}^{K-30}\,[\mu\mathrm{K}^2]$"); 
ax10.yaxis.labelpad = 10*width/17.; ax10.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax10.set_xscale("log")
ax10.set_yscale("log")
plt.ylim([0.003, 30]); plt.xlim([0, 16]);

plt.xticks([0,4,8,12, 16], [r"$2$", r"$6$", r"$10$", r"$14$", 18])
#plt.yticks([0.0,1,100], [r"$10^{-2}$", r"$10^0$", r"$10^2$"])
#plt.setp( ax10.get_xticklabels(), visible=False)
plt.setp( ax10.get_yticklabels(), visible=False)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

plt.text(15,8,r"$(W1+W2)-(W3+W4)$", fontsize=12, ha='right')
#plt.text(16,0.9,r"30-$Ka$, BB", fontsize=12)



# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax3.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax9.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")    


#ax2.errorbar(ls, binned1, yerr=rms1, fmt='.', color='red')

# save to pdf with right bounding box
plt.savefig("cls_cg_WMAP_lowl.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

# Make table

