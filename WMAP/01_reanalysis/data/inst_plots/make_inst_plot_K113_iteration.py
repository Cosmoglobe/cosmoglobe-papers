# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as np
import scipy.stats
import healpy as hp
import cosmoglobe as cg

import plotly.colors as pcol
import matplotlib as mpl

DIR = '/mn/stornext/d5/data/duncanwa/WMAP'


cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 8.8

# Load data
chain = cg.Chain(f'{DIR}/chains_CG_b_230203/chain_c0001.h5')

chisq = chain.get('tod/023-WMAP_K/chisq')
xi_n = chain.get('tod/023-WMAP_K/xi_n')
baseline = chain.get('tod/023-WMAP_K/baseline')
gain = chain.get('tod/023-WMAP_K/gain')
x_im = chain.get('tod/023-WMAP_K/x_im')

pid = 533

sigma = xi_n[2:,0,0,pid]
fknee = xi_n[2:,1,0,pid]*1e3
alpha = xi_n[2:,2,0,pid]
base = baseline[2:,0,0,pid]
slope = baseline[2:,1,0,pid]
gain = gain[2:,0,pid]
chisq = chisq[2:,0,pid]

print(x_im.shape)


base -= base.mean()
slope -= slope.mean()

sigma = sigma / gain * np.sqrt(1.536/12) / np.sqrt(2.)

samps = np.arange(len(sigma)) + 2



vmin = -110
vmax =  160


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 2*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



###############
#   gain
###############


ax1 = plt.subplot2grid((7, 1), (0, 0))
plt.plot(samps,gain, linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,1.25,r"K113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax1.yaxis.labelpad = 10*width/17.

###############
#   x_im1
###############

ax2 = plt.subplot2grid((7, 1), (1, 0))
plt.plot(samps,x_im[2:,0], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,70,r"K113", fontsize=10)
plt.ylabel(r"$x_\mathrm{im,1}$")
ax1.yaxis.labelpad = 10*width/17.
#plt.ylim([-100, 100]);
#plt.yticks([-60,0,60], [r"$-60$", r"$0$", r"$60$"])
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])


###############
#   x_im2
###############

ax3 = plt.subplot2grid((7, 1), (2, 0))
plt.plot(samps,x_im[2:,1], linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,2.8,r"K113", fontsize=10)
plt.ylabel(r"$x_\mathrm{im,2}$")
ax1.yaxis.labelpad = 10*width/17.
#plt.ylim([-1, 1]);
#plt.yticks([-0.5,0,0.5], [r"$-0.5$", r"$0$", r"$0.5$"])


###############
#   sigma0
###############

ax4 = plt.subplot2grid((7, 1), (3, 0))
plt.plot(samps,sigma, linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.ylim([0.7015, 0.70615]); #plt.text(52200,0.835,r"K113", fontsize=10)
plt.ylabel(r"$\sigma_{0}$ [mK\,$\mathrm{s}^{\frac{1}{2}}$]"); ax1.yaxis.labelpad = 10*width/17.
#plt.yticks([0.67,0.70,0.73], [r"$0.67$", r"$0.70$", r"$0.73$"])


###############
#   fknee
###############

ax5 = plt.subplot2grid((7, 1), (4, 0))
plt.plot(samps,fknee, linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,4.5,r"K113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([0.6, 1.1]);
plt.yticks([0.75, 1.00], [r"$0.75$", r"$1.00$"])


###############
#   alpha
###############

ax6 = plt.subplot2grid((7, 1), (5, 0))
plt.plot(samps,alpha, linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,-0.6,r"K113", fontsize=10)
plt.ylabel(r"$\alpha$");
ax1.yaxis.labelpad = 10*width/17.
#plt.ylim([-1.4, -0.6]);
#plt.yticks([-1.3,-1.0,-0.7], [r"$-1.3$", r"$-1.0$", r"$-0.7$"])

###############
#   chisq
###############

ax7 = plt.subplot2grid((7, 1), (6, 0))
plt.plot(samps,chisq, linewidth=1, color='black', label='CG')
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=True)
plt.setp( ax7.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,11.2,r"K113", fontsize=10)
plt.ylabel(r"$\chi^2$ [$\sigma$]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-3, -0.5]);
#plt.yticks([-6,-3,0], [r"$-6$", r"$-3$", r"$0$"])
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
# labels
plt.xlabel(r"Gibbs Iteration");
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
plt.savefig('../../figures/instpar_CG_K113_samples_v1.pdf', bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
