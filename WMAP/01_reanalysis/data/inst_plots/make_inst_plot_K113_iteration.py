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
chain_a = cg.Chain(f'{DIR}/chains_CG_a_230206/chain_c0001.h5')
chain_b = cg.Chain(f'{DIR}/chains_CG_b_230203/chain_c0001.h5')

chisq_a = chain_a.get('tod/023-WMAP_K/chisq')
xi_n_a = chain_a.get('tod/023-WMAP_K/xi_n')
baseline_a = chain_a.get('tod/023-WMAP_K/baseline')
gain_a = chain_a.get('tod/023-WMAP_K/gain')
x_im_a = chain_a.get('tod/023-WMAP_K/x_im')

chisq_b = chain_b.get('tod/023-WMAP_K/chisq')
xi_n_b = chain_b.get('tod/023-WMAP_K/xi_n')
baseline_b = chain_b.get('tod/023-WMAP_K/baseline')
gain_b = chain_b.get('tod/023-WMAP_K/gain')
x_im_b = chain_b.get('tod/023-WMAP_K/x_im')



pid = 533
pid = 50

mjds = chain_a.get('tod/023-WMAP_K/MJD')[-1]
print(mjds[pid], mjds[pid+1])

sigma_a = xi_n_a[2:,0,0,pid]
fknee_a = xi_n_a[2:,1,0,pid]*1e3
alpha_a = xi_n_a[2:,2,0,pid]
base_a = baseline_a[2:,0,0,pid]
slope_a = baseline_a[2:,1,0,pid]
gain_a = gain_a[2:,0,pid]
chisq_a = chisq_a[2:,0,pid]

sigma_b = xi_n_b[2:,0,0,pid]
fknee_b = xi_n_b[2:,1,0,pid]*1e3
alpha_b = xi_n_b[2:,2,0,pid]
base_b = baseline_b[2:,0,0,pid]
slope_b = baseline_b[2:,1,0,pid]
gain_b = gain_b[2:,0,pid]
chisq_b = chisq_b[2:,0,pid]



base_a -= base_a.mean()
slope_a -= slope_a.mean()
base_b -= base_b.mean()
slope_b -= slope_b.mean()

#sigma_a = sigma_a / gain_a * np.sqrt(1.536/12) / np.sqrt(2.)
#sigma_b = sigma_b / gain_b * np.sqrt(1.536/12) / np.sqrt(2.)

samps_a = np.arange(len(sigma_a)) + 2
samps_b = np.arange(len(sigma_b)) + 2



vmin = -110
vmax =  160


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 2*cm2inch(width)*5/7))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



###############
#   gain
###############


ax1 = plt.subplot2grid((5, 1), (0, 0))
plt.plot(samps_a,gain_a, linewidth=1, color='C0', rasterized=True)
plt.plot(samps_b,gain_b, linewidth=1, color='C1', rasterized=True)
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,1.25,r"K113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([0.949, 0.9549])
plt.yticks([0.95, 0.953], ['$0.950$', '$0.953$'])

###############
#   x_im1
###############

#ax2 = plt.subplot2grid((7, 1), (1, 0))
#plt.plot(samps_a,x_im_a[2:,0], linewidth=1, color='C0', label='CG')
#plt.plot(samps_b,x_im_b[2:,0], linewidth=1, color='C1', label='CG')
#plt.grid(False, which="major", axis="both")
#plt.setp( ax2.get_xticklabels(), visible=False)
#plt.setp( ax2.get_yticklabels(), visible=True)
#plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
##plt.text(52200,70,r"K113", fontsize=10)
#plt.ylabel(r"$x_\mathrm{im,1}$")
#ax1.yaxis.labelpad = 10*width/17.
#plt.ylim([-5e-4, 1.5e-3]);
#plt.subplots_adjust(left=0, right=1, top=1, bottom=0)


###############
#   x_im2
###############

#ax3 = plt.subplot2grid((7, 1), (2, 0))
#plt.plot(samps_a,x_im_a[2:,1], linewidth=1, color='C0', label='CG')
#plt.plot(samps_b,x_im_b[2:,1], linewidth=1, color='C1', label='CG')
#plt.grid(False, which="major", axis="both")
#plt.setp( ax3.get_xticklabels(), visible=False)
#plt.setp( ax3.get_yticklabels(), visible=True)
##plt.text(52200,2.8,r"K113", fontsize=10)
#plt.ylabel(r"$x_\mathrm{im,2}$")
#ax1.yaxis.labelpad = 10*width/17.
##plt.yticks([-0.5,0,0.5], [r"$-0.5$", r"$0$", r"$0.5$"])
 
 
###############
#   sigma0
###############

ax4 = plt.subplot2grid((5, 1), (1, 0))
plt.plot(samps_a,sigma_a, linewidth=1, color='C0', label='CG', rasterized=True)
plt.plot(samps_b,sigma_b, linewidth=1, color='C1', label='CG', rasterized=True)
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.ylim([0.7015, 0.70615]); #plt.text(52200,0.835,r"K113", fontsize=10)
plt.ylabel(r"$\sigma_{0}$ [du]"); ax1.yaxis.labelpad = 10*width/17.
plt.ylim([2.65347, 2.655])
plt.ylim(  [2.65347, 2.6545])
plt.yticks([2.654,   2.6543])
plt.yticks([2.6537,   2.6542], ['$2.6537$', '$2.6542$'])


###############
#   fknee
###############

ax5 = plt.subplot2grid((5, 1), (2, 0))
plt.plot(samps_a,fknee_a, linewidth=1, color='C0', label='CG', rasterized=True)
plt.plot(samps_b,fknee_b, linewidth=1, color='C1', label='CG', rasterized=True)
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,4.5,r"K113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([0.65, 1.1]);
plt.yticks([0.75, 1.00], [r"$0.75$", r"$1.00$"])


###############
#   alpha
###############

ax6 = plt.subplot2grid((5, 1), (3, 0))
plt.plot(samps_a,alpha_a, linewidth=1, color='C0', label='CG', rasterized=True)
plt.plot(samps_b,alpha_b, linewidth=1, color='C1', label='CG', rasterized=True)
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,-0.6,r"K113", fontsize=10)
plt.ylabel(r"$\alpha$");
ax1.yaxis.labelpad = 10*width/17.
plt.ylim([-1.29, -0.9]);
#plt.yticks([-1.3,-1.0,-0.7], [r"$-1.3$", r"$-1.0$", r"$-0.7$"])

###############
#   chisq
###############

ax7 = plt.subplot2grid((5, 1), (4, 0))
plt.plot(samps_a,chisq_a, linewidth=1, color='C0', label='CG', rasterized=True)
plt.plot(samps_b,chisq_b, linewidth=1, color='C1', label='CG', rasterized=True)
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=True)
plt.setp( ax7.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
#plt.text(52200,11.2,r"K113", fontsize=10)
plt.ylabel(r"$\chi^2$ [$\sigma$]");
ax1.yaxis.labelpad = 10*width/17.
#plt.ylim([-3, -0.5]);
plt.ylim([-8.7,-6.8])
#plt.yticks([-6,-3,0], [r"$-6$", r"$-3$", r"$0$"])
#plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])

#plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
# labels
plt.xlabel(r"Gibbs Iteration");


# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

#for ticklabel in ax2.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")
#
#for ticklabel in ax3.yaxis.get_ticklabels():
#    ticklabel.set_rotation("vertical")

for ticklabel in ax4.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax6.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax7.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

# save to pdf with right bounding box
plt.savefig('../../figures/instpar_CG_K113_samples_v1.pdf', bbox_inches='tight',
    bbox_extra_artists=[],pad_inches=0.03, dpi=100)
