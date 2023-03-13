# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
import cosmoglobe

import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 19.0

ind2chan = {}
ind2name = {}

for i in range(40):
    if i < 4:
        name = 'K'
        ind2chan[i] = '023-WMAP_K'
    elif i < 8:
        name = 'Ka'
        ind2chan[i] = '030-WMAP_Ka'
    elif i < 12:
        name = 'Q1'
        ind2chan[i] = '040-WMAP_Q1'
    elif i < 16:
        name = 'Q2'
        ind2chan[i] = '040-WMAP_Q2'
    elif i < 20:
        name = 'V1'
        ind2chan[i] = '060-WMAP_V1'
    elif i < 24:
        name = 'V2'
        ind2chan[i] = '060-WMAP_V2'
    elif i < 28:
        name = 'W1'
        ind2chan[i] = '090-WMAP_W1'
    elif i < 32:
        name = 'W2'
        ind2chan[i] = '090-WMAP_W2'
    elif i < 36:
        name = 'W3'
        ind2chan[i] = '090-WMAP_W3'
    else:
        name = 'W4'
        ind2chan[i] = '090-WMAP_W4'
    if i % 4 == 0:
        ind2name[i] = f'{name}13'
    if i % 4 == 1:
        ind2name[i] = f'{name}14'
    if i % 4 == 2:
        ind2name[i] = f'{name}23'
    if i % 4 == 3:
        ind2name[i] = f'{name}24'

# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 1.35*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)

mjd_wmap = [52130, 52477]
mjd_gsfc = [52130, 55414]

plot_texts = [r"K113", r"K114", r"K123", r"K124", r"Ka113", r"Ka114", r"Ka123", r"Ka124", r"Q113", r"Q114", r"Q123", r"Q124", r"Q213", r"Q214", r"Q223", r"Q224", r"V113", r"V114", r"V123", r"V124", r"V213", r"V214", r"V223", r"V224", r"W113", r"W114", r"W123", r"W124", r"W213", r"W214", r"W223", r"W224", r"W313", r"W314", r"W323", r"W324", r"W413", r"W414", r"W423", r"W424"]
plot_text_coords = [[52200, -0.6]] * 40  # all channels

chain_fname1 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_a_230206/chain_c0001.h5'
chain_fname2 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203/chain_c0001.h5'
chain1 = cosmoglobe.h5.chain.Chain(chain_fname1)
chain2 = cosmoglobe.h5.chain.Chain(chain_fname2)
axs = []
for i in range(40):
    if i % 4 == 0:
        currax = plt.subplot2grid((10, 4), (i//4, i % 4))
        firstax = currax
    else:
        currax = plt.subplot2grid((10, 4), (i//4, i%4), sharey=firstax)
    axs.append(currax)
    channel = ind2chan[i]
    samples1 = chain1.get(f'tod/{channel}/xi_n')[:, 2, :, :]
    samples2 = chain2.get(f'tod/{channel}/xi_n')[:, 2, :, :]
    tot_samples = np.concatenate((samples1[50:],  samples2[50:]), axis=0)
    mean = np.mean(tot_samples, axis=0)
    std = np.std(tot_samples, axis=0)
    accept = chain1.get(f'tod/{channel}/accept')[-1, i % 4, :].astype(bool)
    mjd = chain1.get(f'tod/{channel}/MJD')[-1, :]
    plt.errorbar(mjd[accept], mean[i%4, :][accept], yerr=std[i%4][accept],
        linewidth=0.5, color='black', zorder=1, label = 'CG', rasterized=True)


    plt.grid(False, which="major", axis="both")
    if i >= 36:
        plt.setp( currax.get_xticklabels(), visible=True)
    else:
        plt.setp( currax.get_xticklabels(), visible=False)
    if i % 4 == 0:
        plt.setp( currax.get_yticklabels(), visible=True)
    else:
        plt.setp( currax.get_yticklabels(), visible=False)
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.text(plot_text_coords[i][0],plot_text_coords[i][1], plot_texts[i], fontsize=10)
    plt.ylim([-1.6, -0.4]);
    if i % 4 == 0:
        plt.ylabel(r"$\alpha$");
        currax.yaxis.labelpad = 10*width/17.
        plt.yticks([-1.4,-1.0,-0.6], [r"$-1.4$", r"$-1.0$", r"$-0.6$"])
    if i >= 36:
        if i % 4 == 0:
            plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
        else:
            plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])
        ## labels
        plt.xlabel(r"MJD");
    if i == 0:
        leg = plt.legend(frameon=True, loc=1, fontsize=8)
        # remove box around legend
        leg.get_frame().set_edgecolor("white")
        leg.get_frame().set_alpha(0)

for i in range(0, len(axs), 4):
    for ticklabel in axs[i].yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

# save to pdf with right bounding box
plt.savefig("../../figures/instpar_CG_alpha_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
