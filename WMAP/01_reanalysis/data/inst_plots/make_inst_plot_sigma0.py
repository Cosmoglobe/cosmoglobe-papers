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
wmap = {}
gsfc = {}

for i in range(40):
    if i < 4:
        name = 'K'
        ind2chan[i] = '023-WMAP_K'
        if i in (0, 1):
            wmap[i] = [0.66, 0.40]
            gsfc[i] = [0.72, 6.13]
        else:
            wmap[i] = [0.75, 0.51]
            gsfc[i] = [0.87, 5.37]
    elif i < 8:
        name = 'Ka'
        ind2chan[i] = '030-WMAP_Ka'
        if i in (4, 5):
            wmap[i] = [0.71, 0.71]
            gsfc[i] = [0.75, 1.66]
        elif i in (6, 7):
            wmap[i] = [0.72, 0.32]
            gsfc[i] = [0.77, 1.29]
    elif i < 12:
        name = 'Q1'
        ind2chan[i] = '040-WMAP_Q1'
        if i in (8, 9):
            wmap[i] = [0.92, 1.09]
            gsfc[i] = [0.99, 3.21]
        elif i in (10, 11):
            wmap[i] = [1.02, 0.35]
            gsfc[i] = [0.95, 3.13]
    elif i < 16:
        name = 'Q2'
        ind2chan[i] = '040-WMAP_Q2'
        if i in (12, 13):
            wmap[i] = [0.85, 5.76]
            gsfc[i] = [0.89, 1.92]
        elif i in (14, 15):
            wmap[i] = [0.99, 8.62]
            gsfc[i] = [1.04, 4.61]
    elif i < 20:
        name = 'V1'
        ind2chan[i] = '060-WMAP_V1'
        if i in (16, 17):
            wmap[i] = [1.22, 0.09]
            gsfc[i] = [1.25, 2.56]
        elif i in (18, 19):
            wmap[i] = [1.11, 1.41]
            gsfc[i] = [1.07, 4.49]
    elif i < 24:
        name = 'V2'
        ind2chan[i] = '060-WMAP_V2'
        if i in (20, 21):
            wmap[i] = [0.97, 0.88]
            gsfc[i] = [1.01, 2.43]
        elif i in (22, 23):
            wmap[i] = [1.10, 8.35]
            gsfc[i] = [1.13, 3.06]
    elif i < 28:
        name = 'W1'
        ind2chan[i] = '090-WMAP_W1'
        if i in (24, 25):
            wmap[i] = [1.35, 7.88]
            gsfc[i] = [1.18, 16.2]
        elif i in (26, 27):
            wmap[i] = [1.61, 0.66]
            gsfc[i] = [1.41, 15.1]
    elif i < 32:
        name = 'W2'
        ind2chan[i] = '090-WMAP_W2'
        if i in (28, 29):
            wmap[i] = [1.61, 9.02]
            gsfc[i] = [1.38, 1.76]
        elif i in (30, 31):
            wmap[i] = [1.72, 7.47]
            gsfc[i] = [1.44, 0.77]
    elif i < 36:
        name = 'W3'
        ind2chan[i] = '090-WMAP_W3'
        if i in (32, 33):
            wmap[i] = [1.65, 0.93]
            gsfc[i] = [1.47, 1.84]
        elif i in (34, 35):
            wmap[i] = [1.86, 0.28]
            gsfc[i] = [1.69, 2.39]
    else:
        name = 'W4'
        ind2chan[i] = '090-WMAP_W4'
        if i in (36, 37):
            wmap[i] = [1.71, 46.6]
            gsfc[i] = [1.60, 8.46]
        elif i in (38, 39):
            wmap[i] = [1.65, 26.0]
            gsfc[i] = [1.43, 5.31]
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
plot_text_coords = ([[52200, 0.835]] * 4 +  # K
    [[52200, 0.755]] * 4 +  # Ka
    [[52200, 1.11]] * 4 + # Q1
    [[52200, 1.005]] * 4 +  # Q2
    [[52200, 1.44]] * 4 +  # V1
    [[52200, 1.375]] * 4 +  # V2
    [[52200, 2.1]] * 4 + # W1
    [[52200, 2.1]] * 4 +  # W2
    [[52200, 2.4]] * 4 +  # W3
    [[52200, 2.35]] * 4 # W4
)

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
    if channel in ('023-WMAP_K', '030-WMAP_Ka'):
        factor = 12
    elif 'Q' in channel:
        factor = 15
    elif 'V' in channel:
        factor = 20
    elif 'WMAP_W' in channel:
        factor = 30
    else:
        print(channel)
        print('OOPs')
    samples1 = chain1.get(f'tod/{channel}/xi_n')[:, 0, :, :] / chain1.get(f'tod/{channel}/gain') * np.sqrt(1.536 / factor) / np.sqrt(2.0)
    samples2 = chain2.get(f'tod/{channel}/xi_n')[:, 0, :, :] / chain2.get(f'tod/{channel}/gain') * np.sqrt(1.536 / factor) / np.sqrt(2.0)
    tot_samples = np.concatenate((samples1[50:],  samples2[50:]), axis=0)
    mean = np.mean(tot_samples, axis=0)
    std = np.std(tot_samples, axis=0)
    accept = chain1.get(f'tod/{channel}/accept')[-1, i % 4, :].astype(bool)
    mjd = chain1.get(f'tod/{channel}/MJD')[-1, :]
    plt.errorbar(mjd[accept], mean[i%4, :][accept], yerr=std[i%4][accept],
        linewidth=0.5, color='black', zorder=1, label = 'CG', rasterized=True)
    plt.plot(mjd_wmap,[wmap[i][0],wmap[i][0]], linewidth=1, color='red',
        linestyle=':', label='WMAP', rasterized=True)
    plt.plot(mjd_gsfc,[gsfc[i][0],gsfc[i][0]], linewidth=1, color='orange',
        linestyle=':', label='GSFC', rasterized=True)


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
    if i % 4 == 0:
        plt.ylabel(r"$\sigma_{0}$ [mK\,$\mathrm{s}^{\frac{1}{2}}$]");
        currax.yaxis.labelpad = 10*width/17.
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
    if i == 4:
        plt.yticks([0.73,0.77], [r"$0.73$", r"$0.77$"])
    elif i == 36:
        plt.yticks([1.6,2.0,2.4], [r"$1.6$", r"$2.0$", r"$2.4$"])

for i in range(0, len(axs), 4):
    for ticklabel in axs[i].yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")

# save to pdf with right bounding box
plt.savefig("../../figures/instpar_CG_sigma0_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
