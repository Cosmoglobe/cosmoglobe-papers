# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp

import plotly.colors as pcol
import matplotlib as mpl

from glob import glob

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 19.0

# Load data
dataK = np.loadtxt('gain_CG_023-WMAP_K_v1.dat')
dataKa = np.loadtxt('gain_CG_030-WMAP_Ka_v1.dat')
dataQ1 = np.loadtxt('gain_CG_040-WMAP_Q1_v1.dat')
dataQ2 = np.loadtxt('gain_CG_040-WMAP_Q2_v1.dat')
dataV1 = np.loadtxt('gain_CG_060-WMAP_V1_v1.dat')
dataV2 = np.loadtxt('gain_CG_060-WMAP_V2_v1.dat')
dataW1 = np.loadtxt('gain_CG_090-WMAP_W1_v1.dat')
dataW2 = np.loadtxt('gain_CG_090-WMAP_W2_v1.dat')
dataW3 = np.loadtxt('gain_CG_090-WMAP_W3_v1.dat')
dataW4 = np.loadtxt('gain_CG_090-WMAP_W4_v1.dat')

fnames = glob('*_g0.txt')
fnames.sort()
thin = 73


vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 1.35*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



###############
#   K-band
###############

axs = []
wmap_data = [dataK, dataKa, dataQ1, dataQ2, dataV1, dataV2, dataW1, dataW2, dataW3, dataW4]
plot_texts = [r"K113", r"K114", r"K123", r"K124", r"Ka113", r"Ka114", r"Ka123", r"Ka124", r"Q113", r"Q114", r"Q123", r"Q124", r"Q213", r"Q214", r"Q223", r"Q224", r"V113", r"V114", r"V123", r"V124", r"V213", r"V214", r"V223", r"V224", r"W113", r"W114", r"W123", r"W124", r"W213", r"W214", r"W223", r"W224", r"W313", r"W314", r"W323", r"W324", r"W413", r"W414", r"W423", r"W424"]
plot_text_coords = ([[52200, 1.25]] * 4 +  # K
    [[52200, 1.16]] * 4 +  # Ka
    [[52200, 0.55]] * 2 +  [[52200, 1.06]] * 2 +  # Q1
    [[52200, 1.10]] * 4 +  # Q2
    [[52200, 0.55]] * 4 +  # V1
    [[52200, 0.46]] * 4 +  # V2
    [[52200, 0.25]] * 2 + [[52200, 0.31]] * 2 + # W1
    [[52200, 0.305]] * 4 +  # W2
    [[52200, 0.282]] * 4 +  # W3
    [[52200, 0.286]] * 2 + [[52200, 0.24]] * 2 # W4
)

for i in range(40):
    print(i)
    print(plot_text_coords[i])
#ax1 = plt.subplot2grid((10, 4), (0, 0))
    if i % 4 == 0:
        currax = plt.subplot2grid((10, 4), (i//4, i % 4))
        firstax = currax
#        axs.append(plt.subplot2grid((10, 4), (i//4, i%4)))
    else:
#        axs.append(plt.subplot2grid((10, 4), (i//4, i%4), sharey=axs[i//4]))
        currax = plt.subplot2grid((10, 4), (i//4, i%4), sharey=firstax)
    axs.append(currax)
    data = np.loadtxt(fnames[i])
    curr_wmap_data = wmap_data[i//4]
    plt.plot(curr_wmap_data[:,0],curr_wmap_data[:,(i % 4) + 1], linewidth=1, color='black')
#    mean = np.mean(abs(data[:, 1]))
#    std = np.std(abs(data[:, 1]))
    plt.plot(data[::thin,0], abs(data[::thin,1]), linewidth=0.5, color='red')
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
        plt.ylabel(r"$g$ [du/mK]");
        currax.yaxis.labelpad = 10*width/17.
    if i >= 36:
        if i % 4 == 0:
            plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
        else:
            plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])
        ## labels
        plt.xlabel(r"MJD");

for i in range(0, len(axs), 4):
    for ticklabel in axs[i].yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")



# save to pdf with right bounding box
plt.savefig("../../figures/instpar_CG_gain_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
