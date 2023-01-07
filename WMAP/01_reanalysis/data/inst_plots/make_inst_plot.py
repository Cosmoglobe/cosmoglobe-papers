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

wmap = np.loadtxt('regressed_gains.txt')


vmin = -110
vmax =  160
#data_Q = N.log10(0.5*(data_Q+N.sqrt(4.+data_Q*data_Q)))
#data_U = N.log10(0.5*(data_U+N.sqrt(4.+data_U*data_U)))
#data_Q = N.minimum(N.maximum(data_Q,vmin),vmax)
#data_U = N.minimum(N.maximum(data_U,vmin),vmax)


# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 1.35*cm2inch(width)))
# this should be changed for making a panel of multiple figures
#ax1 = fig.add_subplot(211)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



###############
#   K-band
###############

mjd_wmap = dataK[0,0] + np.arange(1,3280+1)/3280. * (dataK[-1,0]-dataK[0,0])

wmap2 = wmap[0,:]
inds = (np.abs(wmap2) < 0.01)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax1 = plt.subplot2grid((10, 4), (0, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=0.5, color='red')
plt.plot(dataK[:,0],dataK[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.25,r"K113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax1.yaxis.labelpad = 10*width/17.

wmap2 = wmap[1,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax2 = plt.subplot2grid((10, 4), (0, 1), sharey=ax1)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=0.5, color='red')
plt.plot(dataK[:,0],dataK[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.25,r"K114", fontsize=10)

wmap2 = wmap[2,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax3 = plt.subplot2grid((10, 4), (0, 2), sharey=ax1)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=0.5, color='red')
plt.plot(dataK[:,0],dataK[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.25,r"K123", fontsize=10)

wmap2 = wmap[3,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax4 = plt.subplot2grid((10, 4), (0, 3), sharey=ax1)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=0.5, color='red')
plt.plot(dataK[:,0],dataK[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.25,r"K124", fontsize=10)


###############
#   Ka-band
###############
mjd_wmap = dataKa[0,0] + np.arange(1,3280+1)/3280. * (dataKa[-1,0]-dataKa[0,0])

wmap2 = wmap[4,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax5 = plt.subplot2grid((10, 4), (1, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataKa[:,0],dataKa[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.16,r"Ka113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax5.yaxis.labelpad = 10*width/17.

wmap2 = wmap[5,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax6 = plt.subplot2grid((10, 4), (1, 1), sharey=ax5)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataKa[:,0],dataKa[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.16,r"Ka114", fontsize=10)

wmap2 = wmap[6,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax7 = plt.subplot2grid((10, 4), (1, 2), sharey=ax5)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataKa[:,0],dataKa[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=False)
plt.setp( ax7.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.16,r"Ka123", fontsize=10)

wmap2 = wmap[7,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax8 = plt.subplot2grid((10, 4), (1, 3), sharey=ax5)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataKa[:,0],dataKa[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax8.get_xticklabels(), visible=False)
plt.setp( ax8.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.16,r"Ka124", fontsize=10)

###############
#   Q1-band
###############

mjd_wmap = dataQ1[0,0] + np.arange(1,3280+1)/3280. * (dataQ1[-1,0]-dataQ1[0,0])
wmap2 = wmap[8,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax9 = plt.subplot2grid((10, 4), (2, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ1[:,0],dataQ1[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax9.get_xticklabels(), visible=False)
plt.setp( ax9.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"Q113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax9.yaxis.labelpad = 10*width/17.

wmap2 = wmap[9,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax10 = plt.subplot2grid((10, 4), (2, 1), sharey=ax9)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ1[:,0],dataQ1[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax10.get_xticklabels(), visible=False)
plt.setp( ax10.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"Q114", fontsize=10)

wmap2 = wmap[10,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax11 = plt.subplot2grid((10, 4), (2, 2), sharey=ax9)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ1[:,0],dataQ1[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax11.get_xticklabels(), visible=False)
plt.setp( ax11.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.06,r"Q123", fontsize=10)

wmap2 = wmap[11,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax12 = plt.subplot2grid((10, 4), (2, 3), sharey=ax9)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ1[:,0],dataQ1[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax12.get_xticklabels(), visible=False)
plt.setp( ax12.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.06,r"Q124", fontsize=10)


###############
#   Q2-band
###############

mjd_wmap = dataQ2[0,0] + np.arange(1,3280+1)/3280. * (dataQ2[-1,0]-dataQ2[0,0])
wmap2 = wmap[12,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
ax13 = plt.subplot2grid((10, 4), (3, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ2[:,0],dataQ2[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax13.get_xticklabels(), visible=False)
plt.setp( ax13.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.10,r"Q211", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax13.yaxis.labelpad = 10*width/17.

wmap2 = wmap[13,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax14 = plt.subplot2grid((10, 4), (3, 1), sharey=ax13)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ2[:,0],dataQ2[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax14.get_xticklabels(), visible=False)
plt.setp( ax14.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.10,r"Q212", fontsize=10)

wmap2 = wmap[14,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax15 = plt.subplot2grid((10, 4), (3, 2), sharey=ax13)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ2[:,0],dataQ2[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax15.get_xticklabels(), visible=False)
plt.setp( ax15.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.10,r"Q223", fontsize=10)

wmap2 = wmap[15,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax16 = plt.subplot2grid((10, 4), (3, 3), sharey=ax13)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataQ2[:,0],dataQ2[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax16.get_xticklabels(), visible=False)
plt.setp( ax16.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.10,r"Q224", fontsize=10)


###############
#   V1-band
###############

mjd_wmap = dataV1[0,0] + np.arange(1,3280+1)/3280. * (dataV1[-1,0]-dataV1[0,0])
wmap2 = wmap[16,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax17 = plt.subplot2grid((10, 4), (4, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV1[:,0],dataV1[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax17.get_xticklabels(), visible=False)
plt.setp( ax17.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"V113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax17.yaxis.labelpad = 10*width/17.

wmap2 = wmap[17,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax18 = plt.subplot2grid((10, 4), (4, 1), sharey=ax17)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV1[:,0],dataV1[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax18.get_xticklabels(), visible=False)
plt.setp( ax18.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"V114", fontsize=10)

wmap2 = wmap[18,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax19 = plt.subplot2grid((10, 4), (4, 2), sharey=ax17)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV1[:,0],dataV1[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax19.get_xticklabels(), visible=False)
plt.setp( ax19.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"V123", fontsize=10)

wmap2 = wmap[19,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]
        
ax20 = plt.subplot2grid((10, 4), (4, 3), sharey=ax17)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV1[:,0],dataV1[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax20.get_xticklabels(), visible=False)
plt.setp( ax20.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.55,r"V124", fontsize=10)


###############
#   V2-band
###############

mjd_wmap = dataV2[0,0] + np.arange(1,3280+1)/3280. * (dataV1[-1,0]-dataV1[0,0])
wmap2 = wmap[20,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax21 = plt.subplot2grid((10, 4), (5, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV2[:,0],dataV2[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax21.get_xticklabels(), visible=False)
plt.setp( ax21.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.46,r"V211", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax21.yaxis.labelpad = 10*width/17.

wmap2 = wmap[21,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax22 = plt.subplot2grid((10, 4), (5, 1), sharey=ax21)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV2[:,0],dataV2[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax22.get_xticklabels(), visible=False)
plt.setp( ax22.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.46,r"V212", fontsize=10)

wmap2 = wmap[22,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax23 = plt.subplot2grid((10, 4), (5, 2), sharey=ax21)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV2[:,0],dataV2[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax23.get_xticklabels(), visible=False)
plt.setp( ax23.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.46,r"V223", fontsize=10)

wmap2 = wmap[23,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax24 = plt.subplot2grid((10, 4), (5, 3), sharey=ax21)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataV2[:,0],dataV2[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax24.get_xticklabels(), visible=False)
plt.setp( ax24.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.46,r"V224", fontsize=10)



###############
#   W1-band
###############

mjd_wmap = dataW1[0,0] + np.arange(1,3280+1)/3280. * (dataW1[-1,0]-dataW1[0,0])
wmap2 = wmap[24,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax25 = plt.subplot2grid((10, 4), (6, 0))
plt.ylim([0.2, 0.4])
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW1[:,0],dataW1[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax25.get_xticklabels(), visible=False)
plt.setp( ax25.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.25,r"W113", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax25.yaxis.labelpad = 10*width/17.

wmap2 = wmap[25,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax26 = plt.subplot2grid((10, 4), (6, 1), sharey=ax25)
plt.ylim([0.2, 0.4])
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW1[:,0],dataW1[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax26.get_xticklabels(), visible=False)
plt.setp( ax26.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.25,r"W114", fontsize=10)

wmap2 = wmap[26,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax27 = plt.subplot2grid((10, 4), (6, 2), sharey=ax25)
plt.ylim([0.2, 0.4])
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW1[:,0],dataW1[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax27.get_xticklabels(), visible=False)
plt.setp( ax27.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.31,r"W123", fontsize=10)

wmap2 = wmap[27,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax28 = plt.subplot2grid((10, 4), (6, 3), sharey=ax25)
plt.ylim([0.2, 0.4])
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW1[:,0],dataW1[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax28.get_xticklabels(), visible=False)
plt.setp( ax28.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.31,r"W124", fontsize=10)


###############
#   W2-band
###############

mjd_wmap = dataW2[0,0] + np.arange(1,3280+1)/3280. * (dataW2[-1,0]-dataW2[0,0])
wmap2 = wmap[28,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax29 = plt.subplot2grid((10, 4), (7, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW2[:,0],dataW2[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax29.get_xticklabels(), visible=False)
plt.setp( ax29.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.305,r"W211", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax29.yaxis.labelpad = 10*width/17.

wmap2 = wmap[29,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax30 = plt.subplot2grid((10, 4), (7, 1), sharey=ax29)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW2[:,0],dataW2[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax30.get_xticklabels(), visible=False)
plt.setp( ax30.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.305,r"W212", fontsize=10)

wmap2 = wmap[30,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax31 = plt.subplot2grid((10, 4), (7, 2), sharey=ax29)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW2[:,0],dataW2[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax31.get_xticklabels(), visible=False)
plt.setp( ax31.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.305,r"W223", fontsize=10)

wmap2 = wmap[31,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax32 = plt.subplot2grid((10, 4), (7, 3), sharey=ax29)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW2[:,0],dataW2[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax32.get_xticklabels(), visible=False)
plt.setp( ax32.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.305,r"W224", fontsize=10)


###############
#   W3-band
###############

mjd_wmap = dataW3[0,0] + np.arange(1,3280+1)/3280. * (dataW3[-1,0]-dataW3[0,0])
wmap2 = wmap[32,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax33 = plt.subplot2grid((10, 4), (8, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW3[:,0],dataW3[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax33.get_xticklabels(), visible=False)
plt.setp( ax33.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.282,r"W311", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax33.yaxis.labelpad = 10*width/17.

wmap2 = wmap[33,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax34 = plt.subplot2grid((10, 4), (8, 1), sharey=ax33)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW3[:,0],dataW3[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax34.get_xticklabels(), visible=False)
plt.setp( ax34.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.282,r"W312", fontsize=10)

wmap2 = wmap[34,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax35 = plt.subplot2grid((10, 4), (8, 2), sharey=ax33)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW3[:,0],dataW3[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax35.get_xticklabels(), visible=False)
plt.setp( ax35.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.282,r"W323", fontsize=10)

wmap2 = wmap[35,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax36 = plt.subplot2grid((10, 4), (8, 3), sharey=ax33)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW3[:,0],dataW3[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax36.get_xticklabels(), visible=False)
plt.setp( ax36.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.282,r"W324", fontsize=10)


###############
#   W4-band
###############

mjd_wmap = dataW4[0,0] + np.arange(1,3280+1)/3280. * (dataW4[-1,0]-dataW4[0,0])
wmap2 = wmap[36,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

ax37 = plt.subplot2grid((10, 4), (9, 0))
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW4[:,0],dataW4[:,1], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax37.get_xticklabels(), visible=True)
plt.setp( ax37.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.286,r"W411", fontsize=10)
plt.ylabel(r"$g$ [du/mK]");
ax37.yaxis.labelpad = 10*width/17.

wmap2 = wmap[37,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax38 = plt.subplot2grid((10, 4), (9, 1), sharey=ax37)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW4[:,0],dataW4[:,2], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax38.get_xticklabels(), visible=True)
plt.setp( ax38.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.286,r"W412", fontsize=10)

wmap2 = wmap[38,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax39 = plt.subplot2grid((10, 4), (9, 2), sharey=ax37)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW4[:,0],dataW4[:,3], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax39.get_xticklabels(), visible=True)
plt.setp( ax39.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.24,r"W423", fontsize=10)

wmap2 = wmap[39,:]
inds = (np.abs(wmap2) < 0.1)
for i in range(len(inds)):
    if inds[i]:
        wmap2[i] = wmap2[i-1]

plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])        
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax40 = plt.subplot2grid((10, 4), (9, 3), sharey=ax37)
plt.plot(mjd_wmap,1./np.abs(wmap2), linewidth=1, color='red')
plt.plot(dataW4[:,0],dataW4[:,4], linewidth=1, color='black')
plt.grid(False, which="major", axis="both")
plt.setp( ax40.get_xticklabels(), visible=True)
plt.setp( ax40.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,0.24,r"W424", fontsize=10)

plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])
# labels
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels




# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax5.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax9.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax13.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax17.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax21.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax25.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax29.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax33.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax37.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")        



# save to pdf with right bounding box
plt.savefig("instpar_CG_gain_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

