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
dataK = np.loadtxt('fknee_CG_023-WMAP_K_v1.dat')
dataKa = np.loadtxt('fknee_CG_030-WMAP_Ka_v1.dat')
dataQ1 = np.loadtxt('fknee_CG_040-WMAP_Q1_v1.dat')
dataQ2 = np.loadtxt('fknee_CG_040-WMAP_Q2_v1.dat')
dataV1 = np.loadtxt('fknee_CG_060-WMAP_V1_v1.dat')
dataV2 = np.loadtxt('fknee_CG_060-WMAP_V2_v1.dat')
dataW1 = np.loadtxt('fknee_CG_090-WMAP_W1_v1.dat')
dataW2 = np.loadtxt('fknee_CG_090-WMAP_W2_v1.dat')
dataW3 = np.loadtxt('fknee_CG_090-WMAP_W3_v1.dat')
dataW4 = np.loadtxt('fknee_CG_090-WMAP_W4_v1.dat')

maskK = np.loadtxt('mask_CG_023-WMAP_K_v1.dat')
maskKa = np.loadtxt('mask_CG_030-WMAP_Ka_v1.dat')
maskQ1 = np.loadtxt('mask_CG_040-WMAP_Q1_v1.dat')
maskQ2 = np.loadtxt('mask_CG_040-WMAP_Q2_v1.dat')
maskV1 = np.loadtxt('mask_CG_060-WMAP_V1_v1.dat')
maskV2 = np.loadtxt('mask_CG_060-WMAP_V2_v1.dat')
maskW1 = np.loadtxt('mask_CG_090-WMAP_W1_v1.dat')
maskW2 = np.loadtxt('mask_CG_090-WMAP_W2_v1.dat')
maskW3 = np.loadtxt('mask_CG_090-WMAP_W3_v1.dat')
maskW4 = np.loadtxt('mask_CG_090-WMAP_W4_v1.dat')

dataK[:,1:5] = 1000*dataK[:,1:5]
dataKa[:,1:5] = 1000*dataKa[:,1:5]
dataQ1[:,1:5] = 1000*dataQ1[:,1:5]
dataQ2[:,1:5] = 1000*dataQ2[:,1:5]
dataV1[:,1:5] = 1000*dataV1[:,1:5]
dataV2[:,1:5] = 1000*dataV2[:,1:5]
dataW1[:,1:5] = 1000*dataW1[:,1:5]
dataW2[:,1:5] = 1000*dataW2[:,1:5]
dataW3[:,1:5] = 1000*dataW3[:,1:5]
dataW4[:,1:5] = 1000*dataW4[:,1:5]

inds = np.where(maskK == 0)
dataK[inds] = np.nan
inds = np.where(maskKa == 0)
dataKa[inds] = np.nan
inds = np.where(maskQ1 == 0)
dataQ1[inds] = np.nan
inds = np.where(maskQ2 == 0)
dataQ2[inds] = np.nan
inds = np.where(maskV1 == 0)
dataV1[inds] = np.nan
inds = np.where(maskV2 == 0)
dataV2[inds] = np.nan
inds = np.where(maskW1 == 0)
dataW1[inds] = np.nan
inds = np.where(maskW2 == 0)
dataW2[inds] = np.nan
inds = np.where(maskW3 == 0)
dataW3[inds] = np.nan
inds = np.where(maskW4 == 0)
dataW4[inds] = np.nan


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
#ax1 = fig.add_subplot(213)

fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)

wmap = N.zeros((20,2))
wmap[0,:] = [0.66, 0.40]
wmap[1,:] = [0.75, 0.51]
wmap[2,:] = [0.71, 0.71]
wmap[3,:] = [0.72, 0.32]
wmap[4,:] = [0.92, 1.09]
wmap[5,:] = [1.02, 0.35]
wmap[6,:] = [0.85, 5.76]
wmap[7,:] = [0.99, 8.62]
wmap[8,:] = [1.22, 0.09]
wmap[9,:] = [1.11, 1.41]
wmap[10,:] = [0.97, 0.88]
wmap[11,:] = [1.10, 8.35]
wmap[12,:] = [1.35, 7.88]
wmap[13,:] = [1.61, 0.66]
wmap[14,:] = [1.61, 9.02]
wmap[15,:] = [1.72, 7.47]
wmap[16,:] = [1.65, 0.93]
wmap[17,:] = [1.86, 0.28]
wmap[18,:] = [1.71, 46.5]
wmap[19,:] = [1.65, 26.0]

gsfc = N.zeros((20,2))
gsfc[0,:] = [0.72, 6.13]
gsfc[1,:] = [0.87, 5.37]
gsfc[2,:] = [0.75, 1.66]
gsfc[3,:] = [0.77, 1.29]
gsfc[4,:] = [0.99, 3.21]
gsfc[5,:] = [0.95, 3.13]
gsfc[6,:] = [0.89, 1.92]
gsfc[7,:] = [1.04, 4.61]
gsfc[8,:] = [1.25, 2.56]
gsfc[9,:] = [1.07, 4.49]
gsfc[10,:] = [1.01, 2.43]
gsfc[11,:] = [1.13, 3.06]
gsfc[12,:] = [1.18, 16.2]
gsfc[13,:] = [1.41, 15.1]
gsfc[14,:] = [1.38, 1.76]
gsfc[15,:] = [1.44, 0.77]
gsfc[16,:] = [1.47, 1.84]
gsfc[17,:] = [1.69, 2.39]
gsfc[18,:] = [1.60, 8.46]
gsfc[19,:] = [1.43, 5.31]

mjd_wmap = [52130, 52477]
mjd_gsfc = [52130, 55414]


###############
#   K-band
###############

ax1 = plt.subplot2grid((10, 4), (0, 0))
plt.plot(dataK[:,0],dataK[:,1], linewidth=1, color='black', label='CG')
plt.plot(mjd_wmap,[wmap[0,1],wmap[0,1]], linewidth=1, color='red', linestyle=':', label='WMAP')
plt.plot(mjd_gsfc,[gsfc[0,1],gsfc[0,1]], linewidth=1, color='orange', linestyle=':', label='GSFC')
plt.grid(False, which="major", axis="both")
plt.setp( ax1.get_xticklabels(), visible=False)
plt.setp( ax1.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,4.5,r"K113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax1.yaxis.labelpad = 10*width/17.
plt.yticks([2,4,6], [r"$2$", r"$4$", r"$6$"])

ax2 = plt.subplot2grid((10, 4), (0, 1), sharey=ax1)
plt.plot(dataK[:,0],dataK[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[0,1],wmap[0,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[1,1],gsfc[1,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax2.get_xticklabels(), visible=False)
plt.setp( ax2.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,4.5,r"K114", fontsize=10)

ax3 = plt.subplot2grid((10, 4), (0, 2), sharey=ax1)
plt.plot(dataK[:,0],dataK[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[1,1],wmap[1,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[1,1],gsfc[1,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax3.get_xticklabels(), visible=False)
plt.setp( ax3.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,4.5,r"K123", fontsize=10)

ax4 = plt.subplot2grid((10, 4), (0, 3), sharey=ax1)
plt.plot(dataK[:,0],dataK[:,4], linewidth=1, color='black', label='CG')
plt.plot(mjd_wmap,[wmap[1,1],wmap[1,1]], linewidth=1, color='red', linestyle=':', label='WMAP')
plt.plot(mjd_gsfc,[gsfc[1,1],gsfc[1,1]], linewidth=1, color='orange', linestyle=':', label='GSFC')
plt.grid(False, which="major", axis="both")
plt.setp( ax4.get_xticklabels(), visible=False)
plt.setp( ax4.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,4.5,r"K124", fontsize=10)

# legend
leg = plt.legend(frameon=True, loc=1, fontsize=8)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)



###############
#   Ka-band
###############
ax5 = plt.subplot2grid((10, 4), (1, 0))
plt.plot(dataKa[:,0],dataKa[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[2,1],wmap[2,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[2,1],gsfc[2,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax5.get_xticklabels(), visible=False)
plt.setp( ax5.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.8,r"Ka113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax5.yaxis.labelpad = 10*width/17.

ax6 = plt.subplot2grid((10, 4), (1, 1), sharey=ax5)
plt.plot(dataKa[:,0],dataKa[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[2,1],wmap[2,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[2,1],gsfc[2,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax6.get_xticklabels(), visible=False)
plt.setp( ax6.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.8,r"Ka114", fontsize=10)

ax7 = plt.subplot2grid((10, 4), (1, 2), sharey=ax5)
plt.plot(dataKa[:,0],dataKa[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[3,1],wmap[3,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[3,1],gsfc[3,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax7.get_xticklabels(), visible=False)
plt.setp( ax7.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.8,r"Ka123", fontsize=10)

ax8 = plt.subplot2grid((10, 4), (1, 3), sharey=ax5)
plt.plot(dataKa[:,0],dataKa[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[3,1],wmap[3,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[3,1],gsfc[3,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax8.get_xticklabels(), visible=False)
plt.setp( ax8.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,1.8,r"Ka124", fontsize=10)

###############
#   Q1-band
###############

ax9 = plt.subplot2grid((10, 4), (2, 0))
plt.plot(dataQ1[:,0],dataQ1[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[4,1],wmap[4,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[4,1],gsfc[4,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax9.get_xticklabels(), visible=False)
plt.setp( ax9.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.6,r"Q113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax9.yaxis.labelpad = 10*width/17.

ax10 = plt.subplot2grid((10, 4), (2, 1), sharey=ax9)
plt.plot(dataQ1[:,0],dataQ1[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[4,1],wmap[4,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[4,1],gsfc[4,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax10.get_xticklabels(), visible=False)
plt.setp( ax10.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.6,r"Q114", fontsize=10)

ax11 = plt.subplot2grid((10, 4), (2, 2), sharey=ax9)
plt.plot(dataQ1[:,0],dataQ1[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[5,1],wmap[5,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[5,1],gsfc[5,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax11.get_xticklabels(), visible=False)
plt.setp( ax11.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.6,r"Q123", fontsize=10)

ax12 = plt.subplot2grid((10, 4), (2, 3), sharey=ax9)
plt.plot(dataQ1[:,0],dataQ1[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[5,1],wmap[5,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[5,1],gsfc[5,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax12.get_xticklabels(), visible=False)
plt.setp( ax12.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.6,r"Q124", fontsize=10)


###############
#   Q2-band
###############

ax13 = plt.subplot2grid((10, 4), (3, 0))
plt.plot(dataQ2[:,0],dataQ2[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[6,1],wmap[6,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[6,1],gsfc[6,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax13.get_xticklabels(), visible=False)
plt.setp( ax13.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,7.5,r"Q213", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax13.yaxis.labelpad = 10*width/17.

ax14 = plt.subplot2grid((10, 4), (3, 1), sharey=ax13)
plt.plot(dataQ2[:,0],dataQ2[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[6,1],wmap[6,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[6,1],gsfc[6,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax14.get_xticklabels(), visible=False)
plt.setp( ax14.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,7.5,r"Q214", fontsize=10)

ax15 = plt.subplot2grid((10, 4), (3, 2), sharey=ax13)
plt.plot(dataQ2[:,0],dataQ2[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[7,1],wmap[7,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[7,1],gsfc[7,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax15.get_xticklabels(), visible=False)
plt.setp( ax15.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,7.5,r"Q223", fontsize=10)

ax16 = plt.subplot2grid((10, 4), (3, 3), sharey=ax13)
plt.plot(dataQ2[:,0],dataQ2[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[7,1],wmap[7,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[7,1],gsfc[7,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax16.get_xticklabels(), visible=False)
plt.setp( ax16.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,7.5,r"Q224", fontsize=10)


###############
#   V1-band
###############

ax17 = plt.subplot2grid((10, 4), (4, 0))
plt.plot(dataV1[:,0],dataV1[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[8,1],wmap[8,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[8,1],gsfc[8,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax17.get_xticklabels(), visible=False)
plt.setp( ax17.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,3.9,r"V113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax17.yaxis.labelpad = 10*width/17.

ax18 = plt.subplot2grid((10, 4), (4, 1), sharey=ax17)
plt.plot(dataV1[:,0],dataV1[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[8,1],wmap[8,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[8,1],gsfc[8,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax18.get_xticklabels(), visible=False)
plt.setp( ax18.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,3.9,r"V114", fontsize=10)

ax19 = plt.subplot2grid((10, 4), (4, 2), sharey=ax17)
plt.plot(dataV1[:,0],dataV1[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[9,1],wmap[9,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[9,1],gsfc[9,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax19.get_xticklabels(), visible=False)
plt.setp( ax19.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,3.9,r"V123", fontsize=10)

ax20 = plt.subplot2grid((10, 4), (4, 3), sharey=ax17)
plt.plot(dataV1[:,0],dataV1[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[9,1],wmap[9,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[9,1],gsfc[9,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax20.get_xticklabels(), visible=False)
plt.setp( ax20.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,3.9,r"V124", fontsize=10)


###############
#   V2-band
###############

ax21 = plt.subplot2grid((10, 4), (5, 0))
plt.plot(dataV2[:,0],dataV2[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[10,1],wmap[10,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[10,1],gsfc[10,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax21.get_xticklabels(), visible=False)
plt.setp( ax21.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,6.9,r"V213", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax21.yaxis.labelpad = 10*width/17.

ax22 = plt.subplot2grid((10, 4), (5, 1), sharey=ax21)
plt.plot(dataV2[:,0],dataV2[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[10,1],wmap[10,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[10,1],gsfc[10,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax22.get_xticklabels(), visible=False)
plt.setp( ax22.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,6.9,r"V214", fontsize=10)

ax23 = plt.subplot2grid((10, 4), (5, 2), sharey=ax21)
plt.plot(dataV2[:,0],dataV2[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[11,1],wmap[11,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[11,1],gsfc[11,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax23.get_xticklabels(), visible=False)
plt.setp( ax23.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,6.9,r"V223", fontsize=10)

ax24 = plt.subplot2grid((10, 4), (5, 3), sharey=ax21)
plt.plot(dataV2[:,0],dataV2[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[11,1],wmap[11,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[11,1],gsfc[11,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax24.get_xticklabels(), visible=False)
plt.setp( ax24.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,6.9,r"V224", fontsize=10)



###############
#   W1-band
###############

ax25 = plt.subplot2grid((10, 4), (6, 0))
plt.plot(dataW1[:,0],dataW1[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[12,1],wmap[12,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[12,1],gsfc[12,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax25.get_xticklabels(), visible=False)
plt.setp( ax25.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200, 13,r"W113", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax25.yaxis.labelpad = 10*width/17.

plt.yticks([5,10,15], [r"$5$", r"$10$", r"$15$"])

ax26 = plt.subplot2grid((10, 4), (6, 1), sharey=ax25)
plt.plot(dataW1[:,0],dataW1[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[12,1],wmap[12,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[12,1],gsfc[12,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax26.get_xticklabels(), visible=False)
plt.setp( ax26.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,13,r"W114", fontsize=10)

ax27 = plt.subplot2grid((10, 4), (6, 2), sharey=ax25)
plt.plot(dataW1[:,0],dataW1[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[13,1],wmap[13,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[13,1],gsfc[13,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax27.get_xticklabels(), visible=False)
plt.setp( ax27.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,13,r"W123", fontsize=10)

ax28 = plt.subplot2grid((10, 4), (6, 3), sharey=ax25)
plt.plot(dataW1[:,0],dataW1[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[13,1],wmap[13,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[13,1],gsfc[13,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax28.get_xticklabels(), visible=False)
plt.setp( ax28.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,13,r"W124", fontsize=10)


###############
#   W2-band
###############

ax29 = plt.subplot2grid((10, 4), (7, 0))
plt.plot(dataW2[:,0],dataW2[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[14,1],wmap[14,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[14,1],gsfc[14,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax29.get_xticklabels(), visible=False)
plt.setp( ax29.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,16.7,r"W213", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax29.yaxis.labelpad = 10*width/17.

ax30 = plt.subplot2grid((10, 4), (7, 1), sharey=ax29)
plt.plot(dataW2[:,0],dataW2[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[14,1],wmap[14,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[14,1],gsfc[14,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax30.get_xticklabels(), visible=False)
plt.setp( ax30.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,16.7,r"W214", fontsize=10)

ax31 = plt.subplot2grid((10, 4), (7, 2), sharey=ax29)
plt.plot(dataW2[:,0],dataW2[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[15,1],wmap[15,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[15,1],gsfc[15,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax31.get_xticklabels(), visible=False)
plt.setp( ax31.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,16.7,r"W223", fontsize=10)

ax32 = plt.subplot2grid((10, 4), (7, 3), sharey=ax29)
plt.plot(dataW2[:,0],dataW2[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[15,1],wmap[15,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[15,1],gsfc[15,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax32.get_xticklabels(), visible=False)
plt.setp( ax32.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,16.7,r"W224", fontsize=10)


###############
#   W3-band
###############

ax33 = plt.subplot2grid((10, 4), (8, 0))
plt.plot(dataW3[:,0],dataW3[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[16,1],wmap[16,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[16,1],gsfc[16,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax33.get_xticklabels(), visible=False)
plt.setp( ax33.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.76,r"W313", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax33.yaxis.labelpad = 10*width/17.

ax34 = plt.subplot2grid((10, 4), (8, 1), sharey=ax33)
plt.plot(dataW3[:,0],dataW3[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[16,1],wmap[16,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[16,1],gsfc[16,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax34.get_xticklabels(), visible=False)
plt.setp( ax34.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.76,r"W314", fontsize=10)

ax35 = plt.subplot2grid((10, 4), (8, 2), sharey=ax33)
plt.plot(dataW3[:,0],dataW3[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[17,1],wmap[17,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[17,1],gsfc[17,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax35.get_xticklabels(), visible=False)
plt.setp( ax35.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.76,r"W323", fontsize=10)

ax36 = plt.subplot2grid((10, 4), (8, 3), sharey=ax33)
plt.plot(dataW3[:,0],dataW3[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[17,1],wmap[17,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[17,1],gsfc[17,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax36.get_xticklabels(), visible=False)
plt.setp( ax36.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,2.76,r"W324", fontsize=10)


###############
#   W4-band
###############

ax37 = plt.subplot2grid((10, 4), (9, 0))
plt.plot(dataW4[:,0],dataW4[:,1], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[18,1],wmap[18,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[18,1],gsfc[18,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax37.get_xticklabels(), visible=True)
plt.setp( ax37.get_yticklabels(), visible=True)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,40,r"W413", fontsize=10)
plt.ylabel(r"$f_{\mathrm{knee}}$ [mHz]");
ax37.yaxis.labelpad = 10*width/17.

plt.xticks([52000,53000,54000,55000], [r"$52\,000$", r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax38 = plt.subplot2grid((10, 4), (9, 1), sharey=ax37)
plt.plot(dataW4[:,0],dataW4[:,2], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[18,1],wmap[18,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[18,1],gsfc[18,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax38.get_xticklabels(), visible=True)
plt.setp( ax38.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,40,r"W414", fontsize=10)

plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax39 = plt.subplot2grid((10, 4), (9, 2), sharey=ax37)
plt.plot(dataW4[:,0],dataW4[:,3], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[19,1],wmap[19,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[19,1],gsfc[19,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax39.get_xticklabels(), visible=True)
plt.setp( ax39.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,40,r"W423", fontsize=10)

plt.xticks([53000,54000,55000], [r"$53\,000$", r"$54\,000$", r"$55\,000$"])        
plt.xlabel(r"MJD");
ax3.yaxis.labelpad = 10*width/17.; ax3.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax40 = plt.subplot2grid((10, 4), (9, 3), sharey=ax37)
plt.plot(dataW4[:,0],dataW4[:,4], linewidth=1, color='black')
plt.plot(mjd_wmap,[wmap[19,1],wmap[19,1]], linewidth=1, color='red', linestyle=':')
plt.plot(mjd_gsfc,[gsfc[19,1],gsfc[19,1]], linewidth=1, color='orange', linestyle=':')
plt.grid(False, which="major", axis="both")
plt.setp( ax40.get_xticklabels(), visible=True)
plt.setp( ax40.get_yticklabels(), visible=False)
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.text(52200,40,r"W424", fontsize=10)

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
plt.savefig("../../figures/instpar_CG_fknee_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

