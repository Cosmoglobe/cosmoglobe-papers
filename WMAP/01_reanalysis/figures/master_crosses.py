import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cosmoglobe as cg
import astropy.units as u

from glob import glob

# Import the NaMaster python wrapper
import pymaster as nmt

#  Simple example showcasing the use of NaMaster to compute the pseudo-Cl
#  estimator of the angular cross-power spectrum of a spin-0 field and a
#  spin-2 field

# HEALPix resolution parameter used here
nside = 512

# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map(
    "/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/data/wmap_kq75_TQU_mask_r9.fits"
)


# bin4 = nmt.NmtBin.from_edges(l_ini, l_end)
logbins = np.unique(np.geomspace(10, 3 * 512 - 1, 100).astype("int"))
l_ini = np.concatenate((np.arange(2, 10), logbins[:-1]))
l_end = np.concatenate((np.arange(2, 10) + 1, logbins[1:]))
#b = nmt.NmtBin.from_nside_linear(nside, 1)
b = nmt.NmtBin.from_edges(l_ini, l_end)
ell_eff = b.get_effective_ells()
print(ell_eff)
print(b)
print(ell_eff.shape)

DIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"

wmap_maps = [
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_K1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Ka1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W3_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W4_v5.fits",
]
CG_DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

cg_maps = [
    f"{CG_DIR}/CG_023-WMAP_K_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_030-WMAP_Ka_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_040-WMAP_Q1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_040-WMAP_Q2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_060-WMAP_V1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_060-WMAP_V2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W3_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W4_IQU_n0512_v1.fits",
]

bands = ["K", "Ka", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]

diff_labs = ['KKa', 'Q', 'V', 'W']

x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip_W = -0.233 * x - 2.222 * y + 2.504 * z

ms = np.zeros((10, 3, hp.nside2npix(512)))
for i in range(10):
    ms[i] = hp.read_map(cg_maps[i], field=(0,1,2))*1e3
    ms[i][0] -= dip_W*1e3
    #cg.plot(ms[i], min=-0.25, max=0.25)


fig1, axes1 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=3)
fig2, axes2 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=3)
fig3, axes3 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=3)
diffs = np.zeros((4, 3, hp.nside2npix(512)))
diffs[0] = (ms[1] - 0.32*ms[0])/2
diffs[1] = (ms[3] - ms[2])/2
diffs[2] = (ms[5] - ms[4])/2
diffs[3] = (ms[9] - ms[8] - (ms[7] - ms[6]))/4

#for i in range(len(diffs)):
#    cg.plot(diffs[i], sig=1, fwhm=2*u.deg)
#plt.show()

fs  = [nmt.NmtField(mask, diffs[i][1:]) for i in range(len(diffs))]
fTs = [nmt.NmtField(mask, [diffs[i][0]]) for i in range(len(diffs))]
Clhat_iis = []
ClhatT_iis = []
for i in range(len(diffs)):
    Clhat_iis.append(nmt.compute_full_master(fs[i], fs[i], b))
    ClhatT_iis.append(nmt.compute_full_master(fTs[i], fTs[i], b))
for i in range(len(diffs)):
    print(i)
    for j in range(1,i):
        print(i,j, 'off')
        #if j > 0: 
        #    axes1[i,j-1].axis('off')
        #    axes2[i,j-1].axis('off')
        #    axes3[i,j-1].axis('off')
    for j in range(i+1, len(diffs)):
        print(i,j)
        ClhatT_ij = nmt.compute_full_master(fTs[i], fTs[j], b)
        Clhat_ij = nmt.compute_full_master(fs[i], fs[j], b)
        print('power spectrum taken')
        axes1[i,j-1].semilogx(ell_eff,
            ell_eff*ClhatT_ij[0]/np.sqrt(ClhatT_iis[i][0]*ClhatT_iis[j][0]))
        axes2[i,j-1].semilogx(ell_eff,
            ell_eff*Clhat_ij[0]/np.sqrt(Clhat_iis[i][0]*Clhat_iis[j][0]))
        axes3[i,j-1].semilogx(ell_eff,
            ell_eff*Clhat_ij[3]/np.sqrt(Clhat_iis[i][3]*Clhat_iis[j][3]))
        axes1[i,j-1].set_title(f'{diff_labs[i]}_{diff_labs[j]}')
        axes2[i,j-1].set_title(f'{diff_labs[i]}_{diff_labs[j]}')
        axes3[i,j-1].set_title(f'{diff_labs[i]}_{diff_labs[j]}')
ms = np.zeros((10, 3, hp.nside2npix(512)))
for i in range(10):
    ms[i] = hp.read_map(wmap_maps[i], field=(0,1,2))*1e3
    ms[i][0] -= dip_W*1e3
    #cg.plot(ms[i], min=-0.25, max=0.25)


diffs = np.zeros((4, 3, hp.nside2npix(512)))
diffs[0] = (ms[1] - 0.32*ms[0])/2
diffs[1] = (ms[3] - ms[2])/2
diffs[2] = (ms[5] - ms[4])/2
diffs[3] = (ms[9] - ms[8] - (ms[7] - ms[6]))/4

#for i in range(len(diffs)):
#    cg.plot(diffs[i], sig=1, fwhm=2*u.deg)
#plt.show()

fs = [nmt.NmtField(mask, diffs[i][1:]) for i in range(len(diffs))]
fTs = [nmt.NmtField(mask, [diffs[i][0]]) for i in range(len(diffs))]
Clhat_iis = []
ClhatT_iis = []
for i in range(len(diffs)):
    Clhat_iis.append(nmt.compute_full_master(fs[i], fs[i], b))
    ClhatT_iis.append(nmt.compute_full_master(fTs[i], fTs[i], b))
for i in range(len(diffs)):
    print(i)
    for j in range(1,i):
        print(i,j, 'off')
        #if j > 0: axes[i,j-1].axis('off')
    for j in range(i+1, len(diffs)):
        print(i,j)
        ClhatT_ij = nmt.compute_full_master(fTs[i], fTs[j], b)
        Clhat_ij = nmt.compute_full_master(fs[i], fs[j], b)
        print('power spectrum taken')
        axes1[i,j-1].semilogx(ell_eff,
            ell_eff*ClhatT_ij[0]/np.sqrt(ClhatT_iis[i][0]*ClhatT_iis[j][0]))
        axes2[i,j-1].semilogx(ell_eff,
            ell_eff*Clhat_ij[0]/np.sqrt(Clhat_iis[i][0]*Clhat_iis[j][0]))
        axes3[i,j-1].semilogx(ell_eff,
            ell_eff*Clhat_ij[3]/np.sqrt(Clhat_iis[i][3]*Clhat_iis[j][3]))
plt.savefig('wmap_cross_specs_BB.png', bbox_inches='tight')
plt.close()
plt.savefig('wmap_cross_specs_EE.png', bbox_inches='tight')
plt.close()
plt.savefig('wmap_cross_specs_TT.png', bbox_inches='tight')
plt.close()
plt.show()
