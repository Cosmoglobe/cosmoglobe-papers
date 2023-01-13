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
    "/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/wmap_kq75_TQU_mask_r9.fits"
)


# bin4 = nmt.NmtBin.from_edges(l_ini, l_end)
logbins = np.geomspace(10, 3 * 512 - 1, 50).astype("int")
l_ini = np.concatenate((np.arange(2, 10), logbins[:-1]))
l_end = np.concatenate((np.arange(2, 10) + 1, logbins[1:]))
# b = nmt.NmtBin.from_nside_linear(nside, 1)
#b = nmt.NmtBin.from_edges(l_ini, l_end)
#ell_eff = b.get_effective_ells()

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
CG_DIR = "/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_857_KKaQVW_a_230110"

cg_maps = [
    f"{CG_DIR}/tod_023-WMAP_K_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_030-WMAP_Ka_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_040-WMAP_Q1_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_040-WMAP_Q2_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_060-WMAP_V1_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_060-WMAP_V2_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_090-WMAP_W1_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_090-WMAP_W2_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_090-WMAP_W3_map_c0001_k000008.fits",
    f"{CG_DIR}/tod_090-WMAP_W4_map_c0001_k000008.fits",
]

bands = ["K", "Ka", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]

x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip_W = -0.233 * x - 2.222 * y + 2.504 * z

ms = np.zeros((10, 3, np.nside2npix(512)))
for i in range(10):
    ms[i] = hp.read_map(cg_maps[i], field=(0,1,2))
    ms[i][0] -= dip_W
    cg.plot(ms[i])
    plt.show()


n = 0
for i in range(len(cg_maps)):
    print(i)
    m_WMAP = hp.read_map(wmap_maps[i]) * 1e3
    m_CG = hp.read_map(cg_maps[i]) * 1e3
    m_CG -= dip_W * 1e3
    mono_CG = hp.fit_monopole(m_CG, gal_cut=50)
    mono_WM = hp.fit_monopole(m_WMAP, gal_cut=50)
    m_WMAP = m_WMAP - mono_WM + mono_CG
    f_WMAP = nmt.NmtField(mask, [m_WMAP])
    f_CG = nmt.NmtField(mask, [m_CG])
    f_diff = nmt.NmtField(mask, [m_CG - m_WMAP])
    Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)[0]
    Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)[0]

n = 0
for i in range(len(cg_maps)):
    m_WMAP = hp.read_map(wmap_maps[i], field=(1, 2)) * 1e3
    m_CG = hp.read_map(cg_maps[i], field=(1, 2)) * 1e3
    f_WMAP = nmt.NmtField(mask, m_WMAP)
    f_CG = nmt.NmtField(mask, m_CG)
    Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)
    Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)

