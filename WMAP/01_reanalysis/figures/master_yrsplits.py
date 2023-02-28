import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import os.path

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


b = nmt.NmtBin.from_nside_linear(nside, 1)

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
CG_DIR = "/mn/stornext/d16/cmbco/bp/dwatts/WMAP/"
CG_DIR = "/mn/stornext/d5/data/duncanwa/WMAP"

cg_maps = [
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_023-WMAP_K_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_030-WMAP_Ka_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_040-WMAP_Q1_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_040-WMAP_Q2_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_060-WMAP_V1_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_060-WMAP_V2_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_090-WMAP_W1_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_090-WMAP_W2_map_c0001_k000001.fits",
    f"{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_090-WMAP_W3_map_c0001_k000001.fits",
]
# f'{CG_DIR}/chains_CG_LFI_KKaQVW_221130/tod_090-WMAP_W4_map_c0001_k000002.fits']


DIRs = [
    "/mn/stornext/d5/data/duncanwa/WMAP/chains_baselinetest",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps",
]
labels = ["CG", "W9"]

bands = [
    "023-WMAP_K",
    "030-WMAP_Ka",
    "040-WMAP_Q1",
    "040-WMAP_Q2",
    "060-WMAP_V1",
    "060-WMAP_V2",
    "090-WMAP_W1",
    "090-WMAP_W2",
    "090-WMAP_W3",
    "090-WMAP_W4",
]

for band in bands:
    fnames = glob(f"{DIRs[0]}/tod_{band}_map_yr000?_c0001_k000001.fits")
    fnames.sort()
    print(fnames)
    if len(fnames) != 9:
        continue
    if os.path.exists("fnames_all_{band}_{label}.npy"):
        continue
    f2s = []
    for f in fnames:
        print(f)
        f2s.append(nmt.NmtField(mask, hp.read_map(f, field=(1, 2))))
    Cls = []
    for i in range(len(f2s)):
        for j in range(i + 1, len(f2s)):
            print(i, j)
            cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
            Cls.append(cl_22)

    Cls = np.array(Cls)

    np.save(f"fnames_all_{band}_{labels[0]}", Cls)

bands = ["K1", "Ka1", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]
for band in bands:
    fnames = glob(f"{DIRs[1]}/wmap_iqusmap_r9_yr?_{band}_v5.fits")
    fnames.sort()
    if len(fnames) != 9:
        continue
    if os.path.exists("fnames_all_{band}_{label}.npy"):
        continue
    f2s = []
    for f in fnames:
        print(f)
        f2s.append(nmt.NmtField(mask, hp.read_map(f, field=(1, 2))))
    Cls = []
    for i in range(len(f2s)):
        for j in range(i + 1, len(f2s)):
            print(i, j)
            cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
            Cls.append(cl_22)

    Cls = np.array(Cls)

    np.save(f"fnames_all_{band}_{labels[1]}", Cls)


fnames1 = glob(f"{DIR}/wmap_iqusmap_r9_yr?_W1_v5.fits")
fnames2 = glob(f"{DIR}/wmap_iqusmap_r9_yr?_W2_v5.fits")
fnames3 = glob(f"{DIR}/wmap_iqusmap_r9_yr?_W3_v5.fits")
fnames4 = glob(f"{DIR}/wmap_iqusmap_r9_yr?_W4_v5.fits")

fnames = fnames1 + fnames2 + fnames3
f2s = []
for f in fnames:
    print(f)
    f2s.append(nmt.NmtField(mask, hp.read_map(f, field=(1, 2))))

Cls = []

for i in range(len(f2s)):
    for j in range(i + 1, len(f2s)):
        print(i, j)
        cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
        Cls.append(cl_22)

Cls = np.array(Cls)

np.save("fnames_all", Cls)

W1 = hp.read_map(f"{DIR}/wmap_iqusmap_r9_9yr_W1_v5.fits", field=(1, 2))
W2 = hp.read_map(f"{DIR}/wmap_iqusmap_r9_9yr_W2_v5.fits", field=(1, 2))
W3 = hp.read_map(f"{DIR}/wmap_iqusmap_r9_9yr_W3_v5.fits", field=(1, 2))
W4 = hp.read_map(f"{DIR}/wmap_iqusmap_r9_9yr_W4_v5.fits", field=(1, 2))

Ws = [W1, W2, W3, W4]
f2s = []
for f in Ws:
    f2s.append(nmt.NmtField(mask, f))

crosses = []
for i in range(len(f2s)):
    for j in range(i, len(f2s)):
        cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
        crosses.append(cl_22)

crosses = np.array(crosses)
np.save("crosses_wmap9", crosses)

DIR = "/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_CG_LFI_KaQVW_221124"
W1 = hp.read_map(f"{DIR}/tod_090-WMAP_W1_map_c0001_k000010.fits", field=(1, 2))
W2 = hp.read_map(f"{DIR}/tod_090-WMAP_W2_map_c0001_k000010.fits", field=(1, 2))
W3 = hp.read_map(f"{DIR}/tod_090-WMAP_W3_map_c0001_k000010.fits", field=(1, 2))
W4 = hp.read_map(f"{DIR}/tod_090-WMAP_W4_map_c0001_k000010.fits", field=(1, 2))

Ws = [W1, W2, W3, W4]
f2s = []
for f in Ws:
    f2s.append(nmt.NmtField(mask, f))

crosses = []
for i in range(len(f2s)):
    for j in range(i, len(f2s)):
        cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
        crosses.append(cl_22)

crosses = np.array(crosses)
np.save("crosses_CG", crosses)

n1 = W1 - (W2 + W3 + W4) / 3
n2 = W2 - (W1 + W3 + W4) / 3
n3 = W3 - (W1 + W2 + W4) / 3
n4 = W4 - (W1 + W2 + W3) / 3
n5 = (W1 + W2) / 2 - (W3 + W4) / 2
n6 = (W1 + W3) / 2 - (W2 + W4) / 2

nulls = [n1, n2, n3, n4, n5, n6]
Cls = []
for null in nulls:
    f2 = nmt.NmtField(mask, null)
    cl_22 = nmt.compute_full_master(f2, f2, b)
    Cls.append(cl_22)
Cls = np.array(Cls)
np.save("null_spectra_CG", Cls)
