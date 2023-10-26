# Configure Matplotlib options
import healpy as hp
from matplotlib.ticker import MaxNLocator
import numpy as np
import sys
import matplotlib.pyplot as plt

import pymaster as nmt
import healpy as hp

from glob import glob
from tqdm import tqdm

nside = 1024

# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map(
        "COM_Mask_CMB-common-Mask-Pol_2048_R3.00.fits"
        )
mask = hp.ud_grade(mask, nside)

# check that I'm using the same mask for the power spectra
# Compute the 2015 synch as well, get the half-mission maps, etc.


l_ini = np.concatenate((
  np.array([2,4,12]),
  np.arange(20, 520, 20),
  np.array([550])))
l_end = np.concatenate((
  np.array([3,11,19]),
  np.arange(40, 520, 20) - 1,
  np.array([549, 599])))

assert len(l_ini) == len(l_end)
for i in range(len(l_ini)):
    print(l_ini[i], l_end[i])

b = nmt.NmtBin.from_edges(l_ini, l_end)
ell_eff = b.get_effective_ells()
Z = ell_eff*(ell_eff+1)/(2*np.pi)

fnames_1 = glob('/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_HM1/synch_c0001_k??????.fits')
fnames_2 = glob('/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_HM2/synch_c0001_k??????.fits')

fnames_1.sort()
fnames_2.sort()

fnames_1 = fnames_1[2:]
fnames_2 = fnames_2[2:]

L = min(len(fnames_1), len(fnames_2))
EEs = np.zeros((L, len(ell_eff)))
BBs = np.zeros((L, len(ell_eff)))

#fnames_1 = glob('*hm1.fits')
#fnames_2 = glob('*hm2.fits')
#L = 1

for i in tqdm(range(L)):
    synch_cg11 = hp.read_map(fnames_1[i],
        field=(0,1,2))
    synch_cg12 = hp.read_map(fnames_2[i],
        field=(0,1,2))
    #synch_cg11 = hp.ud_grade(hp.read_map(fnames_1[i],
    #    field=(0,1)), nside)
    #synch_cg12 = hp.ud_grade(hp.read_map(fnames_2[i],
    #    field=(0,1)), nside)
    
    nmt_field1 = nmt.NmtField(mask, synch_cg11[1:])
    nmt_field2 = nmt.NmtField(mask, synch_cg12[1:])
    #nmt_field1 = nmt.NmtField(mask, synch_cg11)
    #nmt_field2 = nmt.NmtField(mask, synch_cg12)
    
    Clhat_cg1 = nmt.compute_full_master(nmt_field1, nmt_field2, b)
    EEs[i] = Clhat_cg1[0]
    BBs[i] = Clhat_cg1[3]
    plt.figure(1)
    plt.loglog(ell_eff, Z*EEs[i], color='k', alpha=0.5)
    plt.figure(2)
    plt.loglog(ell_eff, Z*BBs[i], color='k', alpha=0.5)
    if ((i+1) % 5) == 0:
        plt.show() 

plt.show()
np.savetxt('ee_spectra.txt', EEs)
np.savetxt('bb_spectra.txt', BBs)
#np.savetxt('ee_spectra_2018.txt', EEs)
#np.savetxt('bb_spectra_2018.txt', BBs)
