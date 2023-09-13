import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
import cosmoglobe as cg

DIR = '/mn/stornext/d5/data/duncanwa/comm1'

bands = ['023', '030', '033', '041_Q1', '041_Q2', '044', '061_V1', '061_V2', '070',
    '094_W1', '094_W2', '094_W3', '094_W4']

band_labels = ['K', '30', 'Ka', 'Q1', 'Q2', '44', 'V1', 'V2', '70','W1', 'W2', 'W3', 'W4']


fnames_res = glob(f'{DIR}/cg_plots/res*.fits')
fnames_res.sort()
fnames_res = fnames_res[:-1]

fnames = glob(f'{DIR}/maps/CG*_sigma0*600arc_uK*.fits')
fnames.sort()

fnames[1], fnames[2] = fnames[2], fnames[1]

plt.figure(figsize=(12, 4))
n = 0
for i in range(len(fnames_res)):
    n += 1
    if band_labels[i] == 'W1':
        n += 1
    m  = hp.read_map(fnames_res[i], field=(0,1,2))
    if 'W' in band_labels[i]:
        cg.plot(m, sig=1, sub=(3,5, n), cbar=False,
                min=-10, max=10, llabel=band_labels[i])
    else:
        cg.plot(m, sig=1, sub=(3,5, n), cbar=False,
                min=-5, max=5, llabel=band_labels[i])
plt.tight_layout()
plt.savefig('../figures/comm1_res_Q.pdf', bbox_inches='tight')

plt.figure(figsize=(12, 4))
n = 0
for i in range(len(fnames_res)):
    n += 1
    if band_labels[i] == 'W1':
        n += 1
    m  = hp.read_map(fnames_res[i], field=(0,1,2))
    if 'W' in band_labels[i]:
        cg.plot(m, sig=2, sub=(3,5, n), cbar=False,
                min=-10, max=10, llabel=band_labels[i])
    else:
        cg.plot(m, sig=2, sub=(3,5, n), cbar=False,
                min=-5, max=5, llabel=band_labels[i])
plt.tight_layout()
plt.savefig('../figures/comm1_res_U.pdf', bbox_inches='tight')
