import numpy as np
import healpy as hp
import cosmoglobe as cg
import h5py
import matplotlib.pyplot as plt
#from setup_matplotlib import *


CGDIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_writeW2_230217'
data = h5py.File(f'{CGDIR}/tod_001304_samp000002.h5', 'r')
chain = cg.Chain(f'{CGDIR}/chain_c0001.h5')

d2 = h5py.File('/mn/stornext/d16/cmbco/ola/wmap/tods/hdf5/calibrated/wmap_W2_001304.h5', 'r')
sl = data['sl'][1]
s_tot = data['s_tot'][1]
tod = data['tod'][1]
s_sky = data['s_sky'][1]
n_corr = data['n_corr'][1]
gain = data['gain_000002'][()]
bpcorr = data['bpcorr'][1]
s_orb = data['s_orb'][1]

mask = data['flag'][0]

res = (tod - n_corr)/gain - s_tot

fig, axes = plt.subplots(sharex=True, nrows=7, figsize=(8, 6))

dt = 1.536/30
t = np.arange(0, dt*tod.size, dt)
inds = (t < 60*100)

axes[0].plot(t[inds], tod[inds]/gain, 'k.', ms=2, alpha=0.5)
axes[0].plot(t[inds], s_sky[inds], 'r')
axes[1].plot(t[inds], tod[inds]/gain - d2['001304/W214/tod'][inds], 'k.', ms=2,
    alpha=0.5)
axes[2].plot(t[inds], n_corr[inds]/gain, 'k')
axes[3].plot(t[inds], s_orb[inds], 'k')
axes[4].plot(t[inds], sl[inds], 'k')
axes[5].plot(t[inds], bpcorr[inds], 'k')

axes[6].plot(t[inds], res[inds], 'k.', ms=2, alpha=0.5)

axes[0].set_ylabel(r"$s_\mathrm{cal}\ \mathrm{[mK]}$")
axes[1].set_ylabel(r"$\Delta s_\mathrm{cal}\ \mathrm{[mK]}$")
axes[2].set_ylabel(r'$s_\mathrm{corr}$')
axes[3].set_ylabel(r'$s_\mathrm{orb}$')
axes[4].set_ylabel(r'$s_\mathrm{sl}$')
axes[5].set_ylabel(r'$s_\mathrm{leak}$')
axes[6].set_ylabel(r'$s_\mathrm{res}$')

#axes[0].set_ylim([-60, 60])
#axes[1].set_ylim([-60, 60])
#axes[2].set_ylim([-60, 60])
#axes[3].set_ylim([-60, 60])
#axes[4].set_ylim([-60, 60])
#axes[5].set_ylim([-60, 60])
#axes[6].set_ylim([-60, 60])

pid = 32

from cosmoglobe.tod_tools import TODLoader
comm_tod = TODLoader("/mn/stornext/d16/cmbco/ola/wmap/tods/hdf5/uncalibrated", "wmap")
comm_tod.init_file('K1', f'{pid:06}')
ztod = comm_tod.load_field(f'{pid:06}/K113/ztod')


baseline = chain.get('tod/023-WMAP_K/baseline')
data = h5py.File(f'{CGDIR}/tod_023-WMAP_K{pid:06}_samp000002.h5', 'r')
d2 = h5py.File(f'/mn/stornext/d16/cmbco/ola/wmap/tods/hdf5/calibrated/wmap_K1_{pid:06}.h5', 'r')
sl = data['sl'][0]
s_tot = data['s_tot'][0]
tod = data['tod'][0]
s_sky = data['s_sky'][0]
n_corr = data['n_corr'][0]
gain = data['gain_000001'][()]
bpcorr = data['bpcorr'][0]
s_orb = data['s_orb'][0]

mask = data['flag'][0]

d_calib = (tod - n_corr)/gain - s_tot + s_sky - bpcorr

res = (tod - n_corr)/gain - s_tot


fig, axes = plt.subplots(sharex=True, nrows=8, figsize=(8, 8))
fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)

dt = 1.536/12
t = np.arange(0, dt*tod.size, dt)
inds = (t < 60*10)

xi_n = chain.get(f'tod/023-WMAP_K/xi_n')
sigma = xi_n[0, 0, 0,pid-1]/gain

axes[0].plot(t[inds], ztod[inds], 'k.', ms=2, alpha=0.5)
axes[1].plot(t[inds], tod[inds]/gain, 'k.', ms=2, alpha=0.5)
axes[1].plot(t[inds], s_sky[inds], 'r')
#axes[2].plot(t[inds], d_calib[inds] - d2[f'{pid:06}/K113/tod'][inds], 'k.', ms=2,
axes[2].plot(t[inds], tod[inds]/gain - sl[inds] - d2[f'{pid:06}/K113/tod'][inds], 'k.', ms=2,
    alpha=0.5)
axes[3].plot(t[inds], n_corr[inds]/gain, 'k')
axes[4].plot(t[inds], s_orb[inds], 'k')
axes[5].plot(t[inds], sl[inds], 'k')
axes[6].plot(t[inds], bpcorr[inds], 'k')

axes[7].plot(t[inds], res[inds]/sigma, 'k.', ms=2, alpha=0.5)

axes[0].set_ylabel(r"$d_\mathrm{raw}\ \mathrm{[du]}$")
axes[1].set_ylabel(r"$d_\mathrm{cal}\ \mathrm{[mK]}$")
axes[2].set_ylabel(r"$\Delta s_\mathrm{cal}\ \mathrm{[mK]}$")
axes[3].set_ylabel(r'$s_\mathrm{corr}\ \mathrm{[mK]}$')
axes[4].set_ylabel(r'$s_\mathrm{orb}\ \mathrm{[mK]}$')
axes[5].set_ylabel(r'$s_\mathrm{sl}\ \mathrm{[mK]}$')
axes[6].set_ylabel(r'$s_\mathrm{leak}\ \mathrm{[mK]}$')
axes[7].set_ylabel(r'$s_\mathrm{res}\ [\sigma]$')

axes[0].set_ylim([-32137-40, -32137+40])
axes[1].set_ylim([-40, 40])
axes[2].set_ylim([-0.01, 0.075])
axes[3].set_ylim([-0.5, 0.5])
axes[4].set_ylim([-1, 1])
axes[5].set_ylim([-0.2, 0.2])
axes[6].set_ylim([-0.25, 0.25])
axes[7].set_ylim([-4, 4])

plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
for ax in axes:
  for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")
    ticklabel.set_verticalalignment("center")

axes[1].set_yticks([-20,0,20])
#axes[2].set_yticks([-0.05, 0, 0.05])
#axes[2].set_yticklabels([-0.05, 0, 0.05])
axes[2].set_yticks([0.0, 0.05])
axes[2].set_yticks([0.0,0.01, 0.02, 0.03,0.04, 0.05,0.06], minor=True)
axes[3].set_yticks([-0.3, 0, 0.3])
axes[3].set_yticklabels([-0.3, 0, 0.3])
axes[4].set_yticks([-0.6, 0, 0.6])
axes[4].set_yticklabels([-0.6, 0, 0.6])
axes[5].set_yticks([-0.1, 0, 0.1])
axes[5].set_yticklabels([-0.1, 0, 0.1])
axes[6].set_yticks([-0.15, 0, 0.15])
axes[6].set_yticklabels([-0.15, 0, 0.15])
axes[7].set_yticks([-3, 0, 3])

axes[7].set_xlabel(r'Time [seconds]')
axes[7].set_xticks(np.arange(0, 660, 60))




plt.savefig('K113_timestreams.pdf', bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
