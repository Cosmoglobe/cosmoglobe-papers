import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import astropy.units as u
import healpy as hp

import cosmoglobe as cg
CGDIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203'

fnames = glob(f'{CGDIR}/tod_023-WMAP_K_*_c0001_k000121.fits')
fnames.sort()

width=3

comps = ['ncorr', 'orb', 'bpcorr', 'sl', 'res']
m = hp.read_map(f'{CGDIR}/tod_023-WMAP_K_{comps[0]}_c0001_k000121.fits',
    field=(0,1,2))
m = hp.smoothing(m, fwhm=2*np.pi/180)
cg.plot(m, min=-1, max=1, scale=1e3, sig=0, cbar=True, width=width,
    extend='both', rlabel=r'T_\mathrm{ncorr}^{2^\circ}')
plt.savefig('K_ncorr_I.pdf', bbox_inches='tight')
cg.plot(m, min=-1, max=1, scale=1e3, sig=1, cbar=True, width=width,
    extend='both', rlabel=r'Q_\mathrm{ncorr}^{2^\circ}')
plt.savefig('K_ncorr_Q.pdf', bbox_inches='tight')
cg.plot(m, min=-1, max=1, scale=1e3, sig=2, cbar=True, width=width,
    extend='both', rlabel=r'U_\mathrm{ncorr}^{2^\circ}')
plt.savefig('K_ncorr_U.pdf', bbox_inches='tight')

m = hp.read_map(f'{CGDIR}/tod_023-WMAP_K_{comps[1]}_c0001_k000121.fits',
    field=(0,1,2))
cg.plot(m, min=-250, max=250,  scale=1e3, sig=0, rlabel='T_\mathrm{orb}',
    width=width, extend='both')
plt.savefig('K_orb_I.pdf', bbox_inches='tight')
cg.plot(m, min=-2.5, max=2.5, scale=1e3, sig=1, rlabel='Q_\mathrm{orb}',
    width=width, extend='both')
plt.savefig('K_orb_Q.pdf', bbox_inches='tight')
cg.plot(m, min=-2.5, max=2.5, scale=1e3, sig=2, rlabel='U_\mathrm{orb}',
    width=width, extend='both')
plt.savefig('K_orb_U.pdf', bbox_inches='tight')

m = hp.read_map(f'{CGDIR}/tod_023-WMAP_K_{comps[2]}_c0001_k000121.fits',
    field=(0,1,2))
cg.plot(m, min=-1, max=1, scale=1e3, sig=0, width=width, extend='both',
    rlabel=r'T_\mathrm{leak}')
plt.savefig('K_leak_I.pdf', bbox_inches='tight')
cg.plot(m, min=-25, max=25, scale=1e3, sig=1, width=width, extend='both',
    rlabel=r'Q_\mathrm{leak}')
plt.savefig('K_leak_Q.pdf', bbox_inches='tight')
cg.plot(m, min=-25, max=25, scale=1e3, sig=2, width=width, extend='both',
    rlabel=r'U_\mathrm{leak}')
plt.savefig('K_leak_U.pdf', bbox_inches='tight')


m = hp.read_map(f'{CGDIR}/tod_023-WMAP_K_{comps[3]}_c0001_k000121.fits',
    field=(0,1,2))
cg.plot(m, scale=1e3, sig=0, min=-50, max=50, width=width,
    rlabel=r'T_\mathrm{sl}', extend='both')
plt.savefig('K_sl_I.pdf', bbox_inches='tight')
cg.plot(m, scale=1e3, sig=1, min=-1, max=1, width=width,
    rlabel=r'Q_\mathrm{sl}', extend='both')
plt.savefig('K_sl_Q.pdf', bbox_inches='tight')
cg.plot(m, scale=1e3, sig=2, min=-1, max=1, width=width,
    rlabel=r'U_\mathrm{sl}', extend='both')
plt.savefig('K_sl_U.pdf', bbox_inches='tight')


m = hp.read_map(f'{CGDIR}/tod_023-WMAP_K_{comps[4]}_c0001_k000121.fits',
    field=(0,1,2))
m = hp.smoothing(m, fwhm=2*np.pi/180)
cg.plot(m, scale=1e3, sig=0, min=-10, max=10, width=width, extend='both',
    rlabel=r'T_\mathrm{res}^{2^\circ}', unit=r'$\mathrm{\mu K}$')
plt.savefig('K_res_I.pdf', bbox_inches='tight')
cg.plot(m, scale=1e3, sig=1, min=-10, max=10, width=width, extend='both',
    rlabel=r'Q_\mathrm{res}^{2^\circ}', unit=r'$\mathrm{\mu K}$')
plt.savefig('K_res_Q.pdf', bbox_inches='tight')
cg.plot(m, scale=1e3, sig=2, min=-10, max=10, width=width, extend='both',
    rlabel=r'U_\mathrm{res}^{2^\circ}', unit=r'$\mathrm{\mu K}$')
plt.savefig('K_res_U.pdf', bbox_inches='tight')
