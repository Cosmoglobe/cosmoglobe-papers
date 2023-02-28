import matplotlib.pyplot as plt
import cosmoglobe as cg
import numpy as np
import healpy as hp

import astropy.units as u

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/v1'

A = 3366.168
lon = 264.081 * np.pi/180 - np.pi
lat =  48.274 * np.pi/180
nside = 1024

d_x = A * np.cos(lon) * np.sin(np.pi/2-lat)
d_y = A * np.sin(lon) * np.sin(np.pi/2-lat)
d_z = A * np.cos(np.pi/2-lon)

x,y,z = hp.pix2vec(nside, np.arange(12*nside**2))

dip = d_x*x + d_y*y + d_z*z


m = hp.read_map(f'{DIR}/CG_cmb_IQU_n1024_v1.fits')
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', comp='cmb', min=-3400, max=3400,
    width=8, xsize=800)
plt.savefig('cmb_I_dipole.pdf', bbox_inches='tight')


cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', comp='cmb', 
    width=3, remove_dip='auto',
    rlabel=r'\langle A_\mathrm{cmb}\rangle')
plt.savefig('cmb_I_nodipole.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', comp='cmb', 
    width=3, sig=1, fwhm=np.sqrt(1-(14/60)**2)*u.deg,
    min=-10, max=10,
    rlabel=r'\langle A_\mathrm{cmb}\rangle')
plt.savefig('cmb_Q.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', comp='cmb', 
    width=3, sig=2, fwhm=np.sqrt(1-(14/60)**2)*u.deg,
    min=-10, max=10,
    rlabel=r'\langle A_\mathrm{cmb}\rangle')
plt.savefig('cmb_U.pdf', bbox_inches='tight')


cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', sig=3, min=0, max=30,
    cmap='binary_r', unit=r'\mathrm{\mu K_{CMB}}', llabel='T',
    rlabel=r'\sigma_\mathrm{cmb}', width=3, extend='both')
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', sig=4, min=0, max=30,
    cmap='binary_r', unit=r'\mathrm{\mu K_{CMB}}', llabel='Q',
    rlabel=r'\sigma_\mathrm{cmb}', width=3, extend='both')
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', sig=5, min=0, max=30,
    cmap='binary_r', unit=r'\mathrm{\mu K_{CMB}}', llabel='U',
    rlabel=r'\sigma_\mathrm{cmb}', width=3, extend='both')

plt.show()

cg.plot(f'{DIR}/CG_synch_IQU_n1024_v1.fits', comp='synch', llabel='T',
    rlabel=r'\langle A_\mathrm s\rangle', scale=1e-6,
    unit=r'$\mathrm{K_{RJ}}$', min=10, max=200, ticks=[10,100,200], width=4)
plt.savefig('synch_I.pdf', bbox_inches='tight')

cg.plot(f'{DIR}/CG_synch_IQU_n1024_v1.fits', comp='synch', sig=1, width=4,
    rlabel=r'\langle A_\mathrm s\rangle', min=-30, max=30, norm='linear')
plt.savefig('synch_Q.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_synch_IQU_n1024_v1.fits', comp='synch', sig=2, width=4,
    rlabel=r'\langle A_\mathrm s\rangle', norm='linear', min=-30, max=30)
plt.savefig('synch_U.pdf', bbox_inches='tight')

cg.plot(f'{DIR}/CG_freefree_I_n1024_v1.fits', comp='ff', width=4,
    rlabel=r'\langle A_\mathrm{ff}\rangle', llabel='T')
plt.savefig('ff_I.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_ame_I_n1024_v1.fits', comp='ame', width=4, llabel='T',
    rlabel=r'\langle A_\mathrm{ame}\rangle')
plt.savefig('ame_I.pdf', bbox_inches='tight')


cg.plot(f'{DIR}/CG_dust_IQU_n1024_v1.fits', comp='dust', min=3, max=300,
    ticks=[3, 30, 300], width=4,llabel='T',
    rlabel=r'\langle A_\mathrm{d}\rangle')
plt.savefig('dust_I.pdf', bbox_inches='tight')

cg.plot(f'{DIR}/CG_dust_IQU_n1024_v1.fits', comp='dust', min=-5, max=5,
    norm='linear', width=4,llabel='Q',sig=1,
    rlabel=r'\langle A_\mathrm{d}\rangle')
plt.savefig('dust_Q.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_dust_IQU_n1024_v1.fits', comp='dust', min=-5, max=5,
    norm='linear', width=4,llabel='U',sig=2,
    rlabel=r'\langle A_\mathrm{d}\rangle')
plt.savefig('dust_U.pdf', bbox_inches='tight')
