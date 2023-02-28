import matplotlib.pyplot as plt
import cosmoglobe as cg
import numpy as np
import healpy as hp

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
    width=8)
cg.plot(f'{DIR}/CG_cmb_IQU_n1024_v1.fits', comp='cmb', 
    width=8, remove_dip='auto')
cg.plot(dip, comp='cmb', min=-3400, max=3400, width=8)
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
