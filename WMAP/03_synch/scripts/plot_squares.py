import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import cosmoglobe as cg
import cmasher as cm
CG_DIR = '/mn/stornext/d5/data/duncanwa/WMAP'
DIR_CG = '/mn/stornext/d16/cmbco/cg/v1'

dpi=150
fontsize = {'llabel':9, 'rlabel':9, 'title':8}


synch_CG, h = hp.read_map(f'{DIR_CG}/CG_synch_IQU_n1024_v1.fits',
    field=(0,1,2,7,8), h=True)
DIR_BP = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'

d_K = hp.read_map(f'{DIR_CG}/CG_023-WMAP_K_IQU_n0512_v1.fits',
    field=(0,1,2,4,5,6))
d_30 = hp.read_map(f'{DIR_CG}/CG_030_IQU_n0512_v1.fits', field=(0,1,2,4,5,6))

# Smooth *to* 72 arcmin
scal = 1e3*0.47
data_dir='/mn/stornext/d5/data/duncanwa/WMAP/data'

bl_W = [f"{data_dir}/WMAP9_K1_beam_ext.fits"]
bl_P = [f"{data_dir}/Bl_TEB_npipe6v19_30GHzx30GHz.fits"]

bl_list = bl_W + bl_P
lmax = 800
bls = [hp.read_cl(bl_f) for bl_f in bl_list]


fwhm = 72

bl_G = hp.gauss_beam(fwhm*np.pi/180/60, lmax=lmax)


#fwhm = np.sqrt((72/60)**2 - 0.88**2)*np.pi/180
#d_K2 = hp.smoothing(d_K[:3], fwhm=fwhm)*scal
alm = hp.map2alm(d_K[:3]*scal, lmax=lmax)
d_K2 = hp.alm2map(np.array([
    hp.almxfl(almi, bl_G/bls[0][:lmax+1]) for almi in alm]), 512)
fwhm = np.sqrt((72/60)**2 - (30/60)**2)*np.pi/180
#d_302 = hp.smoothing(d_30[:3], fwhm=fwhm)
alm = hp.map2alm(d_30[:3], lmax=lmax)
d_302 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[1][i][:lmax+1]) for i in range(3)]), 512)

cg_P = np.hypot(synch_CG[1], synch_CG[2])


K_P = np.hypot(d_K2[1], d_K2[2])
P30_P = np.hypot(d_302[1], d_302[2])


plt.figure()
hp.gnomview(P30_P, rot=(-100, 53, 0), min=0, max=10,
        sub=(1,3,1), cmap=cm.swamp, reso = 10)
plt.title('30 GHz')
hp.gnomview(K_P, rot=(-100, 53, 0), min=0, max=10,
        sub=(1,3,2), cmap=cm.swamp, reso = 10)
plt.title('K')
hp.gnomview(cg_P, rot=(-100, 53, 0), min=0, max=10,
        sub=(1,3,3), cmap=cm.swamp, reso = 10)
plt.title('Synch map')
plt.tight_layout()

plt.figure()
cg.gnom(P30_P, lat=53, lon=-100, min=0, max=10, sub=(1,3,1), cmap='swamp',
        cbar=True, llabel=r'\mathrm{30\ GHz}', ticks=[0,5,10])
cg.gnom(K_P, lat=53, lon=-100, min=0, max=10, sub=(1,3,2), cmap='swamp',
        cbar=True, llabel=r'\mathit{K}', ticks=[0,5,10], unit=r'\mathrm{\mu K}')
plt.title('Polarization amplitude at $(l,b)=(-100^\circ,53^\circ)$')
cg.gnom(cg_P, lat=53, lon=-100, min=0, max=10, sub=(1,3,3), cmap='swamp',
        cbar=True, llabel=r'\mathrm{CG\ Synch}', ticks=[0,5,10])
plt.tight_layout()

plt.figure()
cg.gnom(K_P, lat=53, lon=-100, min=0, max=10, cmap='swamp',
        cbar=False, llabel=r'\mathit{K}', ticks=[0,5,10], unit=r'\mathrm{\mu K}',
        figsize=(2, 8/3), fontsize=fontsize)
plt.title('$(l,b)=(-100^\circ,53^\circ)$')
plt.savefig('../figures/K_square.pdf', bbox_inches='tight', dpi=dpi)

plt.figure()
cg.gnom(P30_P, lat=53, lon=-100, min=0, max=10, cmap='swamp',
        cbar=False, llabel=r'\mathrm{30}', ticks=[0,5,10], figsize=(2,8/3),
        fontsize=fontsize)
plt.savefig('../figures/30_square.pdf', bbox_inches='tight', dpi=dpi)
plt.figure()
cg.gnom(cg_P, lat=53, lon=-100, min=0, max=10, sub=(1,3,3), cmap='swamp',
        cbar=True, llabel=r'\mathrm{Cosmoglobe}', ticks=[0,5,10],
        unit=r'$\mathrm{\mu K}$', figsize=(2, 8/3), fontsize=fontsize)
plt.savefig('../figures/CG_square.pdf', bbox_inches='tight', dpi=dpi)
#plt.tight_layout()
#plt.show()
