import numpy as np
import matplotlib.pyplot as plt
import cosmoglobe as cg
import healpy as hp
from glob import glob
import astropy.units as u
import h5py

fontsize = {'llabel':8, 'rlabel':8, 'llabel_align':'right'}

# Need to smooth all to common resolution...
data_dir='/mn/stornext/d5/data/duncanwa/WMAP/data'

bl_W = [f"{data_dir}/WMAP9_K1_beam_ext.fits", f"{data_dir}/WMAP9_Ka1_beam_ext.fits"]
bl_P = [f"{data_dir}/Bl_TEB_npipe6v19_30GHzx30GHz.fits", f"{data_dir}/Bl_TEB_npipe6v19_30GHzx30GHz.fits"]

bl_list = bl_W + bl_P
lmax = 800
bls = [hp.read_cl(bl_f) for bl_f in bl_list]


fwhm = 300

bl_G = hp.gauss_beam(fwhm*np.pi/180/60, lmax=lmax)


CG_DIR = '/mn/stornext/d5/data/duncanwa/WMAP'

DIR_CG = '/mn/stornext/d16/cmbco/cg/v1'

synch_CG, h = hp.read_map(f'{DIR_CG}/CG_synch_IQU_n1024_v1.fits',
    field=(0,1,2), h=True)

synch_CG = hp.smoothing(synch_CG, fwhm=np.sqrt(fwhm**2-72**2)*np.pi/180/60)

synch_CG = hp.ud_grade(synch_CG, 512)

model = cg.sky_model_from_chain(f'{DIR_CG}/CG_c0001_v1.h5', nside=512,
        components=['synch', 'dust'])


CG_K = hp.read_map(f'{DIR_CG}/CG_023-WMAP_K_IQU_n0512_v1.fits',
        field=(0,1,2))*1e3
CG_Ka = hp.read_map(f'{DIR_CG}/CG_030-WMAP_Ka_IQU_n0512_v1.fits',
        field=(0,1,2))*1e3


WDIR = '/mn/stornext/d16/cmbco/ola/wmap/freq_maps'
WMAP_K = hp.read_map(f'{WDIR}/wmap_band_iqusmap_r9_9yr_K_v5.fits',
        field=(0,1,2))*1e3
WMAP_Ka = hp.read_map(f'{WDIR}/wmap_band_iqusmap_r9_9yr_Ka_v5.fits',
        field=(0,1,2))*1e3

PR3 = '/mn/stornext/d16/cmbco/ola/planck_products/dr3'
PR4 = '/mn/stornext/d16/cmbco/ola/planck_products/dr4'

PR3_30 = hp.read_map(f'{PR3}/LFI_SkyMap_030-BPassCorrected-field-IQU_1024_R3.00_full.fits',
        field=(0,1,2))*1e6
PR3_44 = hp.read_map(f'{PR3}/LFI_SkyMap_044-BPassCorrected-field-IQU_1024_R3.00_full.fits',
        field=(0,1,2))*1e6
PR4_30 = hp.read_map(f'{PR4}/LFI_SkyMap_030-BPassCorrected-field-IQU_1024_R4.00_full.fits',
        field=(0,1,2))*1e6
PR4_44 = hp.read_map(f'{PR4}/LFI_SkyMap_044-BPassCorrected-field-IQU_1024_R4.00_full.fits',
        field=(0,1,2))*1e6


CG_30 = hp.read_map(f'{DIR_CG}/CG_030_IQU_n0512_v1.fits', field=(0,1,2))
CG_44 = hp.read_map(f'{DIR_CG}/CG_044_IQU_n0512_v1.fits', field=(0,1,2))

DIR_BP = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'
BP_30 = hp.read_map(f'{DIR_BP}/BP_030_IQU_n0512_v2.fits', field=(0,1,2))
BP_44 = hp.read_map(f'{DIR_BP}/BP_044_IQU_n0512_v2.fits', field=(0,1,2))

alm = hp.map2alm(CG_K, lmax=lmax)
CG_K = hp.alm2map(np.array([
    hp.almxfl(almi, bl_G/bls[0][:lmax+1]) for almi in alm]), 512)
alm = hp.map2alm(CG_Ka, lmax=lmax)
CG_Ka = hp.alm2map(np.array([
    hp.almxfl(almi, bl_G/bls[1][:lmax+1]) for almi in alm]), 512)
alm = hp.map2alm(CG_30, lmax=lmax)
CG_30 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[2][i][:lmax+1]) for i in range(3)]), 512)
alm = hp.map2alm(CG_44, lmax=lmax)
CG_44 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[3][i][:lmax+1]) for i in range(3)]), 512)


alm = hp.map2alm(BP_30, lmax=lmax)
BP_30 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[2][i][:lmax+1]) for i in range(3)]), 512)
alm = hp.map2alm(BP_44, lmax=lmax)
BP_44 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[3][i][:lmax+1]) for i in range(3)]), 512)



alm = hp.map2alm(PR3_30, lmax=lmax)
PR3_30 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[2][i][:lmax+1]) for i in range(3)]), 512)
alm = hp.map2alm(PR3_44, lmax=lmax)
PR3_44 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[3][i][:lmax+1]) for i in range(3)]), 512)
alm = hp.map2alm(PR4_30, lmax=lmax)
PR4_30 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[2][i][:lmax+1]) for i in range(3)]), 512)
alm = hp.map2alm(PR4_44, lmax=lmax)
PR4_44 = hp.alm2map(np.array([
    hp.almxfl(alm[i], bl_G/bls[3][i][:lmax+1]) for i in range(3)]), 512)

alm = hp.map2alm(WMAP_K, lmax=lmax)
WMAP_K = hp.alm2map(np.array([
    hp.almxfl(almi, bl_G/bls[0][:lmax+1]) for almi in alm]), 512)
alm = hp.map2alm(WMAP_Ka, lmax=lmax)
WMAP_Ka = hp.alm2map(np.array([
    hp.almxfl(almi, bl_G/bls[1][:lmax+1]) for almi in alm]), 512)


#cg.plot(CG_K , sig=1, min=-25, max=25)
#cg.plot(CG_30 , sig=1, min=-25, max=25)
#cg.plot(CG_Ka , sig=1, min=-25, max=25)
#cg.plot(CG_44 , sig=1, min=-25, max=25)

data = h5py.File(f'{data_dir}/WMAP_instrument_v14.h5')
nu = data['023-WMAP_K/bandpassx'][()]
weights = data['023-WMAP_K/bandpass'][()]
nu *= u.GHz
weights *= u.Unit("K_RJ")
CG_model_K = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value

nu = data['030-WMAP_Ka/bandpassx'][()]
weights = data['030-WMAP_Ka/bandpass'][()]
nu *= u.GHz
weights *= u.Unit("K_RJ")
CG_model_Ka = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value


data = h5py.File(f'{data_dir}/LFI_instrument_v8.h5')
nu = data['030/bandpassx'] * u.GHz
weights = data['030/bandpass'] * u.Unit("K_RJ")
CG_model_30 = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value

nu = data['044/bandpassx'] * u.GHz
weights = data['044/bandpass'] * u.Unit("K_RJ")
CG_model_44 = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value


#plt.figure(figsize=(8, 8*2/3))
cg.plot(WMAP_K - CG_model_K, sig=1, min=-5, max=5, cbar=False,   sub=(5, 4, 1),
    llabel=r'\mathit{WMAP}\ K', fontsize=fontsize, rlabel='Q')
cg.plot(WMAP_K - CG_model_K, sig=2, min=-5, max=5, cbar=False,   sub=(5, 4, 2),
        fontsize=fontsize, rlabel='U')
cg.plot(CG_K - CG_model_K, sig=1, min=-5, max=5, cbar=False,     sub=(5, 4, 5),
        llabel=r'\mathrm{CG}\ K', fontsize=fontsize)
cg.plot(CG_K - CG_model_K, sig=2, min=-5, max=5, cbar=False,     sub=(5, 4, 6),
        fontsize=fontsize)

cg.plot(WMAP_Ka - CG_model_Ka, sig=1, min=-5, max=5, cbar=False, sub=(5, 4, 9),
    llabel=r'\mathit{WMAP}\ \mathit{Ka}', fontsize=fontsize)
cg.plot(WMAP_Ka - CG_model_Ka, sig=2, min=-5, max=5, cbar=False, sub=(5, 4, 10),
        fontsize=fontsize)
cg.plot(CG_Ka - CG_model_Ka, sig=1, min=-5, max=5, cbar=False,   sub=(5, 4, 13),
        llabel=r'\mathrm{CG}\ \mathit{Ka}', fontsize=fontsize)
cg.plot(CG_Ka - CG_model_Ka, sig=2, min=-5, max=5, cbar=False,   sub=(5, 4, 14),
        fontsize=fontsize)


cg.plot(PR3_30 - CG_model_30, sig=1, min=-5, max=5, cbar=False,  sub=(5, 4, 3),
    llabel=r'\mathrm{PR3}\ 30', fontsize=fontsize, rlabel='Q')
cg.plot(PR3_30 - CG_model_30, sig=2, min=-5, max=5, cbar=False,  sub=(5, 4, 4),
        fontsize=fontsize, rlabel='U')
cg.plot(CG_30 - CG_model_30, sig=1, min=-5, max=5, cbar=False,   sub=(5, 4, 7),
    llabel=r'\mathrm{CG}\ 30', fontsize=fontsize)
cg.plot(CG_30 - CG_model_30, sig=2, min=-5, max=5, cbar=False,   sub=(5, 4, 8),
        fontsize=fontsize)

cg.plot(PR3_44 - CG_model_44, sig=1, min=-5, max=5, cbar=False,  sub=(5, 4, 11),
    llabel=r'\mathrm{PR3}\ 44', fontsize=fontsize)
cg.plot(PR3_44 - CG_model_44, sig=2, min=-5, max=5, cbar=False,  sub=(5, 4, 12),
        fontsize=fontsize)
cg.plot(CG_44 - CG_model_44, sig=1, min=-5, max=5, cbar=False,   sub=(5, 4, 15),
    llabel=r'\mathrm{CG}\ 44', fontsize=fontsize)
cg.plot(CG_44 - CG_model_44, sig=2, min=-5, max=5, cbar=False,   sub=(5, 4, 16),
        fontsize=fontsize)
#cg.plot(PR4_30 - CG_model_30, sig=1, min=-5, max=5, cbar=False,  sub=(5, 4, 17),
#    llabel=r'\mathrm{PR4}\ 30', fontsize=fontsize)
#cg.plot(PR4_30 - CG_model_30, sig=2, min=-5, max=5, cbar=False,  sub=(5, 4, 18),
#        fontsize=fontsize)
#cg.plot(PR4_44 - CG_model_44, sig=1, min=-5, max=5, cbar=False,  sub=(5, 4, 19),
#    llabel=r'\mathrm{PR4}\ 44', fontsize=fontsize)
#cg.plot(PR4_44 - CG_model_44, sig=2, min=-5, max=5, cbar=False,  sub=(5, 4, 20),
#        fontsize=fontsize)

plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../figures/CG_DR1_residuals.pdf', bbox_inches='tight', dpi=150)
plt.close('all')


cg.plot(CG_30 - CG_model_30, sig=1, min=-5, max=5, cbar=False,
    llabel=r'\mathrm{CG}\ 30')
cg.plot(CG_30 - CG_model_30, sig=2, min=-5, max=5, cbar=False)
cg.plot(CG_44 - CG_model_44, sig=1, min=-5, max=5, cbar=False,
    llabel=r'\mathrm{CG}\ 44')
cg.plot(CG_44 - CG_model_44, sig=2, min=-5, max=5, cbar=False)


'''
model = cg.sky_model_from_chain(f'{DIR_BP}/BP_c0001_v2.h5', nside=512,
        components=['synch', 'dust'])

data = h5py.File(f'{data_dir}/LFI_instrument_v8.h5')
nu = data['030/bandpassx'] * u.GHz
weights = data['030/bandpass'] * u.Unit("K_RJ")
BP_model_30 = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value

nu = data['044/bandpassx'] * u.GHz
weights = data['044/bandpass'] * u.Unit("K_RJ")
BP_model_44 = model(nu, weights, fwhm=fwhm*u.arcmin, output_unit='uK_RJ').value

cg.plot(BP_30 - BP_model_30, sig=1, min=-5, max=5, cbar=False,
    llabel=r'\mathrm{BP}\ 30')
cg.plot(BP_30 - BP_model_30, sig=2, min=-5, max=5, cbar=False)
cg.plot(BP_44 - BP_model_44, sig=1, min=-5, max=5, cbar=False,
    llabel=r'\mathrm{BP}\ 44')
cg.plot(BP_44 - BP_model_44, sig=2, min=-5, max=5, cbar=False)

plt.show()
'''
