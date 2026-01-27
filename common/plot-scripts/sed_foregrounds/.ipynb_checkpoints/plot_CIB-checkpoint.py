import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
import healpy as hp
import cosmoglobe as cg

# Theoretical predictions

data = np.loadtxt('EBL_intensity_total_z0.00.dat')
wav = data[:,0]
nu_bla2 = (c.c/(wav*u.angstrom)).to('GHz')
nuInu = data[:,1]*u.nW/u.m**2/u.sr
I_nu = (nuInu/nu_bla2).to('MJy sr-1')

wav = (wav*u.angstrom).to('micron')

#C = (2*c.k_B*nu_bla2**2/c.c**2).to('MJy/uK')
#C_bla = (c.c**2/(2*c.k_B*nu_bla2**2)).to('uK/MJy')

plt.loglog(wav, I_nu, label=r'CIB & COB (Finke et al. 2022)', color='k',
        linestyle='--')
plt.fill_between(wav.value, I_nu.value/2, 2*I_nu.value, color='k', alpha=0.25)
'''
plt.plot(lamb_finke, nuInu_finke, 'k--', label='Finke et al.')


# FIRAS measurements, Planck model
plt.plot(lamb_planck, nuInu_planck, label=r'Planck')


'''

lamb_dirbe = np.array([1.24, 2.2, 3.5, 4.9, 12, 25, 60, 100, 140, 240])*u.micron

data = np.loadtxt('/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.5/chains_v0.5/md_c0001_k000019.dat',
        usecols=(1))

mono_DIRBE = data[::2][:10]
plt.plot(lamb_dirbe, mono_DIRBE, 'ro', label='Cosmoglobe - DIRBE v0.5')


DIR = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.2/chains_v0.2'
OLA = '/mn/stornext/d16/cmbco/ola/dirbe'


pix_256 = np.arange(12*256**2)
lon, lat = hp.pix2ang(256, pix_256, lonlat=True)

e2g = hp.Rotator(coord=['E', 'G'])
lat_ecl = e2g.rotate_map_alms(lat)
mask_256 = (abs(lat) < 45) | (abs(lat_ecl) < 30)

comps = ['co_tot', 'dust', 'dust_cii', 'ff', 'hotPAH', 'stars']

mono_DIRBE_ZSMA = []
for b in range(1, 11):
    m = np.zeros(12*512**2)
    samp = 10
    band = f'{b:02}b'
    for comp in comps:
        m += hp.read_map(f'{DIR}/{comp}_{band}_c0001_k{samp:06}.fits')
    
    
    m = hp.ma(m)
    m_ZSMA = hp.read_map(f'{OLA}/DIRBE_ZSMA_{b:02}_1_256.fits')
    
    m = hp.ud_grade(m, 256)
    diff = hp.ma(m_ZSMA - m)
    diff[mask_256] = hp.UNSEEN
    #cg.plot(diff)
    mono_DIRBE_ZSMA.append(hp.remove_monopole(diff, fitval=True)[1])
#plt.show()

print(mono_DIRBE)
print(mono_DIRBE_ZSMA)
plt.plot(lamb_dirbe, mono_DIRBE_ZSMA, 'ro', markerfacecolor='white', label='DIRBE ZSMA')

mono_FIRAS = data[26:]
nu_FIRAS = np.array([108, 149, 217, 353, 544, 857, 1251, 1904, 2135,
    2802])*u.GHz
lamb_FIRAS = (c.c/nu_FIRAS).to('micron')


plt.plot(lamb_FIRAS, mono_FIRAS, 'bo', label='Cosmoglobe - FIRAS v0.5')



nu = np.array([100, 143, 217, 353, 545, 857])*u.GHz
wav = (c.c/nu).to('micron')

I_nu = np.array([0.007, 0.01, 0.06, 0.149, 0.371, 0.576])*u.MJy/u.sr
sigma = np.array([0.014, 0.019, 0.023, 0.017, 0.018, 0.034])*u.MJy/u.sr

plt.plot(wav[-5:], I_nu[-5:], 'o-', label='Planck (Odegard et al. 2019)')

plt.xlim([1, 2500])
plt.ylim([1e-3, 2.5])
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='lower right', framealpha=0)

plt.xlabel(r'Wavelength ($\mathrm{\mu m}$)')
plt.ylabel(r'CIB Monopole $I_\nu$ (MJy/sr)')
plt.savefig('CIB_mono.png', bbox_inches='tight')

plt.show()
