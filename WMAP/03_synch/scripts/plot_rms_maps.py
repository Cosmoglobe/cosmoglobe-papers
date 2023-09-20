import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cosmoglobe as cg
from tqdm import tqdm

w = 3.5

fontsize = {'llabel':9, 'rlabel':9}
dpi = 150

def update(existingAggregate, newValue):
    (count, mean, M2) = existingAggregate
    count += 1
    delta = newValue - mean
    mean += delta / count
    delta2 = newValue - mean
    M2 += delta * delta2
    return (count, mean, M2)

# Retrieve the mean, variance and sample variance from an aggregate
def finalize(existingAggregate):
    (count, mean, M2) = existingAggregate
    if count < 2:
        return float("nan")
    else:
        (mean, variance, sampleVariance) = (mean, M2 / count, M2 / (count - 1))
        return (mean, variance, sampleVariance)

CG_DIR = '/mn/stornext/d16/cmbco/cg/v1'


d_K = hp.read_map(f'{CG_DIR}/CG_023-WMAP_K_IQU_n0512_v1.fits',
    field=(0,1,2,4,5,6))
d_30 = hp.read_map(f'{CG_DIR}/CG_030_IQU_n0512_v1.fits', field=(0,1,2,4,5,6))

# Smooth *to* 1 degree
fwhm = np.sqrt((72/60)**2 - 0.88**2)*np.pi/180
d_K2 = hp.smoothing(d_K[:3], fwhm=fwhm)
fwhm = np.sqrt((72/60)**2 - (30/60)**2)*np.pi/180
d_302 = hp.smoothing(d_30[:3], fwhm=fwhm)

P_K = hp.ud_grade(np.hypot(d_K2[1], d_K2[2])*1e3*0.47, 128)
P_30 = hp.ud_grade(np.hypot(d_302[1], d_302[2]), 128)

cg.plot(P_K, min=0, max=50, comp='synch', 
        norm='linear',
        fontsize=fontsize,
    llabel=r'K',
    unit=r'\mathrm{\mu K\,@\,30\,GHz}', rlabel='P', cbar=False, width=w)
plt.savefig('../figures/Kband_polint.pdf', bbox_inches='tight', dpi=dpi)
cg.plot(P_30, min=0, max=50, comp='synch',
        norm='linear',
        fontsize=fontsize,
    llabel=r'30',
    unit=r'\mathrm{\mu K\,@\,30\,GHz}', rlabel='P', cbar=False, width=w)
plt.savefig('../figures/30GHz_polint.pdf', bbox_inches='tight', dpi=dpi)


sigmaQU_K = hp.read_map('CG_023_72_n0128.fits', field=(1,2))
sigmaP_K = np.hypot(sigmaQU_K[0], sigmaQU_K[1])*1e3*(23/30)**3.1

sigmaQU_30 = hp.read_map('CG_030_72_n0128.fits', field=(1,2))
sigmaP_30 = np.hypot(sigmaQU_30[0], sigmaQU_30[1])*(28/30)**3.1

cg.plot(sigmaP_K, min=0, max=7.5, cmap='binary_r', unit=r'\mathrm{\mu K}',
        fontsize=fontsize,
    rlabel=r'\sigma_P', llabel=r'K\ \mathrm{(synch. scaled)}', cbar=False,
    width=w)
plt.savefig('../figures/Kband_sigmaP.pdf', bbox_inches='tight', dpi=dpi)
cg.plot(sigmaP_30, cmap='binary_r', unit=r'\mathrm{\mu K}',
        fontsize=fontsize,
    rlabel=r'\sigma_P', llabel=r'30', min=0, max=7.5, cbar=False,
    width=w)
plt.savefig('../figures/30GHz_sigmaP.pdf', bbox_inches='tight', dpi=dpi)
plt.close('all')

cg.plot(P_K/sigmaP_K, min=0, max=10, cmap='bone')
cg.plot(P_30/sigmaP_30, min=0, max=10, cmap='bone')

DIR_CG = '/mn/stornext/d16/cmbco/cg/v1'
synch_CG, h = hp.read_map(f'{DIR_CG}/CG_synch_IQU_n1024_v1.fits',
        field=(0,1,2,7,8), h=True)

P_synch = np.hypot(synch_CG[1], synch_CG[2])
sigmaP_synch = np.hypot(synch_CG[3], synch_CG[4])
cg.plot(P_synch/sigmaP_synch, min=0, max=10, cmap='bone', width=w)

sigmaP_synch = np.hypot(synch_CG[3], synch_CG[4])
cg.plot(sigmaP_synch, min=0, max=7.5, cmap='binary_r',
        fontsize=fontsize,
    llabel=r'\mathrm{Cosmoglobe}',
    rlabel='\sigma_P', 
    unit=r'\mathrm{\mu K}',
    width=w)
plt.savefig('../figures/polint_CG_sigma.pdf', bbox_inches='tight', dpi=dpi)

print(sigmaP_K.mean())
print(sigmaP_30.mean())
print(sigmaP_synch.mean())

#plt.show()
