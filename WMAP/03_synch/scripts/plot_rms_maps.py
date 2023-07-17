import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cosmoglobe as cg

CG_DIR = '/mn/stornext/d16/cmbco/cg/v1'


d_K = hp.read_map(f'{CG_DIR}/CG_023-WMAP_K_IQU_n0512_v1.fits',
    field=(0,1,2,4,5))*1e3*0.47
d_30 = hp.read_map(f'{CG_DIR}/CG_030_IQU_n0512_v1.fits', field=(0,1,2,4,5))

# Smooth *to* 1 degree
fwhm = np.sqrt(1**2 - 0.88**2)*np.pi/180
d_K2 = hp.smoothing(d_K[:3], fwhm=fwhm)
fwhm = np.sqrt(1**2 - (30/60)**2)*np.pi/180
d_302 = hp.smoothing(d_30[:3], fwhm=fwhm)

P_K = np.hypot(d_K2[1], d_K2[2])
P_30 = np.hypot(d_302[1], d_302[2])

cg.plot(P_K, min=0, max=50, comp='synch', 
    llabel=r'K\ \mathrm{(synch.\ scaled})',
    unit=r'\mathrm{\mu K\,@\,30\,GHz}', rlabel='P')
plt.savefig('../figures/Kband_polint.pdf', bbox_inches='tight')
cg.plot(P_30, min=0, max=50, comp='synch',
    llabel=r'30',
    unit=r'\mathrm{\mu K\,@\,30\,GHz}', rlabel='P')
plt.savefig('../figures/30GHz_polint.pdf', bbox_inches='tight')

sigmaP_K = np.hypot(d_K[3], d_K[4])
sigmaP_30 = np.hypot(d_30[3], d_30[4])

cg.plot(sigmaP_K, min=0, max=70, cmap='binary_r', unit=r'\mathrm{\mu K\,@\,30\,GHz}',
    rlabel=r'\sigma_P', llabel=r'K\ \mathrm{(synch. scaled)}')
plt.savefig('../figures/Kband_sigmaP.pdf', bbox_inches='tight')
cg.plot(sigmaP_30, min=0, max=70, cmap='binary_r', unit=r'\mathrm{\mu K}',
    rlabel=r'\sigma_P', llabel=r'30\,\mathrm{GHz}')
plt.savefig('../figures/30GHz_sigmaP.pdf', bbox_inches='tight')
plt.show()
