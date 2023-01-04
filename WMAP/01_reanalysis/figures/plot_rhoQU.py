import healpy as hp
import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KaQVW_221201'
d_Ka = hp.read_map(f'{DIR}/tod_030-WMAP_Ka_rms_c0001_k000001.fits', field=(1,2,3))
d_30 = hp.read_map(f'{DIR}/tod_030_rms_c0001_k000001.fits', field=(1,2,3))

rho_Ka = d_Ka[2]/np.sqrt(d_Ka[0]*d_Ka[1])
rho_30 = d_30[2]/np.sqrt(d_30[0]*d_30[1])

cmap = 'RdBu_r'
cg.plot(rho_Ka, min=-0.5, max=0.5, llabel='\mathit{Ka}', rlabel=r'\rho_{QU}',
    sub=(2,1,1), cmap=cmap, width=4, xsize=1000)
cg.plot(rho_30, min=-0.1, max=0.1, llabel='30', rlabel=r'\rho_{QU}',
    sub=(2,1,2), cmap=cmap, width=4, xsize=1000)
plt.savefig('rho_QU.pdf')
