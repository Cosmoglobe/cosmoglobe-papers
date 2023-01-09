import healpy as hp
import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_b_230105'

width = 4
xsize = 1200

fnames = glob(f'{DIR}/tod*WMAP*map*k000019.fits')
fnames.sort()
print(fnames)

K = hp.read_map(fnames[0], field=(0,1,2))
Ka = hp.read_map(fnames[1], field=(0,1,2))
Q1 = hp.read_map(fnames[2], field=(0,1,2))
Q2 = hp.read_map(fnames[3], field=(0,1,2))
V1 = hp.read_map(fnames[4], field=(0,1,2))
V2 = hp.read_map(fnames[5], field=(0,1,2))
W1 = hp.read_map(fnames[6], field=(0,1,2))
W2 = hp.read_map(fnames[7], field=(0,1,2))
W3 = hp.read_map(fnames[8], field=(0,1,2))
W4 = hp.read_map(fnames[9], field=(0,1,2))


d1 = hp.smoothing(Ka - 0.32*K, fwhm = 2*np.pi/180)*1e3
d2 = hp.smoothing(Q1 - Q2, fwhm = 2*np.pi/180)*1e3
d3 = hp.smoothing(V1 - V2, fwhm = 2*np.pi/180)*1e3
d4 = hp.smoothing(((W1-W2)-(W3-W4))/4, fwhm = 2*np.pi/180)*1e3


cg.plot(d1, sig=1, llabel=r'\mathit K-\mathit{Ka}', rlabel='Q', cbar=False,
    min=-10, max=10, xsize=xsize, width=width)
plt.savefig('KKa_deltaQ.pdf', bbox_inches='tight')
cg.plot(d1, sig=2, rlabel='U', cbar=False, min=-10, max=10, xsize=xsize,
    width=width)
plt.savefig('KKa_deltaU.pdf', bbox_inches='tight')
cg.plot(d2, sig=1, llabel=r'\Delta Q', min=-10, max=10, cbar=False, rlabel='Q',
    xsize=xsize, width=width)
plt.savefig('Q_deltaQ.pdf', bbox_inches='tight')
cg.plot(d2, sig=2, rlabel='U', cbar=False, min=-10, max=10, xsize=xsize,
    width=width)
plt.savefig('Q_deltaU.pdf', bbox_inches='tight')
cg.plot(d3, sig=1, llabel=r'\Delta V', min=-10, max=10, cbar=False, rlabel='Q',
    xsize=xsize, width=width)
plt.savefig('V_deltaQ.pdf', bbox_inches='tight')
cg.plot(d3, sig=2, rlabel='U', cbar=False, min=-10, max=10, xsize=xsize,
    width=width)
plt.savefig('V_deltaU.pdf', bbox_inches='tight')
cg.plot(d4, sig=1, unit=r'\mathrm{\mu K}', llabel=r'\Delta W', rlabel = 'Q',
    min=-10, max=10, xsize=xsize, width=width)
plt.savefig('W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d4, sig=2, rlabel='U', unit=r'\mathrm{\mu K}', min=-10, max=10,
    xsize=xsize, width=width)
plt.savefig('W_deltaU.pdf', bbox_inches='tight')
plt.show()
