import healpy as hp
import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_b_230105'

width = 2
xsize = 1200
fontsize = {'llabel': 5,
            'rlabel': 5,}

fnames = glob(f'{DIR}/tod*WMAP*map*k000019.fits')
fnames.sort()
fname_30 = glob(f'{DIR}/tod_030_*map*k000019.fits')[0]
d_30 = hp.read_map(fname_30, field=(0,1,2))*1e-3
fwhm = 5*np.pi/180

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


d0 = hp.ud_grade(hp.smoothing(0.47*K - d_30, fwhm = fwhm)*1e3, 256)
d1 = hp.ud_grade(hp.smoothing(0.63*d_30 - Ka, fwhm=fwhm)*1e3, 256)
d2 = hp.ud_grade(hp.smoothing(Q1 - Q2, fwhm = fwhm)*1e3, 256)
d3 = hp.ud_grade(hp.smoothing(V1 - V2, fwhm = fwhm)*1e3, 256)
d4 = hp.ud_grade(hp.smoothing(((W1-W2)-(W3-W4))/4, fwhm = fwhm)*1e3, 256)


cg.plot(d0, sig=1, llabel=r'\mathit{K}-30', rlabel='Q, \mathrm{CG}', cbar=False,
    min=-10, max=10, xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('K30_deltaQ.pdf', bbox_inches='tight')
cg.plot(d0, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('K30_deltaU.pdf', bbox_inches='tight')

cg.plot(d1, sig=1, llabel=r'30-\mathit{Ka}', rlabel='Q, \mathrm{CG}', cbar=False,
    min=-10, max=10, xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('30K_deltaQ.pdf', bbox_inches='tight')
cg.plot(d1, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('30K_deltaU.pdf', bbox_inches='tight')

cg.plot(d2, sig=1, llabel=r'\Delta Q', min=-10, max=10, cbar=False, rlabel='Q, \mathrm{CG}',
    xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('Q_deltaQ.pdf', bbox_inches='tight')
cg.plot(d2, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('Q_deltaU.pdf', bbox_inches='tight')
cg.plot(d3, sig=1, llabel=r'\Delta V', min=-10, max=10, cbar=False, rlabel='Q, \mathrm{CG}',
    xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('V_deltaQ.pdf', bbox_inches='tight')
cg.plot(d3, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('V_deltaU.pdf', bbox_inches='tight')
cg.plot(d4, sig=1, unit=r'\mathrm{\mu K}', llabel=r'\Delta W', rlabel = 'Q, \mathrm{CG}',
    min=-10, max=10, xsize=xsize, width=width, extend='both', fontsize=fontsize)
plt.savefig('W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d4, sig=2, rlabel='U, \mathrm{CG}', unit=r'\mathrm{\mu K}', min=-10, max=10,
    xsize=xsize, width=width, extend='both', fontsize=fontsize)
plt.savefig('W_deltaU.pdf', bbox_inches='tight')
plt.close('all')

DIR = '/mn/stornext/d16/cmbco/ola/wmap/freq_maps'
fnames = glob(f'{DIR}/wmap_iqusmap_r9_9yr*.fits')
fnames.sort()

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


d0 = hp.ud_grade(hp.smoothing(0.47*K - d_30, fwhm = fwhm)*1e3, 256)
d1 = hp.ud_grade(hp.smoothing(0.63*d_30 - Ka, fwhm=fwhm)*1e3, 256)
d2 = hp.ud_grade(hp.smoothing(Q1 - Q2, fwhm = fwhm)*1e3, 256)
d3 = hp.ud_grade(hp.smoothing(V1 - V2, fwhm = fwhm)*1e3, 256)
d4 = hp.ud_grade(hp.smoothing(((W1-W2)-(W3-W4))/4, fwhm = fwhm)*1e3, 256)


cg.plot(d0, sig=1, llabel=r'\mathit{K}-30', rlabel='Q, \mathrm{CG}', cbar=False,
    min=-10, max=10, xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('K30_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d0, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('K30_W_deltaU.pdf', bbox_inches='tight')

cg.plot(d1, sig=1, llabel=r'30-\mathit{Ka}', rlabel='Q, \mathrm{CG}', cbar=False,
    min=-10, max=10, xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('30K_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d1, sig=2, rlabel='U, \mathrm{CG}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('30K_W_deltaU.pdf', bbox_inches='tight')

cg.plot(d1, sig=1, llabel=r'\mathit K-\mathit{Ka}', rlabel='Q, \mathit{WMAP}', cbar=False,
    min=-10, max=10, xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('KKa_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d1, sig=2, rlabel='U, \mathit{WMAP}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('KKa_W_deltaU.pdf', bbox_inches='tight')
cg.plot(d2, sig=1, llabel=r'\Delta Q', min=-10, max=10, cbar=False, rlabel='Q, \mathit{WMAP}',
    xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('Q_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d2, sig=2, rlabel='U, \mathit{WMAP}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('Q_W_deltaU.pdf', bbox_inches='tight')
cg.plot(d3, sig=1, llabel=r'\Delta V', min=-10, max=10, cbar=False, rlabel='Q, \mathit{WMAP}',
    xsize=xsize, width=width, fontsize=fontsize)
plt.savefig('V_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d3, sig=2, rlabel='U, \mathit{WMAP}', cbar=False, min=-10, max=10, xsize=xsize,
    width=width, fontsize=fontsize)
plt.savefig('V_W_deltaU.pdf', bbox_inches='tight')
cg.plot(d4, sig=1, unit=r'\mathrm{\mu K}', llabel=r'\Delta W', rlabel = 'Q, \mathit{WMAP}',
    min=-10, max=10, xsize=xsize, width=width, extend='both', fontsize=fontsize)
plt.savefig('W_W_deltaQ.pdf', bbox_inches='tight')
cg.plot(d4, sig=2, rlabel='U, \mathit{WMAP}', unit=r'\mathrm{\mu K}', min=-10, max=10,
    xsize=xsize, width=width, extend='both', fontsize=fontsize)
plt.savefig('W_W_deltaU.pdf', bbox_inches='tight')
