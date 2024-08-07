import healpy as hp
import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

from setup_matplotlib import *


dpi = 100

DIR = "/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_a_230206"

width = 2/3
xsize = 1200
fontsize = {
    "llabel": 8,
    "rlabel": 8,
}
width = cm2inch(9)
width = 6

fnames = glob(f"{DIR}/tod*WMAP*res*k??????.fits")
fnames.sort()
fwhm = 5 * np.pi / 180

bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

fnames_i = []
for b in bands:
  for f in fnames:
    if b + '_' in f:
      fi = f
  fnames_i.append(fi)
fnames = fnames_i
K  = hp.read_map(fnames[0], field=(0, 1, 2))*1e3
Ka = hp.read_map(fnames[1], field=(0, 1, 2))*1e3
Q1 = hp.read_map(fnames[2], field=(0, 1, 2))*1e3
Q2 = hp.read_map(fnames[3], field=(0, 1, 2))*1e3
V1 = hp.read_map(fnames[4], field=(0, 1, 2))*1e3
V2 = hp.read_map(fnames[5], field=(0, 1, 2))*1e3
W1 = hp.read_map(fnames[6], field=(0, 1, 2))*1e3
W2 = hp.read_map(fnames[7], field=(0, 1, 2))*1e3
W3 = hp.read_map(fnames[8], field=(0, 1, 2))*1e3
W4 = hp.read_map(fnames[9], field=(0, 1, 2))*1e3

maps = [K, Ka, Q1, Q2, V1, V2, W1, W2, W3, W4]


for i,m in enumerate(maps):
    cg.plot(m, sig=0, min=-10, max=10, fwhm=5*u.deg, cbar=False,
        llabel=bands[i], rlabel='T', sub=(1,3,1), width=width,
        fontsize=fontsize)
    cg.plot(m, sig=1, min=-10, max=10, fwhm=5*u.deg, cbar=False, rlabel='Q',
        sub=(1,3,2), width=width, fontsize=fontsize)
    cg.plot(m, sig=2, min=-10, max=10, fwhm=5*u.deg, cbar=False, rlabel='U',
        sub=(1,3,3), width=width, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f'tod_res_{bands[i]}_IQU.pdf', bbox_inches='tight', dpi=dpi)
    plt.close('all')


fnames = glob(f"{DIR}/res_*WMAP*k??????.fits")
fnames.sort()

fnames_i = []
for b in bands:
  for f in fnames:
    if b + '_' in f:
      fi = f
  fnames_i.append(fi)
fnames = fnames_i
K  = hp.read_map(fnames[0], field=(0, 1, 2))*1e3
Ka = hp.read_map(fnames[1], field=(0, 1, 2))*1e3
Q1 = hp.read_map(fnames[2], field=(0, 1, 2))*1e3
Q2 = hp.read_map(fnames[3], field=(0, 1, 2))*1e3
V1 = hp.read_map(fnames[4], field=(0, 1, 2))*1e3
V2 = hp.read_map(fnames[5], field=(0, 1, 2))*1e3
W1 = hp.read_map(fnames[6], field=(0, 1, 2))*1e3
W2 = hp.read_map(fnames[7], field=(0, 1, 2))*1e3
W3 = hp.read_map(fnames[8], field=(0, 1, 2))*1e3
W4 = hp.read_map(fnames[9], field=(0, 1, 2))*1e3

maps = [K, Ka, Q1, Q2, V1, V2, W1, W2, W3, W4]

for i,m in enumerate(maps):
    cg.plot(m, sig=0, min=-10, max=10, fwhm=5*u.deg, cbar=False,
        llabel=bands[i], rlabel='T', sub=(1,3,1), width=width,
        fontsize=fontsize)
    cg.plot(m, sig=1, min=-10, max=10, fwhm=5*u.deg, cbar=False, rlabel='Q',
        sub=(1,3,2), width=width, fontsize=fontsize)
    cg.plot(m, sig=2, min=-10, max=10, fwhm=5*u.deg, cbar=False, rlabel='U',
        sub=(1,3,3), width=width, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f'compsep_res_{bands[i]}_IQU.pdf', bbox_inches='tight', dpi=dpi)
    cg.plot(m, sig=0, min=-5, max=5, fwhm=5*u.deg, cbar=False,
        llabel=bands[i], rlabel='T', width=2)
    plt.savefig(f'compsep_res_{bands[i]}_I.pdf', bbox_inches='tight', dpi=dpi)
    cg.plot(m, sig=1, min=-5, max=5, fwhm=5*u.deg, cbar=False, rlabel='Q',
        width=2)
    plt.savefig(f'compsep_res_{bands[i]}_Q.pdf', bbox_inches='tight', dpi=dpi)
    cg.plot(m, sig=2, min=-5, max=5, fwhm=5*u.deg, cbar=False, rlabel='U',
        width=2)
    plt.savefig(f'compsep_res_{bands[i]}_U.pdf', bbox_inches='tight', dpi=dpi)
    plt.close('all')


cg.standalone_colorbar("planck", ticks=[-5,0,5], extend='both',
        unit=r"$\mathrm{\mu K_{CMB}}$", fontsize=18, width=6)

plt.savefig('cbar_5uK.pdf', bbox_inches='tight')

cg.standalone_colorbar("planck", ticks=[-3,0,3], extend='both',
        unit=r"$\mathrm{\mu K_{CMB}}$", fontsize=18, width=6)
plt.savefig('cbar_3uK.pdf', bbox_inches='tight')

cg.standalone_colorbar("planck", ticks=[-10,0,10], extend='both',
        unit=r"$\mathrm{\mu K_{CMB}}$", fontsize=18, width=6)
plt.savefig('cbar_10uK.pdf', bbox_inches='tight')

cg.standalone_colorbar("planck", ticks=[-10,0,10], extend='both',
        unit=r"$\mathrm{\mu K_{CMB}}$", fontsize=18, width=3)
plt.savefig('cbar_10uK_4in.pdf', bbox_inches='tight')
