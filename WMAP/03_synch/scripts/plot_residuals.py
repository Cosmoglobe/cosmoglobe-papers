import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
import cosmoglobe as cg

fontsize = {'cbar_tick_label':8, 'cbar_label':10}

orp = {'cbar_label_pad':-10}

DIR = '/mn/stornext/d5/data/duncanwa/comm1'

bands = ['023', '030', '033', '041_Q1', '041_Q2', '044', '061_V1', '061_V2', '070',
    '094_W1', '094_W2', '094_W3', '094_W4', '353']

band_labels = ['K', '30', 'Ka', 'Q1', 'Q2', '44', 'V1', 'V2', '70','W1', 'W2',
        'W3', 'W4', '353']

wmap = [True, False, True, True, True, False, True, True, False,
        True, True, True, True, False]


fnames_res = glob(f'{DIR}/cg_plots/res*.fits')
fnames_res.sort()
fnames_res = fnames_res[:-1]

fnames = glob(f'{DIR}/maps/CG*_sigma0*600arc_uK*.fits')
fnames.sort()

fnames[1], fnames[2] = fnames[2], fnames[1]

cg.plot(np.arange(12.), cbar=False, sub=(2,4,1))
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
n = 0
rlabQ = 'Q'
rlabU = 'U'
for i in range(len(fnames_res)):
    cbar = False
    if wmap[i]:
        continue
    if band_labels[i] == '70':
        cbar = True
    lim = 2.5
    cg.plot(fnames_res[i], sig=1, sub=(2,4,2*n+1), llabel=band_labels[i],
            min=-lim, max=lim, cbar=cbar, rlabel=rlabQ, fontsize=fontsize,
            extend='both', unit=r'$\mathrm{\mu K}$',
            override_plot_properties=orp)
    cg.plot(fnames_res[i], sig=2, sub=(2,4,2*n+2),
            min=-lim, max=lim, cbar=cbar, rlabel=rlabU, fontsize=fontsize,
            extend='both', unit=r'$\mathrm{\mu K}$',
            override_plot_properties=orp)
    n += 1
    if n == 2:
        rlabQ = ''
        rlabU = ''
plt.savefig('../figures/comm1_res_QU_LFI.pdf', bbox_inches='tight', dpi=150)
plt.close('all')

n = 0
cg.plot(np.arange(12.), cbar=False, sub=(3,4,1))
plt.subplots_adjust(wspace=0.1, hspace=-0.1)

for i in range(len(fnames_res)):
    if band_labels[i] == 'W3' or band_labels[i] == 'W4':
        cbar = True
    elif band_labels[i] == 'V1' or band_labels[i] == 'V2':
        cbar = True
    else:
        cbar = False
    if wmap[i]:
        pass
    else:
        continue
    if 'W' in band_labels[i]:
        lim = 10
    else:
        lim = 5
    if 'W' in band_labels[i]:
        cg.plot(fnames_res[i], sig=1, sub=(2,4,2*n+1), llabel=band_labels[i],
                min=-lim, max=lim, cbar=cbar, fontsize=fontsize, extend='both',
                unit=r'$\mathrm{\mu K}$', override_plot_properties=orp)
        cg.plot(fnames_res[i], sig=2, sub=(2,4,2*n+2),
                min=-lim, max=lim, cbar=cbar, fontsize=fontsize, extend='both',
                unit=r'$\mathrm{\mu K}$', override_plot_properties=orp)
    else:
        cg.plot(fnames_res[i], sig=1, sub=(3,4,2*n+1), llabel=band_labels[i],
                min=-lim, max=lim, cbar=cbar, fontsize=fontsize, extend='both',
                unit=r'$\mathrm{\mu K}$', override_plot_properties=orp)
        cg.plot(fnames_res[i], sig=2, sub=(3,4,2*n+2),
                min=-lim, max=lim, cbar=cbar, fontsize=fontsize, extend='both',
                unit=r'$\mathrm{\mu K}$', override_plot_properties=orp)
    n += 1
    if band_labels[i] == 'V2':
        plt.savefig('../figures/comm1_res_QU_K-V.pdf', bbox_inches='tight',
                dpi=300)
        plt.close('all') 
        n = 0
        cg.plot(np.arange(12.), cbar=False, sub=(2,4,1))
        plt.subplots_adjust(wspace=0.1, hspace=-0.1)
    if band_labels[i] == 'W4':
        plt.savefig('../figures/comm1_res_QU_W.pdf', bbox_inches='tight',
                dpi=300)
        plt.close('all')
        n = 0

cg.plot(np.arange(12.), cbar=False, sub=(1,4,1))
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
m = hp.read_map(f'{DIR}/chisq_PQU.fits', field=(1,2))
mu = 14
sd = (2*mu)**0.5
cg.plot((m-mu)/sd, sig=0, cmap='RdBu_r',
        llabel=r'\chi^2', sub=(1,4,1), cbar=True, extend='both',
        fontsize=fontsize, override_plot_properties=orp, ticks=[-2,0,2])
cg.plot((m-mu)/sd, sig=1, cmap='RdBu_r',
        sub=(1,4,2), cbar=True, extend='both', fontsize=fontsize,
        override_plot_properties=orp, ticks=[-2,0,2])
#plt.tight_layout()
plt.savefig('../figures/comm1_res_QU_chisq.pdf', bbox_inches='tight', dpi=150)
