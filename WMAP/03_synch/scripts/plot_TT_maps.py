import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

pts_to_inches = 1/72
columnwidth = 256.0748*pts_to_inches

DIR = '/mn/stornext/d5/data/duncanwa/comm1'

m = hp.read_map(f'{DIR}/UF_sindex_regions_full.fits')
inds = (m == hp.UNSEEN)

var_beta = hp.read_map(f'{DIR}/var_beta_batmask.fits')
var_beta_persamp = hp.read_map(f'{DIR}/var_beta_persamp_batmask.fits')
mu = hp.read_map(f'{DIR}/mean_beta_batmask.fits')


sd_beta= var_beta**0.5
sd_beta_persamp= var_beta_persamp**0.5

sd_beta[inds] = np.nan
sd_beta_persamp[inds] = np.nan
mu[inds] = np.nan


'''
cg.plot(mu, min=-3.5, max=-2.5, cmap='jet')
cg.plot(sd_beta, min=0, max=0.15, cmap='viridis')
cg.plot(sd_beta_persamp, min=0, max=0.15, cmap='viridis')
'''

fnames = glob('TT*.fits')
fnames.sort()

labels = [r'\mathrm{Cosmoglobe}', 'CG K/CG Ka', r'\mathit{WMAP9}/\mathrm{PR3}', 
        r'\mathit{WMAP9}/\mathrm{PR4}', r'\mathit{WMAP9}/\mathrm{BP}']
fignames = ['CG_K30', 'CG_KKa', 'W9PR3_K30', 'W9PR4_K30', 'W9BP_K30']
for i, f in enumerate(fnames):
    if labels[i] == labels[2]:
        llabel = r'\beta_{\mathrm{s}}'
    else:
        llabel = ''
    cg.plot(f,
            min=-5, max=-1, llabel=llabel, width=columnwidth,
            cbar=False,  rlabel=labels[i], 
            fontsize={'rlabel':8, 'llabel':10}
            )
    plt.savefig(f'TT_map_{fignames[i]}.pdf', bbox_inches='tight', dpi=150)
    plt.close()

cg.standalone_colorbar("planck", 
        ticks=np.arange(-5, -0.5, 0.5),
        ticklabels=['$-5$','','','', '$-3$','','','', '$-1$'], 
        extend='both',
        width=0.8*columnwidth)
plt.savefig('cbar_beta_wide.pdf', bbox_inches='tight', dpi=150)
