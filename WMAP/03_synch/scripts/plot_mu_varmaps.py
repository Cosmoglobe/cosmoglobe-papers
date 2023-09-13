import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

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

var_beta = hp.read_map(f'{DIR}/var_beta.fits')
var_beta_persamp = hp.read_map(f'{DIR}/var_beta_persamp.fits')
mu = hp.read_map(f'{DIR}/mean_beta.fits')


sd_beta= var_beta**0.5
sd_beta_persamp= var_beta_persamp**0.5

#sd_beta[inds] = np.nan
#sd_beta_persamp[inds] = np.nan
#mu[inds] = np.nan

cg.plot(mu, min=-3.5, max=-2.5, llabel=r'\beta_{\mathrm{s}}', width=4,
        extend='both')
plt.savefig('../figures/beta_n0016_mu.pdf', bbox_inches='tight')
cg.plot(sd_beta, cmap='bone', min=0., max=0.15,
    llabel=r'\sigma_{\beta_{\mathrm{s}}}^{\mathrm{sys+stat}}', width=4, extend='both')
plt.savefig('../figures/beta_n0016_sd_tot.pdf', bbox_inches='tight')
cg.plot(sd_beta_persamp, cmap='bone', min=0., max=0.1,
    llabel=r'\sigma_{\beta_{\mathrm{s}}}^{\mathrm{stat}}', width=4, extend='both')
plt.savefig('../figures/beta_n0016_sd_samp.pdf', bbox_inches='tight')
