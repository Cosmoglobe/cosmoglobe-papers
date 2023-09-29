import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.colors import ListedColormap


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


mu_30 = hp.read_map(f'{DIR}/mean_beta_-3.0.fits')

prior_var = abs(mu-mu_30)


sd_beta= var_beta**0.5
sd_beta_persamp= var_beta_persamp**0.5

#sd_beta[inds] = np.nan
#sd_beta_persamp[inds] = np.nan
#mu[inds] = np.nan

cg.plot(mu, min=-3.5, max=-2.5, llabel=r'\beta_{\mathrm{s}}', width=4,
        extend='both')
plt.savefig('../figures/beta_n0016_mu.pdf', bbox_inches='tight')
cg.plot(sd_beta, cmap='bone', min=0., max=0.15,
    llabel=r'\sigma_{\beta_{\mathrm{s}}}^{\mathrm{sys+stat}}', width=4,
    extend='both', cbar=False)
plt.savefig('../figures/beta_n0016_sd_stat_inst.pdf', bbox_inches='tight')


cg.plot(np.hypot(prior_var, sd_beta), cmap='bone', min=0., max=0.15,
    llabel=r'\sigma_{\beta_{\mathrm{s}}}^{\mathrm{sys+stat+prior}}', width=4, extend='both')
plt.savefig('../figures/beta_n0016_sd_stat_inst_prior.pdf', bbox_inches='tight')
cg.plot(sd_beta_persamp, cmap='bone', min=0., max=0.15,
    llabel=r'\sigma_{\beta_{\mathrm{s}}}^{\mathrm{stat}}', width=4,
    extend='both', cbar=False)
plt.savefig('../figures/beta_n0016_sd_samp.pdf', bbox_inches='tight')
plt.close('all')

# Transparency



#sky_mean = mu + 3.1
#
#cg.plot((mu-mu.mean())/sd_beta, min=-6, max=6)
#plt.show()

#ret = cg.plot(mu, return_figure=True, return_only_data=True)
#plt.close('all')
#lon, lat, mu = ret[0]
#clip = 0.15/2
#sd_beta[sd_beta < clip] = clip
#weights = 1/sd_beta**2
#weights /= weights.max()
#cg.plot(weights, min=0, max=1)
#plt.show()
##weights = weights*0 + 1
#ret = cg.plot(weights, return_figure=True, return_only_data=True)
#lon, lat, alph = ret[0]
#
#fig = plt.figure(figsize=(5, 0.6*5))
#ax = fig.add_subplot(1,1,1, projection='mollweide')
#plt.pcolormesh(lon, lat, 0*mu+1, cmap='binary', vmin=0, vmax=1)
#
#cmap_path =  cmap_path = Path(cg.__path__[0]) / "data/planck_cmap.dat"
#cmap = ListedColormap(np.loadtxt(cmap_path) / 255.0, 'planck')
#plt.pcolormesh(lon, lat, mu, alpha=alph, cmap=cmap, vmin=-3.5, vmax=-2.5)
#
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticklabels([])
#ax.tick_params(axis="both", which="both", length=0)
#plt.savefig('test.pdf', bbox_inches='tight')
