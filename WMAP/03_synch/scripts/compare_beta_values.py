import numpy as np
import healpy as hp
import cosmoglobe as cg

import matplotlib.pyplot as plt

from glob import glob

from tqdm import tqdm

from cycler import cycler
default_cycler = (cycler(color=[plt.cm.cividis(i) for i in np.linspace(0,1,5)]))

ms=2.5
elinewidth=1

jitter = 0.12

DIR = '/mn/stornext/d5/data/duncanwa/comm1'

#plt.rc('lines', linewidth=4)
#plt.rc('axes', prop_cycle=default_cycler)

def get_mu_sd(fnames, mu_map=False):
    M2n = np.zeros(12*64**2)
    mu = np.zeros(12*64**2)
    mu_old = np.zeros(12*64**2)
    
    N = 0
    for i in tqdm(range(len(fnames))):
        N += 1
        m = hp.read_map(fnames[i])
        mu = mu_old + (m - mu_old)/N
        M2n = M2n + (m - mu_old)*(m - mu)
        mu_old = mu
    
    var = M2n/(N-1)
    
    mu_vals = np.zeros(24)
    sd_vals = np.zeros(24)
    
    
    for i in reg_inds:
        inds = (regs == i)
        mu_vals[i-1] = mu[inds].mean()
        sd_vals[i-1] = var[inds].mean()**0.5
    if mu_map:
        return mu_vals, sd_vals, mu, var**0.5
    else:
        return mu_vals, sd_vals

reg_inds = np.arange(1, 25)

regs = hp.read_map(f'{DIR}/UF_sindex_regions_full.fits')


fnames = glob(f'{DIR}/cg_chain_batmask_000024_0/synch_beta*k*.fits')

'''
mu, sd = get_mu_sd(fnames)
plt.errorbar(reg_inds, mu, sd, fmt='.')

fnames = glob('cg_chain_batmask_??????_?/synch_beta*k*100.fits')

mu, sd = get_mu_sd(fnames)
plt.errorbar(reg_inds, mu, sd, fmt='.')
'''


fnames = glob(f'{DIR}/cg_chain_batmask_-2.7_000024_0/synch_beta*.fits')
mu1, sd1 = get_mu_sd(fnames)
fnames = glob(f'{DIR}/cg_chain_batmask_-3.5_000024_0/synch_beta*.fits')
mu2, sd2 = get_mu_sd(fnames)
'''
plt.figure()
plt.errorbar(reg_inds, mu1, sd1, fmt='.')
plt.errorbar(reg_inds, mu2, sd2, fmt='.')

plt.figure()
plt.errorbar(reg_inds, mu1 - mu2, np.sqrt(sd1**2+sd2**2), fmt='.')
plt.ylabel(r'$\beta_s$ with $P(-2.7) - P(-3.5)$')
plt.xlabel('region')
'''



fnames = glob(f'{DIR}/cg_chain_n16_??????_?/synch_beta*k*[1-9]?.fits')
mu, sd, mu_16, sd_16= get_mu_sd(fnames, mu_map=True)
#plt.errorbar(reg_inds, mu, sd, fmt='.', label=r'$N_\mathrm{side}=16$')

fnames = glob(f'{DIR}/cg_chain_batmask_??????_?/synch_beta*k*?[1-9]?.fits')
#fnames = glob('cg_chain_batmask_??????_?/synch_beta*k*?[1-9]?.fits')

mu, sd, mu_bat, sd_bat = get_mu_sd(fnames, mu_map=True)
#plt.errorbar(reg_inds, mu, sd, fmt='.')


plt.figure(figsize=(8, 4))
sd_prior = (mu1 - mu2)/4
sd_tot = np.hypot(sd, sd_prior)
plt.errorbar(reg_inds-2*jitter, mu, sd_tot, fmt='.', label=r'Commander1 WMAP+LFI',
        ms=ms, elinewidth=elinewidth)


d_K30 = np.loadtxt('/mn/stornext/d5/unnif/sindex_bp/coswmap23_cos30_500s/combab_011-250/ut_sample_betas_invvar.txt')
d_KKa = np.loadtxt('/mn/stornext/d5/unnif/sindex_bp/coswmap23_coswmap33_500s/combab_011-250/ut_sample_betas_invvar.txt')

plt.errorbar(d_K30[:,0]-jitter, d_K30[:,1], d_K30[:,2], fmt='.', label='TT K/30', color='C3', ms=ms, elinewidth=elinewidth)
plt.errorbar(d_KKa[:,0], d_KKa[:,1], d_KKa[:,2], fmt='.', label='TT K/Ka',
        color='C4', ms=ms, elinewidth=elinewidth)


beta_class = hp.ma(hp.read_map('/mn/stornext/d16/cmbco/ola/class/class_dr1_40GHz_beta_s_nside32-d0.fits'))
var_class = hp.ma(hp.read_map('/mn/stornext/d16/cmbco/ola/class/class_dr1_40GHz_beta_s_nside32-d0.fits', field=1))**2

mu_vals = np.zeros(24)
sd_vals = np.zeros(24)

regs = hp.ud_grade(regs, 32)
#hp.mollview(regs)
#hp.mollview(beta_class)

print(regs.size)
print(beta_class.size)

alphas = []
for i in reg_inds:
    inds = (regs == i)
    ok = (beta_class[inds] != hp.UNSEEN)
    alphas.append(ok.sum()/inds.sum())
    mu_vals[i-1] = beta_class[inds].mean()
    sd_vals[i-1] = np.hypot(var_class[inds].mean()**0.5/ok.sum()**0.5, beta_class[inds].std())
    print(ok.sum(), var_class[inds].mean()**0.5)
    #sd_vals[i-1] = var_class[inds].mean()**0.5
    if ok.sum() > 0.5:
        plt.errorbar(i+jitter, mu_vals[i-1], sd_vals[i-1], color='C1', fmt='.', #label='CLASS',
                #alpha=alphas[-1])
                ms=ms, elinewidth=elinewidth)
plt.errorbar([], [], [], color='C1', fmt='.', label='CLASS 40 + K', ms=ms, elinewidth=elinewidth)


beta_QUI, sigma_QUI = hp.read_map('/mn/stornext/d16/cmbco/ola/quijote/compsep/pol//betas_quijote_mfi_cs_pol_64_dr1.fits', field=(0,1))

beta_QUI = hp.ma(beta_QUI)
sigma_QUI = hp.ma(sigma_QUI)

#print(hp.npix2nside(len(beta_QUI)))
#cg.plot(sigma_QUI,min=0, max=0.4)
#plt.show()

regs = hp.read_map(f'{DIR}/UF_sindex_regions_full.fits')
regs = hp.ud_grade(regs, 64)

alphas = []
for i in reg_inds:
    inds = (regs == i)
    ok = (beta_QUI[inds] != hp.UNSEEN)
    alphas.append(ok.sum()/inds.sum())
    mu_vals[i-1] = beta_QUI[inds].mean()
    sd_vals[i-1] = np.hypot((sigma_QUI[inds]**2).mean()**0.5/ok.sum()**0.5, beta_QUI[inds].std())
    #sd_vals[i-1] = var_class[inds].mean()**0.5
    print(ok.sum(), (sigma_QUI[inds]**2).mean()**0.5)
    if ok.sum() > 0.5:
        plt.errorbar(i+2*jitter, mu_vals[i-1], sd_vals[i-1], color='C2', fmt='.', #label='CLASS',
                #alpha=alphas[-1])
                ms=ms, elinewidth=elinewidth)
plt.errorbar([], [], [], color='C2', fmt='.', label='MFI+K/Ka+PR4', ms=ms,
        elinewidth=elinewidth)


plt.legend(loc='best')


plt.ylim([-4.75, -2.1])
plt.xlabel('Region number')
plt.ylabel(r'Spectral index, $\beta_{\mathrm{s}}$')
plt.xticks([1,5,10,15,20,24])
plt.xticks(np.arange(1,25), minor=True)
plt.yticks(np.arange(-4.7, -2.1, 0.1), minor=True)
plt.savefig('../figures/compare_betas.pdf', bbox_inches='tight')

'''
cg.plot(mu_16, min=-3.5, max=-2.5)
cg.plot(beta_QUI, min=-3.5, max=-2.5)

cg.plot('/mn/stornext/d16/cmbco/ola/class/class_dr1_40GHz_beta_s_nside32-d0.fits',
    min=-3.5, max=-2.5)
'''
#
#
#cg.plot(sd_bat)
#cg.plot(sd_16)
#cg.plot('/mn/stornext/d16/cmbco/ola/class/class_dr1_40GHz_beta_s_nside32-d0.fits',
#    sig=1)
#
#beta_class = hp.read_map('/mn/stornext/d16/cmbco/ola/class/class_dr1_40GHz_beta_s_nside32-d0.fits')
#beta_class = hp.ud_grade(beta_class, 8)
#beta_16 = hp.ud_grade(mu_16, 8)
#cg.plot(beta_class - beta_16, min=-0.5, max=0.5, llabel=r'\beta^\mathrm{CLASS}-\beta^\mathrm{CG}')




#plt.show()
