import numpy as np
import healpy as hp
import cosmoglobe as cg

import matplotlib.pyplot as plt

from matplotlib import rcParams, rc

from setup_matplotlib import *

# common setup for matplotlib
params = {'backend': 'pdf',
          'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 10,
          }

def cm2inch(cm):
    """Centimeters to inches"""
    return cm *0.393701

from matplotlib.ticker import MaxNLocator
width = 8.8
# Create the plot
fig = plt.figure(figsize=(1.25*cm2inch(width), 6./8.*cm2inch(width)))
ax = fig.add_subplot(111)

y_est=-3.25796758275099
y_min=-3.33676787269876
y_max=-3.18442467987845
#plt.axhline(y=y_est, linestyle='-', linewidth=0.8)
plt.plot((0.8,12.2),(y_est,y_est), label=r"Planck 2018 likelihood", linestyle='-', linewidth=0.8)
plt.fill_between(np.linspace(0.8,12.2,11), y_min, y_max, alpha=0.3)


marker='x'
linewidth=0.5
markersize=3
capsize=2
ls='none'
alpha=0.4



from glob import glob

from tqdm import tqdm

from cycler import cycler
default_cycler = (cycler(color=[plt.cm.cividis(i) for i in np.linspace(0,1,5)]))

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



fnames = glob(f'{DIR}/cg_chain_n16_????[2,4,6,8,0]0_1/synch_beta*k*[1-9]?.fits')
mu, sd, mu_16, sd_16= get_mu_sd(fnames, mu_map=True)
#plt.errorbar(reg_inds, mu, sd, fmt='.', label=r'$N_\mathrm{side}=16$')

fnames = glob(f'{DIR}/cg_chain_batmask_????[2,4,6,8,0]0_1/synch_beta*k*?[1-9]?.fits')
#fnames = glob('cg_chain_batmask_??????_?/synch_beta*k*?[1-9]?.fits')

mu, sd, mu_bat, sd_bat = get_mu_sd(fnames, mu_map=True)
#plt.errorbar(reg_inds, mu, sd, fmt='.')


#plt.figure(figsize=(8, 4))
sd_prior = (mu1 - mu2)/4
sd_tot = np.hypot(sd, sd_prior)


d_K30 = np.loadtxt('/mn/stornext/d5/unnif/sindex_bp/coswmap23_cos30_500s/combab_011-250/ut_sample_betas_invvar.txt')
d_KKa = np.loadtxt('/mn/stornext/d5/unnif/sindex_bp/coswmap23_coswmap33_500s/combab_011-250/ut_sample_betas_invvar.txt')

#plt.errorbar(d_K30[:,0]-jitter, d_K30[:,1], d_K30[:,2], fmt='.', label='TT K/30', color='C3', ms=ms, elinewidth=elinewidth)
plt.errorbar(d_KKa[:,0]-jitter, d_KKa[:,1], d_KKa[:,2], label='Cosmoglobe K/Ka TT',
        color='red', ms=markersize, elinewidth=elinewidth, linewidth=linewidth,
        capsize=capsize, marker=marker, ls=ls, alpha=0.4)


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
        plt.errorbar(i, mu_vals[i-1], sd_vals[i-1], color='blue',
                alpha=0.4,
                ms=markersize, elinewidth=elinewidth, linewidth=linewidth,
                capsize=capsize, marker=marker, ls=ls)
plt.errorbar([-5], [-5], [2], color='blue', label='CLASS',
        ms=markersize, elinewidth=elinewidth, linewidth=linewidth,
        capsize=capsize, marker=marker, ls=ls, alpha=0.4)


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
        plt.errorbar(i+jitter, mu_vals[i-1], sd_vals[i-1], color='brown',
                #alpha=alphas[-1])
                alpha=0.4,
                ms=markersize, elinewidth=elinewidth, linewidth=linewidth,
                capsize=capsize, marker=marker, ls=ls)
plt.errorbar([-5], [-5], [0], color='brown', label='QUIJOTE', ms=markersize,
        elinewidth=elinewidth, linewidth=linewidth, capsize=capsize,
        marker=marker, ls=ls, alpha=0.4)

plt.errorbar(reg_inds+2*jitter, mu, sd_tot, label=r'Cosmoglobe DR1',
        ms=markersize, elinewidth=elinewidth, color='k', linewidth=linewidth,
        capsize=capsize, marker=marker, ls=ls)


plt.ylim([-4.1, -2.3])
plt.xlim(0.2,24.8)
plt.xticks([1,5,10,15,20,24])
plt.yticks([-4, -3])
plt.xticks(np.arange(1,25), minor=True)
plt.yticks(np.arange(-4.1, -2.3, 0.1), minor=True)

ax.minorticks_on()
plt.ylabel(r"Spectral index, $\beta_\mathrm{s}$", fontsize=10);
plt.xlabel(r"Region number", fontsize=10);
leg = plt.legend(frameon=True, loc='lower right', prop={'size':7})
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(.0)


plt.savefig('../figures/compare_betas.pdf', bbox_inches='tight',
        bbox_extra_artists=[],pad_inches=0.03)

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
