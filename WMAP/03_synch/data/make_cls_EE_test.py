# Configure Matplotlib options
import healpy as hp
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as np
import sys

import pymaster as nmt
import healpy as hp

from scipy.optimize import minimize
from corner import corner

import emcee

l_fit_min = 2
l_fit_max = 140
l_plot_max = 140
l_plot_min = 2

# Using the same arguments as in Planck X 2015, taking power spectrum within
# bins.


def get_limits(ellb, Dl_avg, Dl_std, inds):
    res = minimize(chisq, np.array([-1, 1]), 
        args=(ellb[inds], Dl_avg[inds], Dl_std[inds]))
    alpha_synch, A_synch = res.x
    
    nwalkers = 32
    ndim = 2
    p0 = np.random.rand(nwalkers, ndim) + res.x
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[ellb[inds],
      Dl_avg[inds], Dl_std[inds]])
    state = sampler.run_mcmc(p0, 100)
    sampler.reset()

    sampler.run_mcmc(state, 100000)
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    meds = []
    for i in range(ndim):
      low, med, hi = np.percentile(flat_samples[:, i], [16, 50, 84])
      up = hi - med
      lo = med - low
      meds.append(med)
      print(f'param {i} = {np.round(med,3)} + {np.round(up,3)} - {np.round(lo,3)}')
    return meds, flat_samples[:,1]

def chisq(x, ell=0, Dl=0, sigmal=0, l_pivot=80):
  alpha, A = x
  model = A*(ell/l_pivot)**alpha
  return (((Dl - model)/sigmal)**2).sum()

def log_prob(x, ell=0, Dl=0, sigmal=0, l_pivot=80):
  alpha, A = x
  model = A*(ell/l_pivot)**alpha
  return -(((Dl - model)/sigmal)**2).sum()

def bin_spec(ell, Dl, l_ini=None, l_end=None):
  bins = []
  binv = []
  for i in range(len(l_ini)):
    bins.append(ell[(ell >= l_ini[i]) & (ell < l_end[i])])
    binv.append(Dl[(ell >= l_ini[i]) & (ell < l_end[i])])

  Dl_avg = np.array([Dli.mean() for Dli in binv])
  Dl_std = np.array([Dli.std()/len(Dli)**0.5 for Dli in binv])
  #ellb   = np.array([elli.prod()**(1/len(elli)) for elli in bins])
  ellb   = np.array([elli.mean() for elli in bins])
  return ellb, Dl_avg, Dl_std

def nmt_xspec(fname1, fname2, label=''):
  synch_cg11 = hp.ud_grade(hp.read_map(fname1,
      field=(0,1,2)), nside)
  synch_cg12 = hp.ud_grade(hp.read_map(fname2,
      field=(0,1,2)), nside)
  
  nmt_field1 = nmt.NmtField(mask, synch_cg11[1:])
  nmt_field2 = nmt.NmtField(mask, synch_cg12[1:])
  
  Clhat_cg1 = nmt.compute_full_master(nmt_field1, nmt_field2, b)
  EE = Clhat_cg1[0]
  BB = Clhat_cg1[3]
  np.save(f'ee{label}_spec.npy', EE)
  np.save(f'bb{label}_spec.npy', BB)
  return EE, BB
def xpol_xspec(fname1, fname2, label=''):
  synch_cg1 = hp.ud_grade(hp.read_map(fname1,
      field=(0,1,2)), nside)
  synch_cg2 = hp.ud_grade(hp.read_map(fname2,
      field=(0,1,2)), nside)

  pcl, cl = xp.get_spectra(synch_cg1, synch_cg2, Dl=False)
  return cl[1], cl[2]


nside = 1024

# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map(
        "COM_Mask_CMB-common-Mask-Pol_2048_R3.00.fits"
        )
mask = hp.ud_grade(mask, nside)

msk_apo = nmt.mask_apodization(mask, 1.0, apotype='C1')
hp.mollview(msk_apo)
hp.mollview(mask)
mask = msk_apo
fsky = sum(mask)/len(mask)

ee_spectra = np.loadtxt('ee_spectra.txt')

l_ini = np.concatenate((
  np.array([2,5,12]),
  #np.array([4,12]),
  np.arange(20, 520, 20),
  np.array([550])))
l_end = np.concatenate((
  np.array([4,11,19]),
  #np.array([11,19]),
  np.arange(40, 520, 20) - 1,
  np.array([549, 599])))

ells = np.arange(2, 600)
sigma = np.sqrt(2/fsky/(2*ells+1))
sigma_bin = np.zeros(len(l_ini))


# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map(
        "COM_Mask_CMB-common-Mask-Pol_2048_R3.00.fits"
        )
mask = hp.ud_grade(mask, nside)


b = nmt.NmtBin.from_nside_linear(nside, 1)
ell_eff = b.get_effective_ells()
Z = ell_eff*(ell_eff+1)/(2*np.pi)

try:
  EE = np.load('ee_spec.npy')
  BB = np.load('bb_spec.npy')
except IOError:
  fname1 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM1_longer/CG_synch_IQU_n1024_CG_HM1_longer.fits'
  fname2 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM2_longer/CG_synch_IQU_n1024_CG_HM2_longer.fits'
  EE, BB = nmt_xspec(fname1, fname2)

try:
  EE_K = np.load('ee_K_spec.npy')*1e6 * ((23/30)**3.1)**2
  BB_K = np.load('bb_K_spec.npy')*1e6 * ((23/30)**3.1)**2
except IOError:
  fname1 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM1_longer/CG_023-WMAP_K_IQU_n0512_CG_HM1_longer.fits'
  fname2 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM2_longer/CG_023-WMAP_K_IQU_n0512_CG_HM2_longer.fits'
  EE_K, BB_K = nmt_xspec(fname1, fname2, label='_K')
  EE_K *= 1e6 * ((23/30)**3.1)**2
  BB_K *= 1e6 * ((23/30)**3.1)**2
try:
  EE_30 = np.load('ee_30_spec.npy') * ((28/30)**3.1)**2
  BB_30 = np.load('bb_30_spec.npy') * ((28/30)**3.1)**2
except IOError:
  fname1 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM1_longer/CG_030_IQU_n0512_CG_HM1_longer.fits'
  fname2 = '/mn/stornext/d5/data/duncanwa/WMAP/CG_HM2_longer/CG_030_IQU_n0512_CG_HM2_longer.fits'
  EE_30, BB_30 = nmt_xspec(fname1, fname2, label='_30')
  EE_K *= 1e6 * ((28/30)**3.1)**2
  BB_K *= 1e6 * ((28/30)**3.1)**2





for i in range(len(l_ini)):
    sigma_in_bin = sigma[(ells >= l_ini[i]) & (ells <= l_end[i])]
    sigma_bin[i] = sum(sigma_in_bin**-2)**-0.5

ell_eff = (l_end + l_ini - 1) / 2
Z = ell_eff*(ell_eff+1)/(2*np.pi)

def a2t(nu, T_cmb = 2.7255):
    "Returns conversion factor between antenna and thermodynamic units"
    h     = 1.0545726691251021e-34 * 2.0*np.pi
    k_b   = 1.3806503e-23
    x     = h*nu/(k_b*T_cmb) * 1e9
    return (np.exp(x)-1.e0)**2 / (x**2 * np.exp(x))

width = 8.8

# Load data
vmin = -110
vmax =  160

#cls_cmb_r0 = np.loadtxt('base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl')
cls_cmb_r0 = np.loadtxt('COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt')
cls_cmb_r2 = np.loadtxt('camb_tau0.06_r0.20_Aprior.dat')

#cls_dust  = np.loadtxt('cl_dust_hm_reformat.dat')
#cls_synch = np.loadtxt('cl_synch_hm_reformat.dat')

#cls_dust  = np.loadtxt('cls_commander_dx12_pol_dust_n2048_v3_hm_EE.dat')
#cls_synch = np.loadtxt('cls_commander_dx12_pol_synch_n2048_v3_hm_EE.dat')
#ls = cls_dust[:,0]
#cls_dust[:,1] = cls_dust[:,1]*ls*(ls+1.)/2./np.pi
#cls_dust[:,1] = cls_dust[:,1]/
#np.exp(-ls*(ls+1)*pow(10.*np.pi/180./60/np.sqrt(8.*np.log(2.)),2))
#cls_synch[:,1] = cls_synch[:,1]*ls*(ls+1.)/2./np.pi
#cls_synch[:,1] = cls_synch[:,1]/
#np.exp(-ls*(ls+1)*pow(40.*np.pi/180./60/np.sqrt(8.*np.log(2.)),2))

scale_synch = (a2t(30.))**2 
scale_dust  = (a2t(353.))**2

scale_all = 1 / 1.2037
scale_comm = 1.1114**2
cls_dust_comm_EE = np.loadtxt('cl_commander_thermaldust_spectra_EE_LR78.txt') #* scale_dust
cls_dust_smica_EE = np.loadtxt('cl_smica_thermaldust_spectra_EE_LR78.txt') #* scale_dust
cls_dust_353_EE = np.loadtxt('cl_353GHz_thermaldust_spectra_EE_LR78.txt') #* scale_dust

cls_synch_comm_EE = np.loadtxt('cl_commander_synchrotron_spectra_EE_LR78.txt') #* scale_dust
cls_synch_smica_EE = np.loadtxt('cl_smica_synchrotron_spectra_EE_LR78.txt')


lmin   = 2
lmax   = 800

dl  = 10


        

A_dust      = 389.*scale_all
alpha_dust  = -0.40 # 0.01
A_synch     = 2.3
alpha_synch = -0.84 # 0.1

b = nmt.NmtBin.from_nside_linear(nside, 1)
ell_eff = b.get_effective_ells()
Z = ell_eff*(ell_eff+1)/(2*np.pi)

# Create the plot
fig = plt.figure(figsize=(cm2inch(width), 6./8.*cm2inch(width)))
# this should be changed for making a panel of multiple figures
ax = fig.add_subplot(111)
ax.set_xscale('log')
#ax.set_yscale('log', nonposy='clip')
ax.set_yscale('log')




b = nmt.NmtBin.from_nside_linear(nside, 1)
ell_eff = b.get_effective_ells()
Z = ell_eff*(ell_eff+1)/(2*np.pi)

ellb, Dl_avg, Dl_std = bin_spec(ell_eff, Z*EE, l_ini=l_ini, l_end=l_end)
Zb = ellb*(ellb+1)/(2*np.pi)
inds = (ellb <= l_fit_max) & (ellb >= l_fit_min)
(alpha_synch, A_synch), E_arr = get_limits(ellb, Dl_avg, Dl_std, inds)

inds = (ellb <= l_plot_max) & (ellb >= l_plot_min)

plt.errorbar(ellb[inds], Dl_avg[inds], Dl_std[inds],
    linestyle='',
    marker='o', ms=2, elinewidth=0.7, capsize=1, label=r'Cosmoglobe',
    color='k')
ellb, Dl_avg, Dl_std = bin_spec(ell_eff, Z*EE_K, l_ini=l_ini, l_end=l_end)
#plt.errorbar(ellb[inds]*0.9, Zb[inds]*Cl_avg[inds], Zb[inds]*Cl_std[inds],
#    linestyle='',
#    marker='o', ms=2, elinewidth=0.7, capsize=1,
#    label='K',
#    color='C1')
#ellb, Cl_avg, Cl_std = bin_spec(ell_eff, EE_30, l_ini=l_ini, l_end=l_end)
#plt.errorbar(ellb[inds]*1.1, Zb[inds]*Cl_avg[inds], Zb[inds]*Cl_std[inds],
#    linestyle='',
#    marker='o', ms=2, elinewidth=0.7, capsize=1,
#    label='30',
#    color='C2')


ellb, Dl_avg, Dl_std = bin_spec(ell_eff, Z*BB, l_ini=l_ini, l_end=l_end)
plt.errorbar(ellb[inds], Dl_avg[inds], Dl_std[inds],
    linestyle='',
    marker='o', ms=2, elinewidth=0.7, capsize=1, 
    color='k', markerfacecolor='w')

inds = (ellb <= l_fit_max) & (ellb >= l_fit_min)

(alpha_synch_BB, A_synch_BB), B_arr = get_limits(ellb, Dl_avg, Dl_std, inds)

lo, med, hi = np.percentile(B_arr/E_arr, [16, 50, 84])
print('Cosmoglobe results')
print(med, hi-med, med-lo)

plt.plot(cls_cmb_r0[:,0], cls_cmb_r0[:,3], label=r"Best-fit $\Lambda$CDM $\mathcal D_\ell^{\mathrm{EE}}$",
    color='black', linewidth=0.5, linestyle='-')

l_pivot     = 80
fit = A_synch * (cls_cmb_r0[:,0]/l_pivot)**alpha_synch
plt.plot(cls_cmb_r0[:,0], fit, color='k', linewidth=1, label='E-modes')

l_pivot     = 80
fit = A_synch_BB * (cls_cmb_r0[:,0]/l_pivot)**alpha_synch_BB
plt.plot(cls_cmb_r0[:,0], fit, color='k', linestyle=':', linewidth=1,
label='B-modes')

EE_18 = np.loadtxt('ee_spectra_unbinned_2018.txt')[0]
BB_18 = np.loadtxt('bb_spectra_unbinned_2018.txt')[0]

ellb, Dl_avg, Dl_std = bin_spec(ell_eff, Z*EE_18, l_ini=l_ini, l_end=l_end)
print(ellb)
inds = (ellb <= l_plot_max) & (ellb >= l_plot_min)
plt.errorbar(ellb[inds]*1.1, Dl_avg[inds], Dl_std[inds],
    linestyle='',
    marker='o', ms=2, elinewidth=0.7, capsize=1,
    color='C0', label=r'PR3')
(alpha_EE, A_EE), E_arr = get_limits(ellb, Dl_avg, Dl_std, inds)

ellb, Dl_avg, Dl_std = bin_spec(ell_eff, Z*BB_18, l_ini=l_ini, l_end=l_end)

plt.errorbar(ellb[inds]*1.1, Dl_avg[inds], Dl_std[inds],
    linestyle='',
    marker='o', ms=2, elinewidth=0.7, capsize=1,
    color='C0', markerfacecolor='w')

inds = (ellb <= l_fit_max) & (ellb >= l_fit_min)
(alpha_BB, A_BB), B_arr = get_limits(ellb, Dl_avg, Dl_std, inds)

lo, med, hi = np.percentile(B_arr/E_arr, [16, 50, 84])
print('PR3 results')
print(med, hi-med, med-lo)

fit = A_EE * (cls_cmb_r0[:,0]/l_pivot)**alpha_EE
plt.plot(cls_cmb_r0[:,0], fit, color='C0', linewidth=1)

l_pivot     = 80
fit = A_BB * (cls_cmb_r0[:,0]/l_pivot)**alpha_BB
plt.plot(cls_cmb_r0[:,0], fit, color='C0', linestyle=':', 
    linewidth=1)


'''
l, Dl, sigmal = np.loadtxt('cl_commander_synchrotron_spectra_EE_LR78.txt').T
inds = (l < 140)
plt.errorbar(l[inds], Dl[inds], sigmal[inds], fmt='.', color='r')
l, Dl, sigmal = np.loadtxt('cl_commander_synchrotron_spectra_BB_LR78.txt').T
plt.errorbar(l[inds], Dl[inds], sigmal[inds], fmt='.', color='r',
        markerfacecolor='w')
'''

# legend
leg = plt.legend(frameon=False, loc=3, fontsize=8)

# labels
plt.xlabel(r"Multipole moment, $\ell$"); plt.ylabel(r"Power spectrum, $\mathcal{D}_{\ell}$ $[\mu\mathrm{K}^2]$")
ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

  
# grid
plt.grid(False, which="major", axis="both")

# axes limits
plt.xlim(1.8, 150); plt.ylim(ymin=0.3)

plt.xticks([3,10,30, 100], [3, 10, r"$30$", r"$100$"])
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(2.2, 400, r"$EE$", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

    
plt.savefig("../figures/cls_synch_ratio.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
