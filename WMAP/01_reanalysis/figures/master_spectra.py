import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cosmoglobe as cg
import astropy.units as u

from glob import glob

# Import the NaMaster python wrapper
import pymaster as nmt

#  Simple example showcasing the use of NaMaster to compute the pseudo-Cl
#  estimator of the angular cross-power spectrum of a spin-0 field and a
#  spin-2 field

# HEALPix resolution parameter used here
nside = 512

# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map("/mn/stornext/d16/cmbco/bp/dwatts/WMAP/data_WMAP/wmap_kq75_TQU_mask_r9.fits")


# bin4 = nmt.NmtBin.from_edges(l_ini, l_end)
logbins = np.geomspace(10, 3*512-1, 50).astype('int')
l_ini = np.concatenate((np.arange(2, 10),     logbins[:-1]))
l_end = np.concatenate((np.arange(2, 10) + 1, logbins[1:]))
#b = nmt.NmtBin.from_nside_linear(nside, 1)
b = nmt.NmtBin.from_edges(l_ini, l_end)
ell_eff = b.get_effective_ells()

DIR = '/mn/stornext/d16/cmbco/ola/wmap/freq_maps'

wmap_maps = [
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_K1_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Ka1_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q2_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V1_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V2_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W1_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W2_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W3_v5.fits',
'/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W4_v5.fits']
CG_DIR = '/mn/stornext/d5/data/duncanwa/WMAP'

cg_maps = [
 f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_023-WMAP_K_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_030-WMAP_Ka_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_040-WMAP_Q1_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_040-WMAP_Q2_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_060-WMAP_V1_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_060-WMAP_V2_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_090-WMAP_W1_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_090-WMAP_W2_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_090-WMAP_W3_map_c0001_k000002.fits',
f'{CG_DIR}/chains_CG_LFI_KKaQVW_c_230103/tod_090-WMAP_W4_map_c0001_k000002.fits']

bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

x,y,z = hp.pix2vec(512, np.arange(12*512**2))
dip_W = -0.233*x -2.222*y + 2.504*z


fig1, axes1 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig2, axes2 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig3, axes3 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)

axs1 = axes1.flatten()
axs2 = axes2.flatten()
axs3 = axes3.flatten()

for i in range(len(cg_maps)):
    print(i)
    m_WMAP = hp.read_map(wmap_maps[i])*1e3
    m_CG = hp.read_map(cg_maps[i])*1e3
    m_CG -= dip_W*1e3
    mono_CG = hp.fit_monopole(m_CG, gal_cut=50)
    mono_WM = hp.fit_monopole(m_WMAP, gal_cut=50)
    m_WMAP = m_WMAP - mono_WM + mono_CG
    f_WMAP = nmt.NmtField(mask, [m_WMAP])
    f_CG = nmt.NmtField(mask, [m_CG])
    f_diff = nmt.NmtField(mask, [m_CG - m_WMAP])
    Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)[0]
    print('ps1')
    Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)[0]
    print('ps2')
    #Clhat_diff = nmt.compute_full_master(f_diff, f_diff, b)[0]
    ell = np.arange(len(Clhat_C))
    print('got the power spectra')


    #l1, = axs1[i].loglog(ell[2:], Clhat_W[2:], label='WMAP')
    #l2, = axs1[i].loglog(ell[2:], Clhat_C[2:], label='Cosmoglobe')
    l1, = axs1[i].loglog(ell_eff, Clhat_W, label='WMAP')
    l2, = axs1[i].loglog(ell_eff, Clhat_C, label='Cosmoglobe')
    axs1[i].text(0.75, 0.75, r'\textit{'+bands[i]+'}', transform=axs1[i].transAxes)
    #plt.ylabel(r'$C_\ell^{TT}\ [\mathrm{\mu K}^2]$')

    inds = (ell_eff > 220)
    #axs2[i].plot(ell[inds], Clhat_W[inds], label='WMAP')
    #axs2[i].plot(ell[inds], Clhat_C[inds], label='Cosmoglobe')
    axs2[i].plot(ell_eff[inds], Clhat_W[inds], label='WMAP')
    axs2[i].plot(ell_eff[inds], Clhat_C[inds], label='Cosmoglobe')
    axs2[i].text(0.75, 0.75,  r'\textit{'+bands[i]+'}', transform=axs2[i].transAxes)
    #axs2[i].set_xlim([225, 1400])
    axs2[i].set_ylim([0.002, 0.1])
    #plt.legend()
    #plt.title(bands[i])
    #plt.ylabel(r'$C_\ell^{TT}\ [\mathrm{\mu K}^2]$')
    #plt.xlabel(r'$\ell$')
    #plt.savefig(f'{bands[i]}_TT_zoom.png', bbox_inches='tight')

    axs3[i].semilogx(ell_eff, Clhat_W/Clhat_C)
    #axs3[i].set_ylabel(r'$C_\ell^\mathit{WMAP}/C_\ell^\mathrm{Cosmoglobe}$')
    #axs3[i].set_xlabel(r'$\ell$')
    axs3[i].text(0.25, 0.75, r'\textit{'+bands[i]+'}', transform=axs3[i].transAxes)
    #plt.title(bands[i])
    #axs3[3].set_ylim([0.8, 1.4])
    #plt.savefig(f'{bands[i]}_TT_ratio.png', bbox_inches='tight')
    #plt.close('all')

axs1[10].axis('off')
axs1[11].axis('off')
axs2[10].axis('off')
axs2[11].axis('off')
axs3[10].axis('off')
axs3[11].axis('off')

axs1[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])
axs2[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])

fig1.supxlabel(r'$\ell$')
fig2.supxlabel(r'$\ell$')
fig3.supxlabel(r'$\ell$')

fig1.supylabel(r'$C_\ell^\mathrm{TT}\ [\mathrm{\mu K}^2]$')
fig2.supylabel(r'$C_\ell^\mathrm{TT}\ [\mathrm{\mu K}^2]$')
fig3.supylabel(r'$C_\ell^\mathit{WMAP}/C_\ell^\mathrm{Cosmoglobe}$')

plt.figure(fig1.number)
plt.savefig(f'TT_spectra.pdf', bbox_inches='tight')
plt.figure(fig2.number)
plt.savefig(f'TT_spectra_zoom.pdf', bbox_inches='tight')
plt.figure(fig3.number)
plt.savefig(f'TT_spectra_ratio.pdf', bbox_inches='tight')
print('here we are')
plt.close('all')



fig1, axes1 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig2, axes2 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig3, axes3 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig4, axes4 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig5, axes5 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig6, axes6 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)

axs1 = axes1.flatten()
axs2 = axes2.flatten()
axs3 = axes3.flatten()
axs4 = axes4.flatten()
axs5 = axes5.flatten()
axs6 = axes6.flatten()


for i in range(len(cg_maps)):
    print(i)
    m_WMAP = hp.read_map(wmap_maps[i], field=(1,2))*1e3
    m_CG = hp.read_map(cg_maps[i], field=(1,2))*1e3
    f_WMAP = nmt.NmtField(mask, m_WMAP)
    f_CG = nmt.NmtField(mask, m_CG)
    Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)
    Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)
    ell = np.arange(len(Clhat_C[0]))

    l1, = axs1[i].loglog(ell_eff, Clhat_W[0], label='WMAP')
    l2, = axs1[i].loglog(ell_eff, Clhat_C[0], label='Cosmoglobe')
    axs2[i].loglog(ell_eff, Clhat_W[3], label='WMAP')
    axs2[i].loglog(ell_eff, Clhat_C[3], label='Cosmoglobe')

    inds = (ell > 400)
    axs3[i].plot(ell_eff[inds], Clhat_W[0][inds], label='WMAP')
    axs3[i].plot(ell_eff[inds], Clhat_C[0][inds], label='Cosmoglobe')
    axs4[i].plot(ell_eff[inds], Clhat_W[3][inds], label='WMAP')
    axs4[i].plot(ell_eff[inds], Clhat_C[3][inds], label='Cosmoglobe')


    axs5[i].semilogx(ell[2:], Clhat_W[0][2:]/Clhat_C[0][2:])
    axs6[i].semilogx(ell[2:], Clhat_W[3][2:]/Clhat_C[3][2:])

axs1[10].axis('off')
axs1[11].axis('off')
axs2[10].axis('off')
axs2[11].axis('off')
axs3[10].axis('off')
axs3[11].axis('off')
axs4[10].axis('off')
axs4[11].axis('off')
axs5[10].axis('off')
axs5[11].axis('off')
axs6[10].axis('off')
axs6[11].axis('off')

axs1[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])
axs2[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])
axs3[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])
axs4[10].legend(handles=[l1, l2], 
    labels=[r'\textit{WMAP}', r'\textsc{Cosmoglobe}'])

fig1.supxlabel(r'$\ell$')
fig2.supxlabel(r'$\ell$')
fig3.supxlabel(r'$\ell$')
fig4.supxlabel(r'$\ell$')
fig5.supxlabel(r'$\ell$')
fig6.supxlabel(r'$\ell$')

fig1.supylabel(r'$C_\ell^\mathrm{EE}\ [\mathrm{\mu K}^2]$')
fig2.supylabel(r'$C_\ell^\mathrm{BB}\ [\mathrm{\mu K}^2]$')
fig3.supylabel(r'$C_\ell^\mathrm{EE}\ [\mathrm{\mu K}^2]$')
fig4.supylabel(r'$C_\ell^\mathrm{BB}\ [\mathrm{\mu K}^2]$')
fig5.supylabel(r'$C_\ell^\mathit{WMAP}/C_\ell^\mathrm{Cosmoglobe}$')
fig6.supylabel(r'$C_\ell^\mathit{WMAP}/C_\ell^\mathrm{Cosmoglobe}$')

plt.figure(fig1.number)
plt.savefig(f'EE_spectra.pdf', bbox_inches='tight')
plt.figure(fig2.number)
plt.savefig(f'BB_spectra.pdf', bbox_inches='tight')
plt.figure(fig3.number)
plt.savefig(f'EE_spectra_zoom.pdf', bbox_inches='tight')
plt.figure(fig1.number)
plt.savefig(f'BB_spectra_zoom.pdf', bbox_inches='tight')
plt.figure(fig2.number)
plt.savefig(f'EE_spectra_ratio.pdf', bbox_inches='tight')
plt.figure(fig3.number)
plt.savefig(f'BB_spectra_ratio.pdf', bbox_inches='tight')
'''
fnames1 = glob(f'{DIR}/wmap_iqusmap_r9_yr?_W1_v5.fits')
fnames2 = glob(f'{DIR}/wmap_iqusmap_r9_yr?_W2_v5.fits')
fnames3 = glob(f'{DIR}/wmap_iqusmap_r9_yr?_W3_v5.fits')
fnames4 = glob(f'{DIR}/wmap_iqusmap_r9_yr?_W4_v5.fits')

fnames = fnames1 + fnames2 + fnames3
f2s = []
for f in fnames:
  print(f)
  f2s.append(nmt.NmtField(mask, hp.read_map(f, field=(1,2))))

Cls = []

for i in range(len(f2s)):
      for j in range(i+1,len(f2s)):
                print(i,j)
                cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
                Cls.append(cl_22)

Cls = np.array(Cls)

np.save('fnames_all', Cls)

W1 = hp.read_map(f'{DIR}/wmap_iqusmap_r9_9yr_W1_v5.fits', field=(1,2))
W2 = hp.read_map(f'{DIR}/wmap_iqusmap_r9_9yr_W2_v5.fits', field=(1,2))
W3 = hp.read_map(f'{DIR}/wmap_iqusmap_r9_9yr_W3_v5.fits', field=(1,2))
W4 = hp.read_map(f'{DIR}/wmap_iqusmap_r9_9yr_W4_v5.fits', field=(1,2))

Ws = [W1, W2, W3, W4]
f2s = []
for f in Ws:
  f2s.append(nmt.NmtField(mask, f))

crosses = []
for i in range(len(f2s)):
      for j in range(i,len(f2s)):
                cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
                crosses.append(cl_22)

crosses = np.array(crosses)
np.save('crosses_wmap9', crosses)

DIR = '/mn/stornext/d16/cmbco/bp/dwatts/WMAP/chains_CG_LFI_KaQVW_221124'
W1 = hp.read_map(f'{DIR}/tod_090-WMAP_W1_map_c0001_k000010.fits', field=(1,2))
W2 = hp.read_map(f'{DIR}/tod_090-WMAP_W2_map_c0001_k000010.fits', field=(1,2))
W3 = hp.read_map(f'{DIR}/tod_090-WMAP_W3_map_c0001_k000010.fits', field=(1,2))
W4 = hp.read_map(f'{DIR}/tod_090-WMAP_W4_map_c0001_k000010.fits', field=(1,2))

Ws = [W1, W2, W3, W4]
f2s = []
for f in Ws:
  f2s.append(nmt.NmtField(mask, f))

crosses = []
for i in range(len(f2s)):
      for j in range(i,len(f2s)):
                cl_22 = nmt.compute_full_master(f2s[i], f2s[j], b)
                crosses.append(cl_22)

crosses = np.array(crosses)
np.save('crosses_CG', crosses)

n1 = W1 - (W2 + W3 + W4)/3
n2 = W2 - (W1 + W3 + W4)/3
n3 = W3 - (W1 + W2 + W4)/3
n4 = W4 - (W1 + W2 + W3)/3
n5 = (W1 + W2)/2 - (W3 + W4)/2
n6 = (W1 + W3)/2 - (W2 + W4)/2

nulls = [n1, n2, n3, n4, n5, n6]
Cls = []
for null in nulls:
     f2 = nmt.NmtField(mask, null)
     cl_22 = nmt.compute_full_master(f2, f2, b)
     Cls.append(cl_22)
Cls = np.array(Cls)
np.save('null_spectra_CG', Cls)
'''
