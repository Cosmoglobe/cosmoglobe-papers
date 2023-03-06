import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cosmoglobe as cg
import astropy.units as u

from setup_matplotlib import *


from glob import glob


width = cm2inch(17)
height = 2**0.5*width

# Import the NaMaster python wrapper
import pymaster as nmt

#  Simple example showcasing the use of NaMaster to compute the pseudo-Cl
#  estimator of the angular cross-power spectrum of a spin-0 field and a
#  spin-2 field

# HEALPix resolution parameter used here
nside = 512

# Read mask and apodize it on a scale of ~1deg
mask = hp.read_map(
    "/mn/stornext/d5/data/duncanwa/WMAP/data/wmap_kq75_TQU_mask_r9.fits"
)


# bin4 = nmt.NmtBin.from_edges(l_ini, l_end)
logbins = np.geomspace(10, 3 * 512 - 1, 50).astype("int")
l_ini = np.concatenate((np.arange(2, 10), logbins[:-1]))
l_end = np.concatenate((np.arange(2, 10) + 1, logbins[1:]))
b = nmt.NmtBin.from_nside_linear(nside, 1)
#b = nmt.NmtBin.from_edges(l_ini, l_end)
ell_eff = b.get_effective_ells()

DIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"

wmap_maps = [
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_K1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Ka1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_Q2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_V2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W1_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W2_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W3_v5.fits",
    "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_W4_v5.fits",
]
CG_DIR = "/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203"
CG_DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

cg_maps = [
    f"{CG_DIR}/CG_023-WMAP_K_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_030-WMAP_Ka_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_040-WMAP_Q1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_040-WMAP_Q2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_060-WMAP_V1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_060-WMAP_V2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W1_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W2_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W3_IQU_n0512_v1.fits",
    f"{CG_DIR}/CG_090-WMAP_W4_IQU_n0512_v1.fits",
]

bands = ["K", "Ka", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]

x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip_W = -0.233 * x - 2.222 * y + 2.504 * z


fig1, axes1 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig2, axes2 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig3, axes3 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)

fig_large, axes_large = plt.subplots(sharex=True, nrows=10,
    ncols=1, figsize=(width/3, height))
plt.subplots_adjust(hspace=0)

axs1 = axes1.flatten()
axs2 = axes2.flatten()
axs3 = axes3.flatten()
axs_large = axes_large.flatten()

n = 0
for i in range(len(cg_maps)):
    print(i)
    try:
      Clhat_W = np.loadtxt(f'fullspec_WM_TT_{bands[i]}.txt')
      Clhat_C = np.loadtxt(f'fullspec_CG_TT_{bands[i]}.txt')
    except IOError:
      m_WMAP = hp.read_map(wmap_maps[i]) * 1e3
      m_CG = hp.read_map(cg_maps[i]) * 1e3
      m_CG -= dip_W * 1e3
      mono_CG = hp.fit_monopole(m_CG, gal_cut=50)
      mono_WM = hp.fit_monopole(m_WMAP, gal_cut=50)
      m_WMAP = m_WMAP - mono_WM + mono_CG
      f_WMAP = nmt.NmtField(mask, [m_WMAP])
      f_CG = nmt.NmtField(mask, [m_CG])
      f_diff = nmt.NmtField(mask, [m_CG - m_WMAP])
      Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)[0]
      print("ps1")
      Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)[0]
      print("ps2")
      # Clhat_diff = nmt.compute_full_master(f_diff, f_diff, b)[0]
      ell = np.arange(len(Clhat_C))
      print("got the power spectra")
      np.savetxt(f'fullspec_WM_TT_{bands[i]}.txt', Clhat_W)
      np.savetxt(f'fullspec_CG_TT_{bands[i]}.txt', Clhat_C)

    (l1,) = axs1[n].loglog(ell_eff, Clhat_W, label="WMAP")
    (l2,) = axs1[n].loglog(ell_eff, Clhat_C, label="Cosmoglobe")
    axs1[n].text(0.75, 0.75, r"\textit{" + bands[i] + "}", transform=axs1[n].transAxes)

    inds = ell_eff > 220
    axs2[n].plot(ell_eff[inds], Clhat_W[inds], label="WMAP")
    axs2[n].plot(ell_eff[inds], Clhat_C[inds], label="Cosmoglobe")
    axs2[n].text(0.75, 0.75, r"\textit{" + bands[i] + "}", transform=axs2[n].transAxes)
    axs2[n].set_ylim([0.002, 0.1])

    axs3[n].semilogx(ell_eff, Clhat_W / Clhat_C)
    axs3[n].text(0.25, 0.75, r"\textit{" + bands[i] + "}", transform=axs3[n].transAxes)

    axs_large[i].semilogx(ell_eff, Clhat_W / Clhat_C, 'k')
    axs_large[i].set_ylim([0.85, 1.15])
    #if i == 5:
    #    n += 3
    if i == 5:
        n += 3
    else:
        n += 1

axs1[6].axis("off")
axs1[7].axis("off")
axs2[6].axis("off")
axs2[7].axis("off")
axs3[6].axis("off")
axs3[7].axis("off")

axs1[6].legend(handles=[l1, l2], labels=[r"\textit{WMAP}", r"\textsc{Cosmoglobe}"])
axs2[6].legend(handles=[l1, l2], labels=[r"\textit{WMAP}", r"\textsc{Cosmoglobe}"])

fig1.supxlabel(r"$\ell$")
fig2.supxlabel(r"$\ell$")
fig3.supxlabel(r"$\ell$")

fig1.supylabel(r"$C_\ell^\mathrm{TT}\ [\mathrm{\mu K}^2]$")
fig2.supylabel(r"$C_\ell^\mathrm{TT}\ [\mathrm{\mu K}^2]$")
fig3.supylabel(r"$C_\ell^\mathit{WMAP}/C_\ell^\mathrm{Cosmoglobe}$")

plt.figure(fig1.number)
plt.savefig(f"TT_spectra.pdf", bbox_inches="tight")
plt.figure(fig2.number)
plt.savefig(f"TT_spectra_zoom.pdf", bbox_inches="tight")
plt.figure(fig3.number)
plt.savefig(f"TT_spectra_ratio.pdf", bbox_inches="tight")
print("here we are")
#plt.close("all")

axs_large[0].set_title(r'$C_\ell^{\mathrm{TT},\mathrm{Cosmoglobe}}/C_\ell^{\mathrm{TT},\mathit{WMAP}}$')
axs_large[9].set_xlabel(r'$\ell$')
plt.figure(fig_large.number)
plt.savefig('TT_ratio.pdf', bbox_inches='tight')

fig1, axes1 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)
fig2, axes2 = plt.subplots(sharex=True, sharey=True, nrows=3, ncols=4)
plt.subplots_adjust(wspace=0, hspace=0)


fig_large, axes_large = plt.subplots(sharex=True, sharey='row', nrows=10,
    ncols=2, figsize=(width/3*2, height))
plt.subplots_adjust(hspace=0, wspace=0)

axs1 = axes1.flatten()
axs2 = axes2.flatten()


n = 0
for i in range(len(cg_maps)):
    try:
      Clhat_W = np.loadtxt(f'fullspec_WM_EB_{bands[i]}.txt')
      Clhat_C = np.loadtxt(f'fullspec_CG_EB_{bands[i]}.txt')
    except IOError:
      m_WMAP = hp.read_map(wmap_maps[i], field=(1, 2)) * 1e3
      m_CG = hp.read_map(cg_maps[i], field=(1, 2)) * 1e3
      f_WMAP = nmt.NmtField(mask, m_WMAP)
      f_CG = nmt.NmtField(mask, m_CG)
      Clhat_W = nmt.compute_full_master(f_WMAP, f_WMAP, b)
      Clhat_C = nmt.compute_full_master(f_CG, f_CG, b)
      np.savetxt(f'fullspec_WM_EB_{bands[i]}.txt', Clhat_W)
      np.savetxt(f'fullspec_CG_EB_{bands[i]}.txt', Clhat_C)

    (l1,) = axs1[n].loglog(ell_eff, Clhat_W[0], label="WMAP")
    (l2,) = axs1[n].loglog(ell_eff, Clhat_C[0], label="Cosmoglobe")
    axs2[n].loglog(ell_eff, Clhat_W[3], label="WMAP")
    axs2[n].loglog(ell_eff, Clhat_C[3], label="Cosmoglobe")

    axs1[n].set_ylim([2e-3, 1e3])
    plt.figure(fig1.number)
    plt.savefig(f"EE_spectra.pdf", bbox_inches="tight")
    plt.figure(fig2.number)
    axs2[n].set_ylim([2e-3, 1e3])
    plt.savefig(f"BB_spectra.pdf", bbox_inches="tight")



    axes_large[i,0].loglog(ell_eff, Clhat_W[0])
    axes_large[i,0].loglog(ell_eff, Clhat_C[0])
    axes_large[i,1].loglog(ell_eff, Clhat_W[3])
    axes_large[i,1].loglog(ell_eff, Clhat_C[3])


    #axes_large[i,1].sharey(axes_large[i,2])

    axs1[n].text(0.75, 0.75, r"\textit{" + bands[i] + "}", transform=axs1[n].transAxes)
    axs2[n].text(0.75, 0.75, r"\textit{" + bands[i] + "}", transform=axs2[n].transAxes)
    axes_large[i,1].text(0.75, 0.75, r"\textit{" + bands[i] + "}", 
        transform=axes_large[i,1].transAxes)
    #if i == 5:
    #    n += 3
    if i == 3:
        n += 3
    else:
        n += 1

axs1[2].axis("off")
axs1[3].axis("off")
axs2[2].axis("off")
axs2[3].axis("off")

axs1[2].legend(handles=[l1, l2], labels=[r"\textit{WMAP}", r"\textsc{Cosmoglobe}"])
axs2[2].legend(handles=[l1, l2], labels=[r"\textit{WMAP}", r"\textsc{Cosmoglobe}"])
axes_large[0,0].legend(handles=[l1, l2], labels=[r"\textit{WMAP}", r"\textsc{Cosmoglobe}"], 
    fontsize=8)

axes_large[0,0].set_title(r'$C_\ell^\mathrm{EE}\ [\mathrm{\mu K}]$')
axes_large[0,1].set_title(r'$C_\ell^\mathrm{BB}\ [\mathrm{\mu K}]$')

axes_large[9,0].set_xlabel(r'$\ell$')
axes_large[9,1].set_xlabel(r'$\ell$')

fig1.supxlabel(r"$\ell$")
fig2.supxlabel(r"$\ell$")

fig1.supylabel(r"$C_\ell^\mathrm{EE}\ [\mathrm{\mu K}^2]$")
fig2.supylabel(r"$C_\ell^\mathrm{BB}\ [\mathrm{\mu K}^2]$")

plt.figure(fig1.number)
plt.savefig(f"EE_spectra.pdf", bbox_inches="tight")
plt.figure(fig2.number)
plt.savefig(f"BB_spectra.pdf", bbox_inches="tight")

plt.figure(fig_large)
plt.savefig('EE_BB_spec.pdf', bbox_inches='tight')
