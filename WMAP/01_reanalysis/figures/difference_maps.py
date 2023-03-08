import numpy as np
import healpy as hp
import cosmoglobe as cg
import matplotlib.pyplot as plt
import astropy.units as u


fwhm = 2  # degrees
fwhm_P = 32 / 60  # degrees
fwhm_W = 0.88  # degrees
width = 4

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203'


try:
    P_18 = hp.read_map("P_18.fits", field=(0, 1, 2))
    P_NPIPE = hp.read_map("P_NPIPE.fits", field=(0, 1, 2))
    P_BP = hp.read_map("P_BP.fits", field=(0, 1, 2))
except FileNotFoundError:

    P_18 = (
        hp.read_map(
            "/mn/stornext/d16/cmbco/ola/planck_products/2018/LFI_SkyMap_030-BPassCorrected-field-IQU_1024_R3.00_full.fits",
            field=(0, 1, 2),
        )
        * 1e6
    )
    P_NPIPE = (
        hp.read_map(
            "/mn/stornext/d16/cmbco/ola/npipe/freqmaps/npipe6v20_030_map_K.fits",
            field=(0, 1, 2),
        )
        * 1e6
    )
    P_BP = hp.read_map(
        "/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2/BP_030_IQU_n0512_v2.fits",
        field=(0, 1, 2),
    )

    P_18 = hp.smoothing(P_18, fwhm=(fwhm**2 - fwhm_P**2) ** 0.5 * np.pi / 180)
    P_NPIPE = hp.smoothing(P_NPIPE, fwhm=(fwhm**2 - fwhm_P**2) ** 0.5 * np.pi / 180)
    P_BP = hp.smoothing(P_BP, fwhm=(fwhm**2 - fwhm_P**2) ** 0.5 * np.pi / 180)
    hp.write_map("P_18.fits", P_18)
    hp.write_map("P_NPIPE.fits", P_NPIPE)
    hp.write_map("P_BP.fits", P_BP)

P_18 = hp.ud_grade(P_18, 512)
P_NPIPE = hp.ud_grade(P_NPIPE, 512)
# cg.plot(P_18, sig=1, min=-20, max=20)
# cg.plot(P_NPIPE, sig=1, min=-20, max=20)
# cg.plot(P_BP, sig=1, min=-20, max=20)

W_DR5 = (
    hp.read_map(
        "/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_band_iqusmap_r9_9yr_K_v5.fits",
        field=(0, 1, 2),
    )
    * 1e3
)
W_CG = (
    hp.read_map(
        f"{DIR}/tod_023-WMAP_K_map_c0001_k000009.fits",
        field=(0, 1, 2),
    )
    * 1e3
)
P_CG = hp.read_map(
    f"{DIR}/tod_030_map_c0001_k000009.fits",
    field=(0, 1, 2),
)

W_DR5 = hp.smoothing(W_DR5, fwhm=(fwhm**2 - fwhm_W**2) ** 0.5 * np.pi / 180)
W_CG = hp.smoothing(W_CG, fwhm=(fwhm**2 - fwhm_W**2) ** 0.5 * np.pi / 180)
P_CG = hp.smoothing(P_CG, fwhm=(fwhm**2 - fwhm_P**2) ** 0.5 * np.pi / 180)


# cg.plot(W_DR5, sig=1, min=-20, max=20)
# cg.plot(W_CG, sig=1, min=-20, max=20)

cg.plot(
    P_18 - 0.495 * W_DR5,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
    rlabel='Q',
    width=width,
)
plt.savefig(f'diff_18_DR5_Q.pdf', bbox_inches='tight')
cg.plot(
    P_18 - 0.495 * W_DR5,
    sig=2,
    min=-10,
    max=10,
    cbar=False,
    rlabel=r"U",
    width=width,
)
plt.text(3.3,  0, r"\textit{Planck} 2018", ha='left', va='center', fontsize=16)
plt.savefig(f'diff_18_DR5_U.pdf', bbox_inches='tight')
cg.plot(
    P_NPIPE - 0.495 * W_DR5,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
    rlabel='Q',
    width=width,
)
plt.savefig(f'diff_NPIPE_DR5_Q.pdf', bbox_inches='tight')
cg.plot(
    P_NPIPE - 0.495 * W_DR5,
    sig=2,
    min=-10,
    max=10,
    cbar=False,
    rlabel="U",
    width=width,
)
plt.text(3.3,  0, r"\textit{Planck} PR4", ha='left', va='center', fontsize=16)
plt.savefig(f'diff_NPIPE_DR5_U.pdf', bbox_inches='tight')
cg.plot(
    P_BP - 0.495 * W_DR5,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
  rlabel='Q',
    width=width,
)
plt.savefig(f'diff_BP_DR5_Q.pdf', bbox_inches='tight')
cg.plot(
    P_BP - 0.495 * W_DR5,
    sig=2,
    min=-10,
    max=10,
    cbar=False,
  rlabel="U",
    width=width,
)
plt.text(3.3,  0, r"\textsc{BeyondPlanck}", ha='left', va='center', fontsize=16)
plt.savefig(f'diff_BP_DR5_U.pdf', bbox_inches='tight')
cg.plot(
    P_CG - 0.495 * W_CG,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
  rlabel='Q',
    width=width,
)
plt.savefig(f'diff_CG_Q.pdf', bbox_inches='tight')
cg.plot(
    P_CG - 0.495 * W_CG,
    sig=2,
    min=-10,
    max=10,
    cbar=False,
  rlabel="U",
    width=width,
)
plt.text(3.3,  0, r"\textsc{Cosmoglobe}", ha='left', va='center', fontsize=16)
plt.savefig(f'diff_CG_U.pdf', bbox_inches='tight')
plt.show()

cg.standalone_colorbar("planck", ticks=[-10,0,10], extend='both',
    unit=r"$\mathrm{\mu K_{CMB}}$", fontsize=18, width=6)

plt.savefig('cbar_10uK.pdf', bbox_inches='tight')
