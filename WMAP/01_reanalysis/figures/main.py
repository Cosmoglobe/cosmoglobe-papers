import cosmoglobe
import astropy.units as u
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

DIR = "/mn/stornext/d16/www_cmb/dwatts/v0"
WDIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"
BPDIR = "/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2"
width = 16
xsize = 800

bands = [
    "023-WMAP_K",
    "030-WMAP_Ka",
    "040-WMAP_Q1",
    "040-WMAP_Q2",
    "060-WMAP_V1",
    "060-WMAP_V2",
    "090-WMAP_W1",
    "090-WMAP_W2",
    "090-WMAP_W3",
    "090-WMAP_W4",
]

v = "version1"
v = "v0"


# Plot of all frequency bands

cosmoglobe.plot(
    f"{DIR}/BP_023-WMAP_K_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 1),
    sig=0,
    llabel="K",
)
cosmoglobe.plot(
    f"{DIR}/BP_023-WMAP_K_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 2),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    f"{DIR}/BP_023-WMAP_K_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 3),
    sig=2,
    fwhm=1 * u.deg,
)

cosmoglobe.plot(
    f"{DIR}/BP_030_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 4),
    sig=0,
    llabel="030",
)
cosmoglobe.plot(
    f"{DIR}/BP_030_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 5),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    f"{DIR}/BP_030_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 6),
    sig=2,
    fwhm=1 * u.deg,
)

cosmoglobe.plot(
    f"{DIR}/BP_030-WMAP_Ka_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 7),
    sig=0,
    llabel="\mathit{Ka}",
)
cosmoglobe.plot(
    f"{DIR}/BP_030-WMAP_Ka_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 8),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    f"{DIR}/BP_030-WMAP_Ka_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 9),
    sig=2,
    fwhm=1 * u.deg,
)

Q1 = hp.read_map(f"{DIR}/BP_040-WMAP_Q1_IQU_n0512_{v}.fits", field=(0, 1, 2))
Q2 = hp.read_map(f"{DIR}/BP_040-WMAP_Q2_IQU_n0512_{v}.fits", field=(0, 1, 2))
sQ1 = hp.read_map(f"{DIR}/BP_040-WMAP_Q1_IQU_n0512_{v}.fits", field=(6, 7, 8))
sQ2 = hp.read_map(f"{DIR}/BP_040-WMAP_Q2_IQU_n0512_{v}.fits", field=(6, 7, 8))
Q = (Q1 / sQ1**2 + Q2 / sQ2**2) / (1 / sQ1**2 + 1 / sQ2**2)
cosmoglobe.plot(
    Q,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 10),
    sig=0,
    llabel="Q",
)
cosmoglobe.plot(
    Q,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 11),
    sig=1,
    fwhm=2 * u.deg,
)
cosmoglobe.plot(
    Q,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 12),
    sig=2,
    fwhm=2 * u.deg,
)

cosmoglobe.plot(
    f"{DIR}/BP_044_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 13),
    sig=0,
    llabel="044",
)
cosmoglobe.plot(
    f"{DIR}/BP_044_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 14),
    sig=1,
    fwhm=2 * u.deg,
)
cosmoglobe.plot(
    f"{DIR}/BP_044_IQU_n0512_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 15),
    sig=2,
    fwhm=2 * u.deg,
)

V1 = hp.read_map(f"{DIR}/BP_060-WMAP_V1_IQU_n0512_{v}.fits", field=(0, 1, 2))
V2 = hp.read_map(f"{DIR}/BP_060-WMAP_V2_IQU_n0512_{v}.fits", field=(0, 1, 2))
sV1 = hp.read_map(f"{DIR}/BP_060-WMAP_V1_IQU_n0512_{v}.fits", field=(6, 7, 8))
sV2 = hp.read_map(f"{DIR}/BP_060-WMAP_V2_IQU_n0512_{v}.fits", field=(6, 7, 8))
V = (V1 / sV1**2 + V2 / sV2**2) / (1 / sV1**2 + 1 / sV2**2)
cosmoglobe.plot(
    V,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 16),
    sig=0,
    llabel="V",
)
cosmoglobe.plot(
    V,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 17),
    sig=1,
    fwhm=3 * u.deg,
)
cosmoglobe.plot(
    V,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 18),
    sig=2,
    fwhm=3 * u.deg,
)

cosmoglobe.plot(
    f"{DIR}/BP_070_IQU_n1024_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 19),
    sig=0,
    llabel="070",
)
cosmoglobe.plot(
    f"{DIR}/BP_070_IQU_n1024_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 20),
    sig=1,
    fwhm=4 * u.deg,
)
cosmoglobe.plot(
    f"{DIR}/BP_070_IQU_n1024_{v}.fits",
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 21),
    sig=2,
    fwhm=4 * u.deg,
)

W1 = hp.read_map(f"{DIR}/BP_090-WMAP_W1_IQU_n0512_{v}.fits", field=(0, 1, 2))
W2 = hp.read_map(f"{DIR}/BP_090-WMAP_W2_IQU_n0512_{v}.fits", field=(0, 1, 2))
W3 = hp.read_map(f"{DIR}/BP_090-WMAP_W3_IQU_n0512_{v}.fits", field=(0, 1, 2))
W4 = hp.read_map(f"{DIR}/BP_090-WMAP_W4_IQU_n0512_{v}.fits", field=(0, 1, 2))
sW1 = hp.read_map(f"{DIR}/BP_090-WMAP_W1_IQU_n0512_{v}.fits", field=(6, 7, 8))
sW2 = hp.read_map(f"{DIR}/BP_090-WMAP_W2_IQU_n0512_{v}.fits", field=(6, 7, 8))
sW3 = hp.read_map(f"{DIR}/BP_090-WMAP_W3_IQU_n0512_{v}.fits", field=(6, 7, 8))
sW4 = hp.read_map(f"{DIR}/BP_090-WMAP_W4_IQU_n0512_{v}.fits", field=(6, 7, 8))
W = (W1 / sW1**2 + W2 / sW2**2 + W3 / sW3**2 + W4 / sW4**2) / (
    1 / sW1**2 + 1 / sW2**2 + 1 / sW3**2 + 1 / sW4**2
)
cosmoglobe.plot(
    W,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 22),
    sig=0,
    llabel="W",
    unit=r"\mathrm{mK}",
    extend="both",
)
cosmoglobe.plot(
    W,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 23),
    sig=1,
    fwhm=5 * u.deg,
    unit=r"\mathrm{mK}",
    extend="both",
)
cosmoglobe.plot(
    W,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 24),
    sig=2,
    fwhm=5 * u.deg,
    unit=r"\mathrm{mK}",
    extend="both",
)

plt.subplots_adjust(wspace=0.04, hspace=-0.5)
plt.savefig("megaplot.pdf", bbox_inches="tight")
plt.close()


# Plot of all official maps
d_x = -0.233
d_y = -2.222
d_z = 2.504
x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip = d_x * x + d_y * y + d_z * z
K_w9 = hp.read_map(f"{WDIR}/wmap_band_iqusmap_r9_9yr_K_v5.fits", field=(0, 1, 2))
Ka_w9 = hp.read_map(f"{WDIR}/wmap_band_iqusmap_r9_9yr_Ka_v5.fits", field=(0, 1, 2))
Q_w9 = hp.read_map(f"{WDIR}/wmap_band_iqusmap_r9_9yr_Q_v5.fits", field=(0, 1, 2))
V_w9 = hp.read_map(f"{WDIR}/wmap_band_iqusmap_r9_9yr_V_v5.fits", field=(0, 1, 2))
W_w9 = hp.read_map(f"{WDIR}/wmap_band_iqusmap_r9_9yr_W_v5.fits", field=(0, 1, 2))

BP_030 = hp.read_map(f"{BPDIR}/BP_030_IQU_n0512_v2.fits", field=(0, 1, 2))
BP_044 = hp.read_map(f"{BPDIR}/BP_044_IQU_n0512_v2.fits", field=(0, 1, 2))
BP_070 = hp.read_map(f"{BPDIR}/BP_070_IQU_n1024_v2.fits", field=(0, 1, 2))

K_w9[0] += dip
Ka_w9[0] += dip
Q_w9[0] += dip
V_w9[0] += dip
W_w9[0] += dip

cosmoglobe.plot(
    K_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 1),
    sig=0,
    llabel="K",
)
cosmoglobe.plot(
    K_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 2),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    K_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 3),
    sig=2,
    fwhm=1 * u.deg,
)

cosmoglobe.plot(
    BP_030,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 4),
    sig=0,
    llabel="030",
)
cosmoglobe.plot(
    BP_030,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 5),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    BP_030,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 6),
    sig=2,
    fwhm=1 * u.deg,
)

cosmoglobe.plot(
    Ka_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 7),
    sig=0,
    llabel="\mathit{Ka}",
)
cosmoglobe.plot(
    Ka_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 8),
    sig=1,
    fwhm=1 * u.deg,
)
cosmoglobe.plot(
    Ka_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 9),
    sig=2,
    fwhm=1 * u.deg,
)

cosmoglobe.plot(
    Q_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 10),
    sig=0,
    llabel="Q",
)
cosmoglobe.plot(
    Q_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 11),
    sig=1,
    fwhm=2 * u.deg,
)
cosmoglobe.plot(
    Q_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 12),
    sig=2,
    fwhm=2 * u.deg,
)

cosmoglobe.plot(
    BP_044,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 13),
    sig=0,
    llabel="044",
)
cosmoglobe.plot(
    BP_044,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 14),
    sig=1,
    fwhm=2 * u.deg,
)
cosmoglobe.plot(
    BP_044,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 15),
    sig=2,
    fwhm=2 * u.deg,
)

cosmoglobe.plot(
    V_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 16),
    sig=0,
    llabel="V",
)
cosmoglobe.plot(
    V_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 17),
    sig=1,
    fwhm=3 * u.deg,
)
cosmoglobe.plot(
    V_w9,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 18),
    sig=2,
    fwhm=3 * u.deg,
)

cosmoglobe.plot(
    BP_070,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-3.4e3,
    max=3.4e3,
    sub=(8, 3, 19),
    sig=0,
    llabel="070",
)
cosmoglobe.plot(
    BP_070,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 20),
    sig=1,
    fwhm=4 * u.deg,
)
cosmoglobe.plot(
    BP_070,
    width=width,
    xsize=xsize,
    cbar=False,
    min=-0.025e3,
    max=0.025e3,
    sub=(8, 3, 21),
    sig=2,
    fwhm=4 * u.deg,
)

cosmoglobe.plot(
    W_w9,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-3.4,
    max=3.4,
    sub=(8, 3, 22),
    sig=0,
    llabel="W",
    unit=r"\mathrm{mK}",
    extend="both",
)
cosmoglobe.plot(
    W_w9,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 23),
    sig=1,
    fwhm=5 * u.deg,
    unit=r"\mathrm{mK}",
    extend="both",
)
cosmoglobe.plot(
    W_w9,
    width=width,
    xsize=xsize,
    cbar=True,
    min=-0.025,
    max=0.025,
    sub=(8, 3, 24),
    sig=2,
    fwhm=5 * u.deg,
    unit=r"\mathrm{mK}",
    extend="both",
)

plt.subplots_adjust(wspace=0.04, hspace=-0.5)
plt.savefig("megaplot_official.pdf", bbox_inches="tight")
plt.close()


# Plot of all diffs

DDIR = f"{DIR}/diffs"
m = hp.read_map(f"{DDIR}/BP_023-WMAP_K_diff_wmap9_{v}.fits", field=(0, 1, 2))
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 1),
    sig=0,
    llabel="K",
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 2),
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 3),
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)

m = hp.read_map(f"{DDIR}/BP_030_diff_BP10_{v}.fits", field=(0, 1, 2))
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 4),
    sig=0,
    llabel="030",
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 5),
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 6),
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)

m = hp.read_map(f"{DDIR}/BP_030-WMAP_Ka_diff_wmap9_{v}.fits", field=(0, 1, 2))
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 7),
    sig=0,
    llabel=r"\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 8),
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 9),
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)

cosmoglobe.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 10),
    sig=0,
    llabel="Q",
    min=-10,
    max=10,
    fwhm=2 * u.deg,
    remove_mono="auto",
    cbar=False,
)
cosmoglobe.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 11),
    sig=1,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 12),
    sig=2,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)

m = hp.read_map(f"{DDIR}/BP_044_diff_BP10_{v}.fits", field=(0, 1, 2))
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 13),
    sig=0,
    llabel="044",
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 14),
    sig=1,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 15),
    sig=2,
    cbar=False,
    min=-10,
    max=10,
)

cosmoglobe.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 16),
    sig=0,
    llabel="V",
    min=-10,
    max=10,
    remove_mono="auto",
    fwhm=2 * u.deg,
    cbar=False,
)
cosmoglobe.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 17),
    sig=1,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 18),
    sig=2,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)

cosmoglobe.plot(
    f"{DDIR}/BP_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 19),
    sig=0,
    llabel="070",
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    f"{DDIR}/BP_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 20),
    sig=1,
    min=-10,
    max=10,
    cbar=False,
)
cosmoglobe.plot(
    f"{DDIR}/BP_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 21),
    sig=2,
    min=-10,
    max=10,
    cbar=False,
)

cosmoglobe.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 22),
    sig=0,
    llabel="W",
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    fwhm=2 * u.deg,
    remove_mono="auto",
    cbar=True,
    ticks=[-10, -7.5, -5, -2.5, -1, 0, 1, 2.5, 5, 7.5, 10],
)
cosmoglobe.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 23),
    sig=1,
    fwhm=5 * u.deg,
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    cbar=True,
    ticks=[-10, -7.5, -5, -2.5, -1, 0, 1, 2.5, 5, 7.5, 10],
)
cosmoglobe.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sub=(8, 3, 24),
    sig=2,
    fwhm=5 * u.deg,
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    cbar=True,
    ticks=[-10, -7.5, -5, -2.5, -1, 0, 1, 2.5, 5, 7.5, 10],
)

plt.subplots_adjust(wspace=0.04, hspace=-0.5)
# plt.tight_layout()
plt.savefig("megadiff.pdf", bbox_inches="tight")
