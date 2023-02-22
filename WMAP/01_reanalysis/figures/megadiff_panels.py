import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"
WDIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"
BPDIR = "/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2"
width = 5
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
v = "v1"


Q1 = hp.read_map(f"{DIR}/CG_040-WMAP_Q1_IQU_n0512_{v}.fits", field=(0, 1, 2))
Q2 = hp.read_map(f"{DIR}/CG_040-WMAP_Q2_IQU_n0512_{v}.fits", field=(0, 1, 2))
sQ1 = hp.read_map(f"{DIR}/CG_040-WMAP_Q1_IQU_n0512_{v}.fits", field=(7,8,9))
sQ2 = hp.read_map(f"{DIR}/CG_040-WMAP_Q2_IQU_n0512_{v}.fits", field=(7,8,9))
Q = (Q1 / sQ1 + Q2 / sQ2) / (1 / sQ1 + 1 / sQ2)

V1 = hp.read_map(f"{DIR}/CG_060-WMAP_V1_IQU_n0512_{v}.fits", field=(0, 1, 2))
V2 = hp.read_map(f"{DIR}/CG_060-WMAP_V2_IQU_n0512_{v}.fits", field=(0, 1, 2))
sV1 = hp.read_map(f"{DIR}/CG_060-WMAP_V1_IQU_n0512_{v}.fits", field=(7,8,9))
sV2 = hp.read_map(f"{DIR}/CG_060-WMAP_V2_IQU_n0512_{v}.fits", field=(7,8,9))
V = (V1 / sV1 + V2 / sV2) / (1 / sV1 + 1 / sV2)

W1 = hp.read_map(f"{DIR}/CG_090-WMAP_W1_IQU_n0512_{v}.fits", field=(0, 1, 2))
W2 = hp.read_map(f"{DIR}/CG_090-WMAP_W2_IQU_n0512_{v}.fits", field=(0, 1, 2))
W3 = hp.read_map(f"{DIR}/CG_090-WMAP_W3_IQU_n0512_{v}.fits", field=(0, 1, 2))
W4 = hp.read_map(f"{DIR}/CG_090-WMAP_W4_IQU_n0512_{v}.fits", field=(0, 1, 2))
sW1 = hp.read_map(f"{DIR}/CG_090-WMAP_W1_IQU_n0512_{v}.fits", field=(7,8,9))
sW2 = hp.read_map(f"{DIR}/CG_090-WMAP_W2_IQU_n0512_{v}.fits", field=(7,8,9))
sW3 = hp.read_map(f"{DIR}/CG_090-WMAP_W3_IQU_n0512_{v}.fits", field=(7,8,9))
sW4 = hp.read_map(f"{DIR}/CG_090-WMAP_W4_IQU_n0512_{v}.fits", field=(7,8,9))
W = (W1 / sW1 + W2 / sW2 + W3 / sW3 + W4 / sW4) / (
    1 / sW1 + 1 / sW2 + 1 / sW3 + 1 / sW4
)

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


# Plot of all diffs

DDIR = f"{DIR}/diffs"
m = hp.read_map(f"{DDIR}/CG_023-WMAP_K_diff_wmap9_{v}.fits", field=(0, 1, 2))
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="K",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_K_I.pdf')
plt.close()
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_K_Q.pdf')
plt.close()
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_K_U.pdf')
plt.close()

m = hp.read_map(f"{DDIR}/CG_030_diff_BP10_{v}.fits", field=(0, 1, 2))
cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="030",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_030_I.pdf')
plt.close()

cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_030_Q.pdf')
plt.close()
cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_030_U.pdf')
plt.close()

m = hp.read_map(f"{DDIR}/CG_030-WMAP_Ka_diff_wmap9_{v}.fits", field=(0, 1, 2))
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel=r"\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_Ka_I.pdf')
plt.close()
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_Ka_Q.pdf')
plt.close()
cg.plot(
    m * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=1 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_Ka_U.pdf')
plt.close()

cg.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="Q",
    min=-10,
    max=10,
    fwhm=2 * u.deg,
    remove_mono="auto",
    cbar=False,
)
plt.savefig('megadiff_Q_I.pdf')
plt.close()
cg.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_Q_Q.pdf')
plt.close()
cg.plot(
    (Q - Q_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_Q_U.pdf')
plt.close()

m = hp.read_map(f"{DDIR}/CG_044_diff_BP10_{v}.fits", field=(0, 1, 2))
cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="044",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_044_I.pdf')
plt.close()
cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_044_Q.pdf')
plt.close()
cg.plot(
    m,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    cbar=False,
    min=-10,
    max=10,
)
plt.savefig('megadiff_044_U.pdf')
plt.close()

cg.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="V",
    min=-10,
    max=10,
    remove_mono="auto",
    fwhm=2 * u.deg,
    cbar=False,
)
plt.savefig('megadiff_V_I.pdf')
plt.close()
cg.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_V_Q.pdf')
plt.close()
cg.plot(
    (V - V_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=2 * u.deg,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_V_U.pdf')
plt.close()

cg.plot(
    f"{DDIR}/CG_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="070",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_070_I.pdf')
plt.close()
cg.plot(
    f"{DDIR}/CG_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_070_Q.pdf')
plt.close()
cg.plot(
    f"{DDIR}/CG_070_diff_BP10_{v}.fits",
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_070_U.pdf')
plt.close()

cg.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=0,
    llabel="W",
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    fwhm=2 * u.deg,
    remove_mono="auto",
    cbar=False,
)
plt.savefig('megadiff_W_I.pdf')
plt.close()
cg.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=1,
    fwhm=5 * u.deg,
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_W_Q.pdf')
plt.close()
cg.plot(
    (W - W_w9) * 1e3,
    norm="symlog2",
    width=width,
    xsize=xsize,
    sig=2,
    fwhm=5 * u.deg,
    unit=r"\mathrm{\mu K}",
    extend="both",
    min=-10,
    max=10,
    cbar=False,
)
plt.savefig('megadiff_W_U.pdf')
plt.close()

    
cg.standalone_colorbar("planck_log",  ticks=[-10, -7.5, -5, -2.5, -1, 0, 1, 2.5, 5, 7.5, 10],
    extend='both',
            unit=r"$\mathrm{\mu K}$", width=width*2)
plt.savefig('cbar_10uK_symlog2.pdf', bbox_inches='tight')

