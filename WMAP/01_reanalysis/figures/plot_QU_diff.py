import healpy as hp
import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

def alpha(nu, nu0, beta):
    S = (nu/nu0)**beta
    return S*g(nu/56.78)/g(nu0/56.78)
def g(x):
    return np.expm1(x)**2/(x**2*np.exp(x))

d_100 = hp.read_map("/mn/stornext/d16/cmbco/ola/planck_products/dr4/HFI_SkyMap_100-BPassCorrected-field-IQU_2048_R4.00_full.fits",field=(0,1,2)) * 1e3
d_100 = hp.ud_grade(d_100,512)

DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

width = 2
xsize = 1200
fontsize = {
    "llabel": 5,
    "rlabel": 5,
}

fnames = glob(f"{DIR}/CG_???-WMAP_*v1.fits")
fnames.sort()
fname_30 = f"{DIR}/CG_030_IQU_n0512_v1.fits"
d_30 = hp.read_map(fname_30, field=(0, 1, 2)) * 1e-3
fname_44 = f"{DIR}/CG_044_IQU_n0512_v1.fits"
d_44 = hp.read_map(fname_44, field=(0, 1, 2)) * 1e-3
#fname_70 = f"{DIR}/CG_070_IQU_n1024_v1.fits"
fname_70 = "/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2/BP_070_IQU_n1024_v2.fits"
d_70 = hp.read_map(fname_70, field=(0, 1, 2)) * 1e-3
d_70 = hp.ud_grade(d_70, 512)

fwhm = 5 * np.pi / 180

print(fnames)

K = hp.read_map(fnames[0], field=(0, 1, 2))
Ka = hp.read_map(fnames[1], field=(0, 1, 2))
Q1 = hp.read_map(fnames[2], field=(0, 1, 2))
Q2 = hp.read_map(fnames[3], field=(0, 1, 2))
V1 = hp.read_map(fnames[4], field=(0, 1, 2))
V2 = hp.read_map(fnames[5], field=(0, 1, 2))
W1 = hp.read_map(fnames[6], field=(0, 1, 2))
W2 = hp.read_map(fnames[7], field=(0, 1, 2))
W3 = hp.read_map(fnames[8], field=(0, 1, 2))
W4 = hp.read_map(fnames[9], field=(0, 1, 2))

Q = (Q1 + Q2)/2
V = (V1 + V2)/2
W = (W1 + W2 + W3 + W4)/4

q44_alpha = (44./40.6)**(-3.1)
K_Ka_alpha = (33./23.)**(-3.1)

K30  = -hp.ud_grade(hp.smoothing(0.495 * K - d_30, fwhm=fwhm) * 1e3, 256)
Ka30 = -hp.ud_grade(hp.smoothing(0.63 * d_30 - Ka, fwhm=fwhm) * 1e3, 256)
Q44  = -hp.ud_grade(hp.smoothing(q44_alpha * Q - d_44, fwhm=fwhm)*1e3, 256)
V70  = hp.ud_grade(hp.smoothing(V-d_70,fwhm=fwhm)*1e3, 256)
W70  = hp.ud_grade(hp.smoothing(W-d_70,fwhm=fwhm)*1e3, 256) 
W100 = hp.ud_grade(hp.smoothing(d_100-W,fwhm=fwhm)*1e3, 256)

deltaKKa = hp.ud_grade(hp.smoothing(K_Ka_alpha*K-Ka, fwhm=fwhm) * 1e3, 256)
deltaQ = hp.ud_grade(hp.smoothing(Q1 - Q2, fwhm=fwhm) * 1e3, 256)
deltaV = hp.ud_grade(hp.smoothing(V1 - V2, fwhm=fwhm) * 1e3, 256)
deltaW = hp.ud_grade(hp.smoothing(((W1 - W2) - (W3 - W4)) / 4, fwhm=fwhm) * 1e3, 256)

Clhat0 = hp.anafast(K30)
Clhat1 = hp.anafast(Ka30)
Clhat2 = hp.anafast(deltaQ)
Clhat3 = hp.anafast(deltaV)
Clhat4 = hp.anafast(deltaW)
Clhat5 = hp.anafast(Q44)
Clhat6 = hp.anafast(V70)
Clhat7 = hp.anafast(W100)

np.savetxt('../data/powspec_intdiff/Clhat_K30_CG.txt', Clhat0)
np.savetxt('../data/powspec_intdiff/Clhat_30Ka_CG.txt', Clhat1)
np.savetxt('../data/powspec_intdiff/Clhat_dQ_CG.txt', Clhat2)
np.savetxt('../data/powspec_intdiff/Clhat_dV_CG.txt', Clhat3)
np.savetxt('../data/powspec_intdiff/Clhat_dW_CG.txt', Clhat4)
np.savetxt('../data/powspec_intdiff/Clhat_44Q_CG.txt', Clhat5)
np.savetxt('../data/powspec_intdiff/Clhat_70V_CG.txt', Clhat6)
np.savetxt('../data/powspec_intdiff/Clhat_70W_CG.txt', Clhat7)

plt.show()

################################################################
# LFI-WMAP Compare -- CG
################################################################

cg.plot(
    K30,
    sig=1,
    llabel=r"30-\mathit{K}",
    rlabel="Q, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    K30,
    sig=2,
    llabel=r"30-\mathit{K}",
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Ka30,
    sig=1,
    llabel=r"\mathit{Ka}-30",
    rlabel="Q, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Ka30,
    sig=2,
    llabel=r"\mathit{Ka}-30",
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Q44,
    sig=1,
    llabel=r"\mathit{Q}-44",
    rlabel="Q, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("44Q_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Q44,
    sig=2,
    llabel=r"\mathit{Q}-44",
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("44Q_deltaU.pdf", bbox_inches="tight", dpi=300
)

cg.plot(
    V70,
    sig=1,
    llabel=r"\mathit{V}-70",
    rlabel="Q, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("70V_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    V70,
    sig=2,
    llabel=r"\mathit{V}-70",
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("70V_deltaU.pdf", bbox_inches="tight", dpi=300
)

cg.plot(
    W100,
    sig=1,
    llabel=r"100-\mathit{W}",
    rlabel="Q, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("100W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    W100,
    sig=2,
    llabel=r"100-\mathit{W}",
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("100W_deltaU.pdf", bbox_inches="tight", dpi=300
)

################################################################
# Intra-channel diffs
################################################################
cg.plot(
    deltaKKa,
    sig=1,
    rlabel="Q, \mathrm{CG}",
    llabel=r"\mathit{K}-\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaKKa,
    sig=2,
    rlabel="U, \mathrm{CG}",
    llabel=r"\mathit{K}-\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaQ,
    sig=1,
    rlabel="Q, \mathrm{CG}",
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaQ,
    sig=2,
    rlabel="U, \mathrm{CG}",
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaV,
    sig=1,
    rlabel="Q, \mathrm{CG}",
    llabel=r"\Delta V",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaV,
    sig=2,
    rlabel="U, \mathrm{CG}",
    llabel=r"\Delta V",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaW,
    sig=1,
    llabel=r"\Delta W",
    rlabel="Q, \mathrm{CG}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaW,
    sig=2,
    llabel=r"\Delta W",
    rlabel="U, \mathrm{CG}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W_deltaU.pdf", bbox_inches="tight", dpi=300)

################################################################
# LFI-WMAP Compare -- WMAP
################################################################

DIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"
fnames = glob(f"{DIR}/wmap_iqusmap_r9_9yr*.fits")
fnames.sort()

BP_DIR = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'

fname_30 = f"{BP_DIR}/BP_030_IQU_n0512_v2.fits"
d_30_bp = hp.read_map(fname_30, field=(0, 1, 2)) * 1e-3
fname_44 = f"{BP_DIR}/BP_044_IQU_n0512_v2.fits"
d_44_bp = hp.read_map(fname_44, field=(0, 1, 2)) * 1e-3
fname_70 = f"{BP_DIR}/BP_070_IQU_n1024_v2.fits"
d_70_bp = hp.read_map(fname_70, field=(0, 1, 2)) * 1e-3
d_70_bp = hp.ud_grade(d_70_bp, 512)

K_wmap = hp.read_map(fnames[0], field=(0, 1, 2))
Ka_wmap = hp.read_map(fnames[1], field=(0, 1, 2))
Q1_wmap = hp.read_map(fnames[2], field=(0, 1, 2))
Q2_wmap = hp.read_map(fnames[3], field=(0, 1, 2))
V1_wmap = hp.read_map(fnames[4], field=(0, 1, 2))
V2_wmap = hp.read_map(fnames[5], field=(0, 1, 2))
W1_wmap = hp.read_map(fnames[6], field=(0, 1, 2))
W2_wmap = hp.read_map(fnames[7], field=(0, 1, 2))
W3_wmap = hp.read_map(fnames[8], field=(0, 1, 2))
W4_wmap = hp.read_map(fnames[9], field=(0, 1, 2))

Q_wmap = (Q1_wmap + Q2_wmap)/2
V_wmap = (V1_wmap + V2_wmap)/2
W_wmap = (W1_wmap + W2_wmap + W3_wmap + W4_wmap)/4

K30_wmap  = -hp.ud_grade(hp.smoothing(0.495 * K_wmap - d_30_bp, fwhm=fwhm) * 1e3, 256)
Ka30_wmap = -hp.ud_grade(hp.smoothing(0.63 * d_30_bp - Ka_wmap, fwhm=fwhm) * 1e3, 256)
Q44_wmap  = -hp.ud_grade(hp.smoothing(q44_alpha * Q_wmap - d_44_bp, fwhm=fwhm)*1e3, 256)
V70_wmap  = hp.ud_grade(hp.smoothing(V_wmap-d_70_bp,fwhm=fwhm)*1e3, 256)
W70_wmap  = hp.ud_grade(hp.smoothing(W_wmap-d_70_bp,fwhm=fwhm)*1e3, 256)
W100_wmap = hp.ud_grade(hp.smoothing(d_100-W_wmap,fwhm=fwhm)*1e3, 256)

deltaKKa_wmap = hp.ud_grade(hp.smoothing(K_Ka_alpha*K_wmap-Ka_wmap, fwhm=fwhm) * 1e3, 256)
deltaQ_wmap = hp.ud_grade(hp.smoothing(Q1_wmap - Q2_wmap, fwhm=fwhm) * 1e3, 256)
deltaV_wmap = hp.ud_grade(hp.smoothing(V1_wmap - V2_wmap, fwhm=fwhm) * 1e3, 256)
deltaW_wmap = hp.ud_grade(hp.smoothing(((W1_wmap - W2_wmap) - (W3_wmap - W4_wmap)) / 4, fwhm=fwhm) * 1e3, 256)

Clhat0 = hp.anafast(K30_wmap)
Clhat1 = hp.anafast(Ka30_wmap)
Clhat2 = hp.anafast(deltaQ_wmap)
Clhat3 = hp.anafast(deltaV_wmap)
Clhat4 = hp.anafast(deltaW_wmap)
Clhat5 = hp.anafast(Q44_wmap)
Clhat6 = hp.anafast(V70_wmap)
Clhat7 = hp.anafast(W100_wmap)

np.savetxt('../data/powspec_intdiff/Clhat_K30_WM.txt', Clhat0)
np.savetxt('../data/powspec_intdiff/Clhat_30Ka_WM.txt', Clhat1)
np.savetxt('../data/powspec_intdiff/Clhat_dQ_WM.txt', Clhat2)
np.savetxt('../data/powspec_intdiff/Clhat_dV_WM.txt', Clhat3)
np.savetxt('../data/powspec_intdiff/Clhat_dW_WM.txt', Clhat4)
np.savetxt('../data/powspec_intdiff/Clhat_44Q_WM.txt', Clhat5)
np.savetxt('../data/powspec_intdiff/Clhat_70V_WM.txt', Clhat6)
np.savetxt('../data/powspec_intdiff/Clhat_70W_WM.txt', Clhat7)


cg.plot(
    K30_wmap,
    sig=1,
    llabel=r"30-\mathit{K}",
    rlabel="Q, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    K30_wmap,
    sig=2,
    llabel=r"30-\mathit{K}",
    rlabel="U, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_W_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Ka30_wmap,
    sig=1,
    llabel=r"\mathit{Ka}-30",
    rlabel="Q, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Ka30_wmap,
    sig=2,
    llabel=r"\mathit{Ka}-30",
    rlabel="U, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_W_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Q44_wmap,
    sig=1,
    llabel=r"\mathit{Q}-44",
    rlabel="Q, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("44Q_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    Q44_wmap,
    sig=2,
    llabel=r"\mathit{Q}-44",
    rlabel="U, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("44Q_W_deltaU.pdf", bbox_inches="tight", dpi=300
)

cg.plot(
    V70_wmap,
    sig=1,
    llabel=r"\mathit{V}-70",
    rlabel="Q, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("70V_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    V70_wmap,
    sig=2,
    llabel=r"\mathit{V}-70",
    rlabel="U, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("70V_W_deltaU.pdf", bbox_inches="tight", dpi=300
)

cg.plot(
    W100_wmap,
    sig=1,
    llabel=r"100-\mathit{W}",
    rlabel="Q, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("100W_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    W100_wmap,
    sig=2,
    llabel=r"100-\mathit{W}",
    rlabel="U, \mathrm{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("100W_W_deltaU.pdf", bbox_inches="tight", dpi=300
)

################################################################
# Intra-channel diff maps
################################################################
cg.plot(
    deltaKKa_wmap,
    sig=1,
    rlabel="Q, \mathrm{WMAP}",
    llabel=r"\mathit{K}-\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaKKa_wmap,
    sig=2,
    rlabel="U, \mathrm{WMAP}",
    llabel=r"\mathit{K}-\mathit{Ka}",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_W_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaQ_wmap,
    sig=1,
    rlabel="Q, \mathrm{WMAP}",
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_W_deltaQ.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    deltaQ_wmap,
    sig=2,
    rlabel="U, \mathrm{WMAP}",
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaV_wmap,
    sig=1,
    rlabel="Q, \mathrm{WMAP}",
    llabel=r"\Delta V",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaV_wmap,
    sig=2,
    rlabel="U, \mathrm{WMAP}",
    llabel=r"\Delta V",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaW_wmap,
    sig=1,
    rlabel="Q, \mathrm{WMAP}",
    llabel=r"\Delta W",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
)
plt.savefig("W_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    deltaW_wmap,
    sig=2,
    rlabel="U, \mathrm{WMAP}",
    llabel=r"\Delta W",
    min=-10,
    max=10,
    cbar=False,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
)
plt.savefig("W_W_deltaU.pdf", bbox_inches="tight", dpi=300)

################################################################
################################################################
