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
fname_70 = f"{DIR}/CG_070_IQU_n1024_v1.fits"
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

d0 = -hp.ud_grade(hp.smoothing(0.495 * K - d_30, fwhm=fwhm) * 1e3, 256)
d1 = -hp.ud_grade(hp.smoothing(0.63 * d_30 - Ka, fwhm=fwhm) * 1e3, 256)
d2 = hp.ud_grade(hp.smoothing(Q1 - Q2, fwhm=fwhm) * 1e3, 256)
d3 = hp.ud_grade(hp.smoothing(V1 - V2, fwhm=fwhm) * 1e3, 256)
d4 = hp.ud_grade(hp.smoothing(((W1 - W2) - (W3 - W4)) / 4, fwhm=fwhm) * 1e3, 256)
d5 = -hp.ud_grade(hp.smoothing(alpha(44.1, 40.6, -3.1)*Q - d_44, fwhm=fwhm), 512)*1e3
d6 = hp.ud_grade(hp.smoothing( alpha(70, 60, -3.1)*V - d_70, fwhm=fwhm), 256)*1e3
d7 = hp.ud_grade(hp.smoothing( alpha(70, 90, -3.1)*W - d_70, fwhm=fwhm), 256)*1e3


plt.show()


cg.plot(
    d0,
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
    d0,
    sig=2,
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
    d1,
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
    d1,
    sig=2,
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
    d2,
    sig=1,
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    rlabel="Q, \mathrm{CG}",
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d2,
    sig=2,
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d3,
    sig=1,
    llabel=r"\Delta V",
    min=-10,
    max=10,
    cbar=False,
    rlabel="Q, \mathrm{CG}",
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d3,
    sig=2,
    rlabel="U, \mathrm{CG}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d4,
    sig=1,
    unit=r"\mathrm{\mu K}",
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
    d4,
    sig=2,
    rlabel="U, \mathrm{CG}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d5,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"44-Q",
    rlabel="Q, \mathrm{CG}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("44Q_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d5,
    sig=2,
    rlabel="U, \mathrm{CG}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("44Q_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"70-V",
    rlabel="Q, \mathrm{CG}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("V70_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=2,
    rlabel="U, \mathrm{CG}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("V70_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d7,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"W-70",
    rlabel="Q, \mathrm{CG}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W70_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=2,
    rlabel="U, \mathrm{CG}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W70_deltaU.pdf", bbox_inches="tight", dpi=300)
plt.close("all")


Clhat0 = hp.anafast(d0)
Clhat1 = hp.anafast(d1)
Clhat2 = hp.anafast(d2)
Clhat3 = hp.anafast(d3)
Clhat4 = hp.anafast(d4)
Clhat5 = hp.anafast(d5)
Clhat6 = hp.anafast(d6)
Clhat7 = hp.anafast(d7)

np.savetxt('Clhat_K30_CG.txt', Clhat0)
np.savetxt('Clhat_30Ka_CG.txt', Clhat1)
np.savetxt('Clhat_dQ_CG.txt', Clhat2)
np.savetxt('Clhat_dV_CG.txt', Clhat3)
np.savetxt('Clhat_dW_CG.txt', Clhat4)
np.savetxt('Clhat_44Q_CG.txt', Clhat5)
np.savetxt('Clhat_70V_CG.txt', Clhat6)
np.savetxt('Clhat_70W_CG.txt', Clhat7)


DIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"
fnames = glob(f"{DIR}/wmap_iqusmap_r9_9yr*.fits")
fnames.sort()

BP_DIR = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'

fname_30 = f"{BP_DIR}/BP_030_IQU_n0512_v2.fits"
d_30 = hp.read_map(fname_30, field=(0, 1, 2)) * 1e-3
fname_44 = f"{BP_DIR}/BP_044_IQU_n0512_v2.fits"
d_44 = hp.read_map(fname_44, field=(0, 1, 2)) * 1e-3
fname_70 = f"{BP_DIR}/BP_070_IQU_n1024_v2.fits"
d_70 = hp.read_map(fname_70, field=(0, 1, 2)) * 1e-3
d_70 = hp.ud_grade(d_70, 512)

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

d0 = -hp.ud_grade(hp.smoothing(0.495 * K - d_30, fwhm=fwhm) * 1e3, 256)
d1 = -hp.ud_grade(hp.smoothing(0.63 * d_30 - Ka, fwhm=fwhm) * 1e3, 256)
d2 = hp.ud_grade(hp.smoothing(Q1 - Q2, fwhm=fwhm) * 1e3, 256)
d3 = hp.ud_grade(hp.smoothing(V1 - V2, fwhm=fwhm) * 1e3, 256)
d4 = hp.ud_grade(hp.smoothing(((W1 - W2) - (W3 - W4)) / 4, fwhm=fwhm) * 1e3, 256)
d5 = -hp.ud_grade(hp.smoothing(alpha(44.1, 40.6, -3.1)*Q - d_44, fwhm=fwhm), 256)*1e3
d6 = hp.ud_grade(hp.smoothing( alpha(70, 60, -3.1)*V - d_70, fwhm=fwhm), 256)*1e3
d7 = hp.ud_grade(hp.smoothing( alpha(70, 90, -3.1)*W - d_70, fwhm=fwhm), 256)*1e3


cg.plot(
    d0,
    sig=1,
    llabel=r"30-\mathit{K}",
    rlabel="Q, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d0,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("K30_W_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    d1,
    sig=1,
    llabel=r"\mathit{Ka}-30",
    rlabel="Q, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d1,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("30Ka_W_deltaU.pdf", bbox_inches="tight", dpi=300)

cg.plot(
    d1,
    sig=1,
    llabel=r"\mathit K-\mathit{Ka}",
    rlabel="Q, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d1,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("KKa_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d2,
    sig=1,
    llabel=r"\Delta Q",
    min=-10,
    max=10,
    cbar=False,
    rlabel="Q, \mathit{WMAP}",
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d2,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("Q_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d3,
    sig=1,
    llabel=r"\Delta V",
    min=-10,
    max=10,
    cbar=False,
    rlabel="Q, \mathit{WMAP}",
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d3,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    cbar=False,
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    fontsize=fontsize,
)
plt.savefig("V_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d4,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"\Delta W",
    rlabel="Q, \mathit{WMAP}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d4,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d5,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"44-Q",
    rlabel="Q, \mathit{WMAP}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("44Q_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d5,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("44Q_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"70-V",
    rlabel="Q, \mathit{WMAP}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("V70_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=2,
    rlabel="U, \mathrm{CG}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("V70_W_deltaU.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d7,
    sig=1,
    unit=r"\mathrm{\mu K}",
    llabel=r"W-70",
    rlabel="Q, \mathit{WMAP}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W70_W_deltaQ.pdf", bbox_inches="tight", dpi=300)
cg.plot(
    d6,
    sig=2,
    rlabel="U, \mathit{WMAP}",
    unit=r"\mathrm{\mu K}",
    min=-10,
    max=10,
    xsize=xsize,
    width=width,
    extend="both",
    fontsize=fontsize,
    cbar=False,
)
plt.savefig("W70_W_deltaU.pdf", bbox_inches="tight", dpi=300)
plt.close("all")

Clhat0 = hp.anafast(d0)
Clhat1 = hp.anafast(d1)
Clhat2 = hp.anafast(d2)
Clhat3 = hp.anafast(d3)
Clhat4 = hp.anafast(d4)
Clhat5 = hp.anafast(d5)
Clhat6 = hp.anafast(d6)
Clhat7 = hp.anafast(d7)

np.savetxt('Clhat_K30_WM.txt', Clhat0)
np.savetxt('Clhat_30Ka_WM.txt', Clhat1)
np.savetxt('Clhat_dQ_WM.txt', Clhat2)
np.savetxt('Clhat_dV_WM.txt', Clhat3)
np.savetxt('Clhat_dW_WM.txt', Clhat4)
np.savetxt('Clhat_44Q_WM.txt', Clhat5)
np.savetxt('Clhat_70V_WM.txt', Clhat6)
np.savetxt('Clhat_70W_WM.txt', Clhat7)
