import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

rng = np.random.default_rng()

width = 6
xsize = 1200

dpi = 150

DIR1 = "/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_a_230206"
DIR2 = "/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203"

RDIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

WDIR = "/mn/stornext/d16/cmbco/ola/wmap/freq_maps"


burn_in = 5
chain1 = cg.Chain(f"{RDIR}/CG_c0001_v1.h5")
chain2 = cg.Chain(f"{RDIR}/CG_c0002_v1.h5")

wbands = [
    "K1", "Ka1", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"]
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

x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip_W = -0.233 * x - 2.222 * y + 2.504 * z



def set_rlabel(rlabel, x=0.995, y=0.925, fontsize=12):
    ax = plt.gca()
    if rlabel != "":
        rlabel = "$" + rlabel + "$"
    plt.text(
        x,
        y,
        rlabel,
        ha="right",
        va="center",
        fontsize=12,
        fontname=None,
        transform=ax.transAxes,
    )


def set_llabel(llabel, x=0.005, y=0.925, fontsize=12):
    ax = plt.gca()
    if llabel != "":
        llabel = "$" + llabel + "$"
    plt.text(
        x,
        y,
        llabel,
        ha="left",
        va="center",
        fontsize=15,
        fontname=None,
        transform=ax.transAxes,
    )


mu_Q = []
rms_Q = []

mu_V = []
rms_V = []

mu_W = []
rms_W = []
for n, b in enumerate(bands):
    print(b)

    ind1 = np.random.randint(burn_in, chain1.nsamples)
    ind2 = np.random.randint(burn_in, chain2.nsamples)
    m1 = chain1.get(f"tod/{b}/map", samples=ind1) * 1e3
    m2 = chain2.get(f"tod/{b}/map", samples=ind2) * 1e3
    #ms = np.concatenate((m1, m2))
    r1 = chain1.get(f"tod/{b}/rms", samples=ind1) * 1e6
    r2 = chain2.get(f"tod/{b}/rms", samples=ind2) * 1e6
    #rs = np.concatenate((r1, r2)) * 1e6
    if np.random.random() < 0.5:
       r = r1
    else:
       r = r2

    mu = hp.read_map(f"{RDIR}/CG_{b}_IQU_n0512_v1.fits", field=(0,1,2)) * 1e3
    rms = hp.read_map(f"{RDIR}/CG_{b}_IQU_n0512_v1.fits", field=(3,4,5,6))
    sd  = hp.read_map(f"{RDIR}/CG_{b}_IQU_n0512_v1.fits", field=(7,8,9,10))

    rms[:3] *= 1e3
    rms[3] *= 1e6
    sd[:3] *= 1e3
    sd[3] *= 1e6

    llabelT = ""
    llabelQ = ""
    llabelU = ""
    llabelQU = ""
    cbar = False
    if "023-WMAP_K" in b:
        llabelT = "T"
        llabelQ = "Q"
        llabelU = "U"
        llabelQU = r"\rho_{QU}"
    # if b == '090-WMAP_W4':
    #   cbar = True

    rlabel = r"\sigma_{" + b.split("_")[1] + r"}"
    if "040-WMAP_Q" in b:
        r /= 1.5**2
        fac = "/1.5"
    elif "060-WMAP_V" in b:
        r /= 2**2
        fac = "/2"
    elif "090-WMAP_W" in b:
        r /= 3**2
        fac = "/3"
    else:
        fac = ""
    cg.plot(
        r[0] ** 0.5,
        cmap="binary_r",
        min=25,
        max=60,
        width=2 * width,
        xsize=xsize,
        sub=(1, 4, 1),
        cbar=cbar,
    )
    set_llabel(llabelT)
    set_rlabel(rlabel + fac)
    cg.plot(
        r[1] ** 0.5,
        cmap="binary_r",
        min=35,
        max=85,
        width=2 * width,
        xsize=xsize,
        sub=(1, 4, 2),
        cbar=cbar,
    )
    set_rlabel(rlabel + fac)
    set_llabel(llabelQ)
    cg.plot(
        r[2] ** 0.5,
        cmap="binary_r",
        min=35,
        max=85,
        width=2 * width,
        xsize=xsize,
        sub=(1, 4, 3),
        cbar=cbar,
    )
    set_rlabel(rlabel + fac)
    set_llabel(llabelU)
    cg.plot(
        r[3] / np.sqrt(r[1] * r[2]),
        cmap="RdBu_r",
        min=-0.5,
        max=0.5,
        width=2 * width,
        xsize=xsize,
        sub=(1, 4, 4),
        cbar=cbar,
    )
    set_rlabel(rlabel)
    set_llabel(llabelQU)
    plt.tight_layout()
    plt.savefig(f"{b}_rms.pdf", bbox_inches="tight", dpi=dpi)
    plt.close("all")

    #mu = m1.mean(axis=0)
    #sd = m1.std(axis=0)
    mu_s = hp.smoothing(mu, fwhm=2 * np.pi / 180)

    rho_QU = sd[3]/sd[1]/sd[2]
    rlabel = r"\sigma_{" + b.split("_")[1] + r"}"
    # Limits:

    # K - (7, 8),    (1,3)
    # Ka - (2, 5),   (2,4)
    # Q  - (2, 5),   (4,8)
    # V - (3, 6),    (5,10)
    # W - (5, 20),   (8, 30)
    fI = ""
    fQ = ""
    fU = ""
    cbar = False
    llabelT = ""
    llabelQ = ""
    llabelU = ""
    llabelQU = ""
    if b == "023-WMAP_K":
        sd[0] /= 4
        f = "/4"
        fI = f
    elif b == "090-WMAP_W4":
        sd[0] /= 6
        sd[1:] /= 8
        fI = "/6"
        fQ = fU = "/8"
    elif b == "090-WMAP_W1":
        sd[0] /= 4
        sd[1] /= 6
        sd[2] /= 6
        fI = "/4"
        fP = "/6"
        fQ = fU = fP
    elif b == "090-WMAP_W2":
        sd[0] /= 4
        sd[1] /= 6
        sd[2] /= 6
        fI = "/4"
        fP = "/6"
        fQ = fU = fP
    elif b == "090-WMAP_W3":
        sd[0] /= 2
        fI = "/2"
        sd[1:] /= 3
        fQ = fU = "/3"
    #elif b == "030-WMAP_Ka":
    #    sd[0] /= 2
    #    f = "/2"
    #    fI = f
    #elif b == "040-WMAP_Q1":
    #    sd /= 2
    #    f = "/2"
    #    fI = fQ = fU = f
    #elif b == "040-WMAP_Q2":
    #    sd /= 2
    #    f = "/2"
    #    fI = fQ = fU = f
    elif "060-WMAP_V" in b:
        sd /= 2
        f = "/2"
        fI = fQ = fU = f

    cg.plot(
        sd,
        sig=0,
        fwhm=2 * u.deg,
        cmap="binary_r",
        width=width * 2,
        xsize=xsize,
        sub=(1, 4, 1),
        min=0,
        max=1,
        cbar=cbar,
        extend="both",
    )
    set_rlabel(rlabel + fI)
    set_llabel(llabelT)
    cg.plot(
        sd,
        sig=1,
        fwhm=2 * u.deg,
        cmap="binary_r",
        width=width * 2,
        xsize=xsize,
        sub=(1, 4, 2),
        cbar=cbar,
        min=0,
        max=1,
        extend="both",
    )
    set_rlabel(rlabel + fQ)
    set_llabel(llabelQ)
    cg.plot(
        sd,
        sig=2,
        fwhm=2 * u.deg,
        cmap="binary_r",
        width=width * 2,
        xsize=xsize,
        sub=(1, 4, 3),
        cbar=cbar,
        min=0,
        max=1,
        extend="both",
    )
    set_rlabel(rlabel + fU)
    set_llabel(llabelU)
    cg.plot(
        rho_QU,
        fwhm=2 * u.deg,
        cmap="RdBu_r",
        cbar=cbar,
        min=-0.5,
        max=0.5,
        sub=(1, 4, 4),
        width=width * 2,
        extend="both",
    )
    set_rlabel(rlabel)
    set_llabel(llabelQU)
    plt.tight_layout()
    plt.savefig(f"{b}_std.pdf", bbox_inches="tight", dpi=dpi)
    plt.close()

    rlabel = r"\langle\textit{" + b.split("_")[1] + r"}\rangle"
    cg.plot(
        mu,
        sig=0,
        rlabel=rlabel,
        llabel="T",
        unit=r"\mathrm{\mu K}",
        min=-3.4e3,
        max=3.4e3,
        width=width,
        xsize=xsize,
        extend="both",
    )
    plt.tight_layout()
    plt.savefig(f"{b}_map_I.pdf", bbox_inches="tight", dpi=300)
    plt.close('all')
    cg.plot(
        mu_s,
        sig=1,
        min=-30,
        max=30,
        width=width,
        xsize=xsize,
        rlabel=rlabel,
        llabel="Q",
        cbar=False,
    )
    plt.tight_layout()
    plt.savefig(f"{b}_map_Q.pdf", bbox_inches="tight", dpi=300)
    plt.close('all')
    cg.plot(
        mu_s,
        sig=2,
        min=-30,
        max=30,
        width=width,
        xsize=xsize,
        rlabel=rlabel,
        llabel="U",
        extend="both",
        unit=r"\mathrm{\mu K}",
    )
    plt.tight_layout()
    plt.savefig(f"{b}_map_U.pdf", bbox_inches="tight", dpi=300)
    plt.close('all')

    if ("023-WMAP_K" in b) or ("030-WMAP_Ka" in b):
        cg.plot(
            mu,
            sig=0,
            rlabel=rlabel,
            llabel="T",
            unit=r"\mathrm{\mu K}",
            min=-3.4e3,
            max=3.4e3,
            width=width,
            xsize=xsize,
            extend="both",
        )
        plt.savefig(f"{b}_mu_I.pdf", bbox_inches="tight", dpi=300)
        plt.close('all')
        cg.plot(
            mu_s,
            sig=1,
            min=-30,
            max=30,
            width=width,
            xsize=xsize,
            rlabel=rlabel,
            llabel="Q",
            cbar=False,
        )
        plt.savefig(f"{b}_mu_Q.pdf", bbox_inches="tight", dpi=300)
        plt.close('all')
        cg.plot(
            mu_s,
            sig=2,
            min=-30,
            max=30,
            width=width,
            xsize=xsize,
            rlabel=rlabel,
            llabel="U",
            extend="both",
            unit=r"\mathrm{\mu K}",
        )
        plt.savefig(f"{b}_mu_U.pdf", bbox_inches="tight", dpi=300)
        plt.close('all')
        plt.close("all")
    elif "040-WMAP_Q" in b:
        mu_Q.append(mu)
        rms_Q.append(r)
    elif "060-WMAP_V" in b:
        mu_V.append(mu)
        rms_V.append(r)
    elif "090-WMAP_W" in b:
        mu_W.append(mu)
        rms_W.append(r)

    #m1 = rng.choice(m1)
    #m2 = rng.choice(m2)
    diff = m1 - m2
    diff = hp.smoothing(diff, fwhm=5 * np.pi / 180)
    cg.plot(
        diff,
        sig=0,
        llabel=b.split("_")[1],
        rlabel=r"\Delta T",
        min=-3,
        max=3,
        remove_mono="true",
        sub=(1, 3, 1),
        cbar=False,
    )
    cg.plot(diff, sig=1, rlabel=r"\Delta Q", min=-3, max=3, sub=(1, 3, 2), cbar=False)
    cg.plot(diff, sig=2, rlabel=r"\Delta U", min=-3, max=3, sub=(1, 3, 3), cbar=False)
    plt.tight_layout()
    plt.savefig(f"{b}_sampdiff.pdf", bbox_inches="tight", dpi=dpi)
    plt.close()

    d_WMAP = hp.read_map(f"{WDIR}/wmap_iqusmap_r9_9yr_{wbands[n]}_v5.fits",
        field=(0,1,2))
    d_WMAP[0] += dip_W
    d_WMAP *= 1e3


Q = (mu_Q[0] / rms_Q[0][:3] + mu_Q[1] / rms_Q[1][:3]) / (
    1 / rms_Q[0][:3] + 1 / rms_Q[1][:3]
)
V = (mu_V[0] / rms_V[0][:3] + mu_V[1] / rms_V[1][:3]) / (
    1 / rms_V[0][:3] + 1 / rms_V[1][:3]
)
W = (
    mu_W[0] / rms_W[0][:3]
    + mu_W[1] / rms_W[1][:3]
    + mu_W[2] / rms_W[2][:3]
    + mu_W[3] / rms_W[3][:3]
) / (1 / rms_W[0][:3] + 1 / rms_W[1][:3] + 1 / rms_W[2][:3] + 1 / rms_W[3][:3])
Q_s = hp.smoothing(Q, fwhm=2 * np.pi / 180)
V_s = hp.smoothing(V, fwhm=2 * np.pi / 180)
W_s = hp.smoothing(W, fwhm=2 * np.pi / 180)


rlabel = r"\langle \textit{Q}\rangle"
cg.plot(
    Q,
    sig=0,
    rlabel=rlabel,
    llabel="T",
    unit=r"\mathrm{\mu K}",
    min=-3.4e3,
    max=3.4e3,
    width=width,
    xsize=xsize,
    extend="both",
)
plt.savefig(f"Q_mu_I.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    Q_s,
    sig=1,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="Q",
    cbar=False,
)
plt.savefig(f"Q_mu_Q.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    Q_s,
    sig=2,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="U",
    extend="both",
    unit=r"\mathrm{\mu K}",
)
plt.savefig(f"Q_mu_U.pdf", bbox_inches="tight", dpi=300)
plt.close("all")


rlabel = r"\langle \textit{V}\rangle"
cg.plot(
    V,
    sig=0,
    rlabel=rlabel,
    llabel="T",
    unit=r"\mathrm{\mu K}",
    min=-3.4e3,
    max=3.4e3,
    width=width,
    xsize=xsize,
    extend="both",
)
plt.savefig(f"V_mu_I.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    V_s,
    sig=1,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="Q",
    cbar=False,
)
plt.savefig(f"V_mu_Q.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    V_s,
    sig=2,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="U",
    extend="both",
    unit=r"\mathrm{\mu K}",
)
plt.savefig(f"V_mu_U.pdf", bbox_inches="tight", dpi=300)
plt.close("all")

rlabel = r"\langle \textit{W}\rangle"
cg.plot(
    W,
    sig=0,
    rlabel=rlabel,
    llabel="T",
    unit=r"\mathrm{\mu K}",
    min=-3.4e3,
    max=3.4e3,
    width=width,
    xsize=xsize,
    extend="both",
)
plt.savefig(f"W_mu_I.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    W_s,
    sig=1,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="Q",
    cbar=False,
)
plt.savefig(f"W_mu_Q.pdf", bbox_inches="tight", dpi=300)
plt.close('all')
cg.plot(
    W_s,
    sig=2,
    min=-30,
    max=30,
    width=width,
    xsize=xsize,
    rlabel=rlabel,
    llabel="U",
    extend="both",
    unit=r"\mathrm{\mu K}",
)
plt.savefig(f"W_mu_U.pdf", bbox_inches="tight", dpi=300)
plt.close('all')


cg.standalone_colorbar("binary_r", ticks=[1,2,3,4,], extend='both',
            unit=r"$\mathrm{\mu K_{CMB}}$",width=4, fontsize=18)
plt.savefig('cbar_std.pdf', dpi=300)
plt.close('all')
cg.standalone_colorbar("binary_r", ticks=[25, 35, 45, 60], extend='both',
            unit=r"$\mathrm{\mu K_{CMB}}$",width=4, fontsize=18)
plt.savefig('cbar_rms_I.pdf', dpi=300)
plt.close('all')
cg.standalone_colorbar("binary_r", ticks=[35, 50, 70, 85], extend='both',
            unit=r"$\mathrm{\mu K_{CMB}}$",width=4, fontsize=18)
plt.savefig('cbar_rms_P.pdf', dpi=300)
plt.close('all')
cg.standalone_colorbar("RdBu_r", ticks=[-0.5, 0,0.5], extend='both',
            width=4, fontsize=18, unit=r'\phantom{$\rho$}')
plt.savefig('cbar_rho.pdf', dpi=300)
plt.close('all')


