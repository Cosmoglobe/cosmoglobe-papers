import healpy as hp
import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt

DIR ="/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_a_230206/"
d_Ka = hp.read_map(f"{DIR}/tod_030-WMAP_Ka_rms_c0001_k000011.fits", field=(1, 2, 3))
d_30 = hp.read_map(f"{DIR}/tod_030_rms_c0001_k000011.fits", field=(1, 2, 3))
d_44 = hp.read_map(f"{DIR}/tod_044_rms_c0001_k000011.fits", field=(1, 2, 3))

rho_Ka = d_Ka[2] / np.sqrt(d_Ka[0] * d_Ka[1])
rho_30 = d_30[2] / np.sqrt(d_30[0] * d_30[1])
rho_44 = d_44[2] / np.sqrt(d_44[0] * d_44[1])

cmap = "RdBu_r"
cg.plot(
    rho_Ka,
    min=-0.5,
    max=0.5,
    llabel="\mathit{Ka}",
    #rlabel=r"\rho_{QU}",
    cmap=cmap,
    width=4,
    xsize=1000,
    cbar=False,
)
plt.savefig("rho_QU_Ka.pdf", dpi=100, bbox_inches='tight')
cg.plot(
    rho_30,
    min=-0.5,
    max=0.5,
    llabel="30",
    #rlabel=r"\rho_{QU}",
    cmap=cmap,
    width=4,
    xsize=1000,
    cbar=True,
    extend='both',
)
plt.savefig("rho_QU_30.pdf", dpi=100, bbox_inches='tight')
cg.plot(
    rho_44,
    min=-0.5,
    max=0.5,
    llabel="44",
    #rlabel=r"\rho_{QU}",
    cmap=cmap,
    width=4,
    xsize=1000,
    cbar=True,
)
plt.savefig("rho_QU_44.pdf", dpi=100, bbox_inches='tight')

#cg.standalone_colorbar("RdBu_r", ticks=[-0.1, 0,0.1], extend='both',
#            width=3, fontsize=18, unit=r'\phantom{$\rho$}')
#plt.savefig("cbar_rho_01.pdf", dpi=300)
