import cosmoglobe
import matplotlib.pyplot as plt

DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

chain = cosmoglobe.Chain(f"{DIR}/CG_c0001_v0.h5")
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
labs = [
    r"\textit K",
    r"\textit{Ka}",
    r"\textit Q1",
    r"\textit Q2",
    r"\textit V1",
    r"\textit V2",
    r"\textit W1",
    r"\textit W2",
    r"\textit W3",
    r"\textit W4",
]
for l, b in zip(labs, bands):
    gain = chain.get(f"tod/{b}/gain")
    xi_n = chain.get(f"tod/{b}/xi_n")
    time = chain.get(f"tod/{b}/MJD")
    accept = chain.get(f"tod/{b}/accept")
    chisq = chain.get(f"tod/{b}/chisq")
    chisq_mu = chisq.mean(axis=0)
    fig, axes = plt.subplots(sharex=True, nrows=5, figsize=(6, 6))
    g_mu = gain.mean(axis=0)
    xi_mu = xi_n.mean(axis=0)
    t = time[-1]
    inds = accept[-1][0] == 1
    for i in range(4):
        axes[4].plot(t[inds], chisq_mu[i, inds], color=plt.cm.viridis(i / 3))
        axes[3].plot(t[inds], xi_mu[2, i, inds], color=plt.cm.viridis(i / 3))
        axes[2].plot(t[inds], 1e3 * xi_mu[1, i, inds], color=plt.cm.viridis(i / 3))
        axes[1].plot(
            t[inds], xi_mu[0, i, inds] / g_mu[i, inds], color=plt.cm.viridis(i / 3)
        )
        axes[0].plot(t[inds], g_mu[i, inds], color=plt.cm.viridis(i / 3))
    axes[0].set_ylabel(r"$g\ [\mathrm{du\,mK^{-1}}]$")
    axes[1].set_ylabel(r"$\sigma_0$ [mK]")
    axes[2].set_ylabel(r"$f_\mathrm{knee}$ [mHz]")
    axes[3].set_ylabel(r"$\alpha$")
    axes[4].set_ylabel(r"$\chi^2$")
    axes[4].set_xlabel(r"Time [MJD]")
    plt.suptitle(l, y=0.95, size=16)
    plt.savefig(f"inst_{b}.pdf", bbox_inches="tight")
plt.show()
