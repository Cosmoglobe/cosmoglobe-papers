import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

rng = np.random.default_rng()

width = 8
xsize = 1200

DIR1 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_c_230103'
DIR2 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_d_230103'

chain1 = cg.Chain(f'{DIR1}/backup.h5')
chain2 = cg.Chain(f'{DIR2}/backup.h5')

bands = ['023-WMAP_K', '030-WMAP_Ka', '040-WMAP_Q1', '040-WMAP_Q2',
         '060-WMAP_V1', '060-WMAP_V2',
         '090-WMAP_W1', '090-WMAP_W2', '090-WMAP_W3', '090-WMAP_W4']

for b in bands:
    rlabel = r'\langle\textit{'+b.split('_')[1]+r'}\rangle'

    m1 = chain1.get(f'tod/{b}/map')
    m2 = chain2.get(f'tod/{b}/map')
    ms = np.concatenate((m1, m2))
    r1 = chain1.get(f'tod/{b}/rms')
    r2 = chain2.get(f'tod/{b}/rms')
    rs = np.concatenate((r1, r2))
    r = rng.choice(rs)


    if len(r) == 3:
        P = np.concatenate((r[1], r[2]))
        minval_P, maxval_P = np.percentile(P*1e3, [5, 95])
        minval_T, maxval_T = np.percentile(r[0]*1e3, [5,95])
        cg.plot(r*1e3, sig=0, cmap='binary_r', min=np.round(minval_T,-1),
            max=np.round(maxval_T, -1), sub=(1,4,1))
        cg.plot(r*1e3, sig=1, cmap='binary_r', min=np.round(minval_P,-1),
            max=np.round(maxval_P, -1), sub=(1,4,2))
        cg.plot(r*1e3, sig=2, cmap='binary_r', min=np.round(minval_P,-1),
            max=np.round(maxval_P, -1), sub=(1,4,3))
    else:
        P = np.concatenate((r[1], r[2]))
        minval_P, maxval_P = np.percentile(P**0.5*1e3, [5, 95])
        minval_T, maxval_T = np.percentile(r[0]**0.5*1e3, [5,95])
        cg.plot(r**0.5*1e3, sig=0, cmap='binary_r', min=np.round(minval_T,-1),
            max=np.round(maxval_T,-1),
            sub=(1,4,1))
        cg.plot(r**0.5*1e3, sig=1, cmap='binary_r', min=np.round(minval_P,-1),
            max=np.round(maxval_P,-1),
            sub=(1,4,2))
        cg.plot(r**0.5*1e3, sig=2, cmap='binary_r', min=np.round(minval_P,-1),
            max=np.round(maxval_P,-1),
            sub=(1,4,3))
        cg.plot(r[3]/np.sqrt(r[1]*r[2]), cmap='RdBu_r', min=-0.5, max=0.5,
            sub=(1,4,4))
    plt.savefig(f'{b}_std.pdf', bbox_inches='tight')
    plt.close('all')
    
    mu = m1.mean(axis=0)
    sd = m1.std(axis=0)
    mu_s = hp.smoothing(mu, fwhm=2*np.pi/180)

    rho_QU = ((ms[:,1]-mu[1])*(ms[:,2]-mu[2])).mean(axis=0)/sd[1]/sd[2]
    cg.plot(sd, sig=0, llabel='\sigma_T', fwhm=2*u.deg, cmap='binary_r',
        width=width, xsize=xsize, sub=(1,4,1), scale=1e3, min=1, max=3)
    cg.plot(sd, sig=1, llabel='\sigma_Q', fwhm=2*u.deg, cmap='binary_r',
        width=width, xsize=xsize, sub=(1,4,2), min=1, max=3, scale=1e3)
    cg.plot(sd, sig=2, llabel='\sigma_U', fwhm=2*u.deg, cmap='binary_r',
        width=width, xsize=xsize, sub=(1,4,3), min=1, max=3, scale=1e3)
    cg.plot(rho_QU, fwhm=2*u.deg, llabel=r'c_{QU}', cmap='RdBu_r',
        min=-1, max=1, sub=(1,4,4))
    plt.savefig(f'{b}_rms.pdf', bbox_inches='tight')
    plt.close()

    cg.plot(mu, sig=0, rlabel=rlabel, llabel='T',
        min=-3.4, max=3.4, width=width, xsize=xsize, sub=(3,1,1))
    cg.plot(mu_s, sig=1, min=-30e-3, max=30e-3, width=width,
        xsize=xsize, rlabel=rlabel, llabel='Q', sub=(3,1,2))
    cg.plot(mu_s, sig=2, min=-30e-3, max=30e-3, width=width,
        xsize=xsize, rlabel=rlabel, llabel='U', sub=(3,1,3))
    plt.tight_layout()
    plt.savefig(f'{b}_map.pdf', bbox_inches='tight')
    plt.close()



    m1 = rng.choice(m1)
    m2 = rng.choice(m2)
    diff = m1 - m2
    diff = hp.smoothing(diff, fwhm=5*np.pi/180)
    cg.plot(diff*1e3, sig=0, llabel=b.split('_')[1], rlabel=r'\Delta T',
        min=-3, max=3, remove_mono='true', sub=(1,3,1), cbar=False)
    cg.plot(diff*1e3, sig=1, rlabel=r'\Delta Q',
        min=-3, max=3, sub=(1,3,2), cbar=False)
    cg.plot(diff*1e3, sig=2, rlabel=r'\Delta U',
        min=-3, max=3, sub=(1,3,3), cbar=False)
    plt.tight_layout()
    plt.savefig(f'{b}_sampdiff.pdf', bbox_inches='tight')
    plt.show()
