import cosmoglobe as cg
import astropy.units as u
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np

rng = np.random.default_rng()

width = 6
xsize = 1200

DIR1 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_a_230104'
DIR2 = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_LFI_KKaQVW_b_230104'

chain1 = cg.Chain(f'{DIR1}/backup.h5', burn_in=2)
chain2 = cg.Chain(f'{DIR2}/backup.h5', burn_in=2)

bands = ['023-WMAP_K', '030-WMAP_Ka', '040-WMAP_Q1', '040-WMAP_Q2',
         '060-WMAP_V1', '060-WMAP_V2',
         '090-WMAP_W1', '090-WMAP_W2', '090-WMAP_W3', '090-WMAP_W4']

def set_rlabel(rlabel,x=0.995, y=0.925, fontsize=12):
    ax = plt.gca()
    if rlabel != '':
      rlabel = '$'+rlabel+'$'
    plt.text(
        x,
        y,
        rlabel,
        ha="right",
        va="center",
        fontsize=15,
        fontname=None,
        transform=ax.transAxes,
    )
def set_llabel(llabel,x=0.005, y=0.925, fontsize=12):
    ax = plt.gca()
    if llabel != '':
      llabel = '$'+llabel+'$'
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

for b in bands:
    print(b)

    m1 = chain1.get(f'tod/{b}/map')*1e3
    m2 = chain2.get(f'tod/{b}/map')*1e3
    ms = np.concatenate((m1, m2))
    r1 = chain1.get(f'tod/{b}/rms')
    r2 = chain2.get(f'tod/{b}/rms')
    rs = np.concatenate((r1, r2))*1e6
    r = rng.choice(rs)

    llabelT = ''
    llabelQ = ''
    llabelU = ''
    llabelQU = ''
    cbar = False
    if '023-WMAP_K' in b:
       llabelT = 'T'
       llabelQ = 'Q'
       llabelU = 'U'
       llabelQU = r'\rho_{QU}'
    #if b == '090-WMAP_W4':
    #   cbar = True


    rlabel=r'A^\mathrm{RMS}_{' + b.split('_')[1] + r'}'
    if '040-WMAP_Q' in b:
        r /= (1.5**2)
        fac = '/1.5'
    elif '060-WMAP_V' in b:
        r /= (2**2)
        fac = '/2'
    elif '090-WMAP_W' in b:
        r /= (3**2)
        fac = '/3'
    else:
        fac = ''
    cg.plot(r[0]**0.5, cmap='binary_r',
        min=25, max=60, width=2*width, xsize=xsize,
        sub=(1,4,1), cbar=cbar)
    set_llabel(llabelT)
    set_rlabel(rlabel+fac)
    cg.plot(r[1]**0.5,  cmap='binary_r',
        min=35, max=85, width=2*width, xsize=xsize,
        sub=(1,4,2), cbar=cbar)
    set_rlabel(rlabel+fac)
    set_llabel(llabelQ)
    cg.plot(r[2]**0.5, cmap='binary_r',
        min=35, max=85, width=2*width, xsize=xsize,
        sub=(1,4,3), cbar=cbar)
    set_rlabel(rlabel+fac)
    set_llabel(llabelU)
    cg.plot(r[3]/np.sqrt(r[1]*r[2]), cmap='RdBu_r', min=-0.5, max=0.5,
        width=2*width, xsize=xsize,
        sub=(1,4,4), cbar=cbar)
    set_rlabel(rlabel)
    set_llabel(llabelQU)
    plt.tight_layout()
    plt.savefig(f'{b}_rms.pdf', bbox_inches='tight')
    plt.close('all')
    
    mu = m1.mean(axis=0)
    sd = m1.std(axis=0)
    mu_s = hp.smoothing(mu, fwhm=2*np.pi/180)

    rho_QU = ((ms[:,1]-mu[1])*(ms[:,2]-mu[2])).mean(axis=0)/sd[1]/sd[2]
    rlabel=r'\sigma_{' + b.split('_')[1] + r'}'
    # Limits:

    # K - (7, 8),    (1,3)
    # Ka - (2, 5),   (2,4)
    # Q  - (2, 5),   (4,8)
    # V - (3, 6),    (5,10)
    # W - (5, 20),   (8, 30)
    fI = ''
    fQ = ''
    fU = ''
    cbar = False
    llabelT = ''
    llabelQ = ''
    llabelU = ''
    llabelQU = ''
    if '023-WMAP_K' in b:
       lim_T = (1,5)
       lim_P = (1,5)
       llabelT = 'T'
       llabelQ = 'Q'
       llabelU = 'U'
       llabelQU = r'\rho_{QU}'
    elif '030-WMAP_Ka' in b:
       lim_T = (1,5)
       lim_P = (1,5)
    elif '040-WMAP_Q' in b:
       lim_T = (1,5)
       lim_P = (1,5)
    elif '060-WMAP_V' in b:
       lim_T = (1,5)
       lim_P = (1,5)
    elif '090-WMAP_W' in b:
       lim_T = (1,5)
       lim_P = (1,5)
    if b == '090-WMAP_W4':
       sd[0] /= 8
       sd[1:] /= 8
       fI = '/8'
       fQ = fU = '/8'
       cbar = True
    elif b == '090-WMAP_W1':
       sd /= 4
       f = '/4'
       fI = fQ = fU = f
    elif b == '090-WMAP_W2':
       sd /= 4
       f = '/4'
       fI = fQ = fU = f
    elif b == '090-WMAP_W3':
       sd[0] /= 2
       fI = '/2'
       sd[1:] /= 3
       fQ = fU = '/3'
    elif b == '023-WMAP_K':
       sd[0] /= 2
       f = '/2'
       fI = f
    elif b == '040-WMAP_Q2':
       sd /= 2
       f = '/2'
       fI = fQ = fU = f
    elif '060-WMAP_V' in b:
       sd /= 2
       f = '/2'
       fI = fQ = fU = f

    cg.plot(sd, sig=0,fwhm=2*u.deg, cmap='binary_r',
        width=width*2, xsize=xsize, sub=(1,4,1),
        min=1, max=4,cbar=cbar,
        )
    set_rlabel(rlabel+fI)
    set_llabel(llabelT)
    cg.plot(sd, sig=1, fwhm=2*u.deg, cmap='binary_r',
        width=width*2, xsize=xsize, sub=(1,4,2), cbar=cbar,
        min=1, max=4)
    set_rlabel(rlabel+fQ)
    set_llabel(llabelQ)
    cg.plot(sd, sig=2, fwhm=2*u.deg, cmap='binary_r',
        width=width*2, xsize=xsize, sub=(1,4,3), cbar=cbar,
        min=1, max=4)
    set_rlabel(rlabel+fU)
    set_llabel(llabelU)
    cg.plot(rho_QU, fwhm=2*u.deg, cmap='RdBu_r',cbar=cbar,
        min=-0.5, max=0.5, sub=(1,4,4), width=width*2)
    set_rlabel(rlabel)
    set_llabel(llabelQU)
    plt.tight_layout()
    plt.savefig(f'{b}_std.pdf', bbox_inches='tight')
    plt.close()

    rlabel = r'\langle\textit{'+b.split('_')[1]+r'}\rangle'
    cg.plot(mu, sig=0, rlabel=rlabel, llabel='T', unit=r'\mathrm{\mu K}',
        min=-3.4e3, max=3.4e3, width=width, xsize=xsize, sub=(3,1,1), extend='both')
    cg.plot(mu_s, sig=1, min=-30, max=30, width=width,
        xsize=xsize, rlabel=rlabel, llabel='Q', sub=(3,1,2), cbar=False)
    cg.plot(mu_s, sig=2, min=-30, max=30, width=width,
        xsize=xsize, rlabel=rlabel, llabel='U', sub=(3,1,3), extend='both',
        unit=r'\mathrm{\mu K}')
    plt.tight_layout()
    plt.savefig(f'{b}_map.pdf', bbox_inches='tight')
    plt.close()



    m1 = rng.choice(m1)
    m2 = rng.choice(m2)
    diff = m1 - m2
    diff = hp.smoothing(diff, fwhm=5*np.pi/180)
    cg.plot(diff, sig=0, llabel=b.split('_')[1], rlabel=r'\Delta T',
        min=-3, max=3, remove_mono='true', sub=(1,3,1), cbar=False)
    cg.plot(diff, sig=1, rlabel=r'\Delta Q',
        min=-3, max=3, sub=(1,3,2), cbar=False)
    cg.plot(diff, sig=2, rlabel=r'\Delta U',
        min=-3, max=3, sub=(1,3,3), cbar=False)
    plt.tight_layout()
    plt.savefig(f'{b}_sampdiff.pdf', bbox_inches='tight')
    plt.close()
