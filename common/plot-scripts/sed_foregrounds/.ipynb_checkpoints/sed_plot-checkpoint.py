import cosmoglobe as cg

import astropy.units as u
import astropy.constants as c
#from astropy import visualization
#visualization.quantity_support()

import numpy as np
import matplotlib.pyplot as plt



def dBdT(nu, T):
    x = c.h*nu/(c.k_B*T)
    A = 2*c.k_B*nu**2/c.c**2
    return (x.si**2*np.exp(x.si)/np.expm1(x.si)**2).to('')


def plot_all(real_units=True,
        plot_ff=True,
        plot_sd=True,
        plot_synch=True,
        plot_dust=True,
        plot_dust2=True,
        plot_zodi=True,
        plot_stars=True,
        plot_cib=True,
        plot_cmb=True):
    nu = np.geomspace(0.4, 3_000_000, 1000)*u.GHz
    #nu = np.geomspace(0.4, 3_000, 1000)*u.GHz
    ff = model.components['ff'].get_freq_scaling(nu, T_e=8000*u.K)*u.uK
    synch = model.components['synch'].get_freq_scaling(nu, -3.1)*u.uK
    sd = model.components['ame'].get_freq_scaling(nu,20*u.GHz)*u.uK
    dust = model.components['dust'].get_freq_scaling(nu, 1.6, 15*u.K)*u.uK
    dust2 = model.components['dust'].get_freq_scaling(nu, 2.0, 30*u.K)*u.uK
    zodi = model.components['dust'].get_freq_scaling(nu, 0, 200*u.K)*u.uK
    cmb = model.components['dust'].get_freq_scaling(nu, 0, 2.7*u.K)*u.uK
    zodi_scat = model.components['dust'].get_freq_scaling(nu, 0, 5700*u.K)*u.uK

    cold_stars = model.components['dust'].get_freq_scaling(nu, 0, 2000*u.K)*u.uK

    if real_units:
       C = (2*c.k_B*nu**2/c.c**2).to('MJy/uK')
    else:
       C = 1
    C_bla = (c.c**2/(2*c.k_B*nu**2)).to('uK/MJy')

    #beta_earth = (29.78*u.km/u.s/c.c).to('')
    #beta_sun = (250*u.km/u.s/c.c).to('')
    #plt.plot(nu, C*cmb[0], label='CMB Monopole')
    #plt.plot(nu, C*(3.5e3*u.uK*dBdT(nu, 2.7275*u.K)), label='Solar Dipole')
    #plt.plot(nu, C*(0.35e3*u.uK*dBdT(nu, 2.725*u.K)), label='Orbital Dipole')
    if plot_cmb: plt.plot(nu, C*dBdT(nu, 2.7275*u.K)*75*u.uK, label='CMB fluctuations', color='k')

    if plot_ff: plt.loglog(nu, (C*ff[0]*25e0*np.exp(-((c.h*nu)/(c.k_B*1e4*u.K)).to('').value)), label='Free-free', color='xkcd:green')
    if plot_synch: plt.loglog(nu, C*synch[0]*5e6, label='Synchrotron', color='xkcd:azure')
    plt.ylim([1e-9, 1e5])
    if C == 1:
        plt.ylim([1e-1, 1e5])


    if plot_sd: plt.plot(nu, C*sd[0]*25, label='AME', color='xkcd:sienna')
    if plot_dust: plt.plot(nu, C*dust[0]*75, label='Cold dust', color='xkcd:salmon')
    if plot_dust2: plt.plot(nu, C*dust2[0]*75, label='Hot dust', color='xkcd:red')
    if plot_stars: plt.plot(nu, C*cold_stars[0]*2e-2, label='Cold stars')
    #todo: replace this with real star SED

    sed_bla1 = np.load('sed_dirbe_zodi.npy')
    wav_bla = np.load('nu_dirbe_zodi.npy')*u.micron
    nu_bla1 = (c.c/wav_bla).to('GHz')
    sed_bla2 =np.load('sed_planck_zodi.npy')
    wav_bla = np.load('nu_planck_zodi.npy')*u.micron
    nu_bla2 = (c.c/wav_bla).to('GHz')

    nu_bla = np.concatenate([nu_bla1, nu_bla2])
    if real_units:
        C = (2*c.k_B*nu_bla**2/c.c**2).to('MJy/uK')
    else:
        C = 1
    C_bla = (c.c**2/(2*c.k_B*nu_bla**2)).to('uK/MJy')
    if plot_zodi: plt.plot(nu_bla, C*C_bla*np.concatenate([sed_bla1, sed_bla2]), label='Zodi')
    data = np.loadtxt('CIB_bethermin.txt')
    inds = np.argsort(data[:,0])
    wav = data[inds,0]
    nu_bla2 = (c.c/(wav*u.micron)).to('GHz')
    nuInu = data[inds,1]*u.nW/u.m**2/u.sr
    I_nu = (nuInu/nu_bla2).to('MJy sr-1')

    if real_units:
        C = (2*c.k_B*nu_bla2**2/c.c**2).to('MJy/uK')
    else:
        C = 1
    C_bla = (c.c**2/(2*c.k_B*nu_bla2**2)).to('uK/MJy')
    if plot_cib: plt.loglog(nu_bla2, C*C_bla*I_nu, label='CIB (Béthermin 2013)', color='b')
    

    #add dirbe bands as vertical lines
    dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25, 60, 100, 140, 240]*u.micron #microns
    plt.vlines(dirbe_bands.to('GHz', equivalencies=u.spectral()), 1e-5, 1e3, color='pink', alpha=0.5)

    dirbe_names = ['1', '2','3', '4', '5', '6', '7', '8', '9', '10']
    locs = [250000, 140000, 90000, 65000, 27000, 12500, 5200, 3100, 2200, 1300]
    for loc, name in zip(locs, dirbe_names):
        plt.text(loc, 400, name, color='pink')

    plt.text(400000, 400, 'DIRBE Bands', color='pink')
    #add planck bands
    hfi_bands = [100, 143, 217, 353, 545, 857]*u.GHz

    hfi_names = ['100', '143', '217' , '353', '545', '857']
    plt.vlines(hfi_bands, 1e-5, 1e3, color='red', alpha=0.3)
    
    hfi_locs = [101, 145, 220, 360, 550, 870]

    for loc, name in zip(hfi_locs, hfi_names):
        plt.text(loc, 450, name, color='red', alpha=0.8, fontsize=6)

    plt.text(20, 400, 'HFI Bands', color='red', alpha=0.6)

    plt.legend(facecolor='white')

    #plt.ylim(ymin=1e-3)
    plt.xlabel(r'GHz')
    if C == 1:
        plt.legend(facecolor='white', loc = 'right')
        plt.ylabel('uK')
        plt.ylim([1e-3, 1e3])
    else:
        plt.legend(facecolor='white')
        plt.ylabel(r'MJy/sr')
        plt.ylim([1e-5, 1e3])
    xlim = np.array([nu[0].value, nu[-1].value])
    xlim2 = (c.c/(xlim*u.GHz)).to('micron').value
    
    ax = plt.gca()
    ax.set_xlim(xlim)
    ax2 = ax.twiny()
    ax2.set_xlim(xlim2)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$\mathrm{\mu m}$')


model = cg.sky_model(nside=1)

plt.figure(figsize=(12, 4))
plot_all(real_units=False)
plt.savefig('all_fgs_uK.pdf', bbox_inches='tight')
plt.figure(figsize=(12, 4))
plot_all(real_units=True,)
plt.savefig('all_fgs.pdf', bbox_inches='tight')
