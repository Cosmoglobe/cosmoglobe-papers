import numpy as np
import matplotlib.pyplot as plt
import cosmoglobe as cg
import healpy as hp
from glob import glob

w = 3.5
fontsize = {'llabel':9, 'rlabel':9}
dpi = 150

CG_DIR = '/mn/stornext/d5/data/duncanwa/WMAP'

#HM1 = hp.read_map(f'{CG_DIR}/CG_HM1/CG_synch_IQU_n1024_CG_HM1.fits',
#    field=(0,1,2))
#HM2 = hp.read_map(f'{CG_DIR}/CG_HM2/CG_synch_IQU_n1024_CG_HM2.fits',
#    field=(0,1,2))

DIR_CG = '/mn/stornext/d16/cmbco/cg/v1'

synch_CG, h = hp.read_map(f'{DIR_CG}/CG_synch_IQU_n1024_v1.fits',
    field=(0,1,2,7,8), h=True)
DIR_BP = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'

synch_BP, h = hp.read_map(f'{DIR_BP}/BP_synch_IQU_n1024_v2.fits',
    field=(0,1,2,7,8), h=True)

DIR_Planck18 = '/mn/stornext/d16/cmbco/ola/planck_products/commander'
DIR_PR4 = '/mn/stornext/d16/cmbco/ola/npipe'

#synch, header = hp.read_map(f'{DIR_PR4}/npipe6v20_comm_synch_n2048_40arc_QU_rc1.fits',
#    field=(1,2), h=True)
#
#synch, header = hp.read_map(f'{DIR_Planck18}/COM_CompMap_QU-synchrotron-commander_2048_R3.00_full.fits',
#    field=(0,1), h=True)
#polint = np.sqrt(HM1[1]*HM2[1] + HM1[2]*HM2[2])
'''
cg.plot(polint, min=0, max=50, comp='synch',
    rlabel=r'P^*', llabel=r'\textsc{Cosmoglobe}')
plt.savefig('../figures/polint_CG_cross.pdf', bbox_inches='tight')

WMAP9 = '/mn/stornext/d16/cmbco/ola/wmap/foreground_products'
#for mcmc in ['c', 'e', 'f', 'g']:
for mcmc in ['c']:
    fnames = glob(f'{WMAP9}/mcmc_{mcmc}/*synch_stk*')
    fnames.sort()
    synch_q = hp.read_map(fnames[0], field=(0,1,2))
    synch_u = hp.read_map(fnames[1], field=(0,1,2))
    cg.plot(1e3*np.hypot(synch_q[0], synch_u[0])*0.42, min=0, max=50, comp='synch',
        sig=0, llabel=r'\mathit{WMAP9}', rlabel='P', unit=r'\mathrm{\mu K}')
    plt.savefig('../figures/polint_WMAP9.pdf', bbox_inches='tight')
    cg.plot(1e3*synch_q**0.5*0.42, sig=2, llabel=r'\sigma_P', min=0, max=5, cmap='binary_r')
    #plt.savefig('polint_WMAP9_sigma.png', bbox_inches='tight')



cg.plot(np.hypot(synch[0], synch[1]), min=0, max=50, comp='synch', 
    llabel=r'\mathrm{PR3}', rlabel='P')
plt.savefig('../figures/polint_PR3.pdf', bbox_inches='tight')


cg.plot(np.hypot(synch[0], synch[1]), min=0, max=50, comp='synch', 
    llabel=r'\mathrm{PR4}', rlabel='P')
plt.savefig('../figures/polint_PR4.pdf', bbox_inches='tight')


cg.plot(np.hypot(synch_BP[1], synch_BP[2]), min=0, max=50, comp='synch', sig=0,
    llabel=r'\textsc{BeyondPlanck}', rlabel='P', unit=r'\mathrm{\mu K}')
plt.savefig('../figures/polint_BP.pdf', bbox_inches='tight')
cg.plot(np.hypot(synch_BP[3], synch_BP[4]), min=0, max=10, cmap='binary_r',
    llabel=r'\textsc{BeyondPlanck}', rlabel='\sigma_P', unit=r'\mathrm{\mu K}')
plt.savefig('../figures/polint_BP_sigma.pdf', bbox_inches='tight')

cg.plot(np.hypot(synch_BP[3], synch_BP[4]), min=0, max=10, cmap='binary_r',
    llabel=r'\textsc{BeyondPlanck}', rlabel='\sigma_P', unit=r'\mathrm{\mu K}')

'''

cg.plot(np.hypot(synch_CG[1], synch_CG[2]), min=0, max=50, comp='synch', sig=0,
        norm='linear',
    llabel=r'\mathrm{Cosmoglobe}', rlabel='P', unit=r'\mathrm{\mu K}', width=w,
    fontsize=fontsize)

fig = plt.gcf()
plt.sca(fig.axes[0])

theta0, phi0 = 53*np.pi/180, (-100)*np.pi/180
topline_phi = np.linspace(phi0-5*np.pi/180, phi0+5*np.pi/180)
topline_theta = np.ones_like(topline_phi)*(theta0-5*np.pi/180)
hp.newprojplot(theta=topline_theta, phi=topline_phi, color='r');

topline_theta = np.ones_like(topline_phi)*(theta0+5*np.pi/180)
hp.newprojplot(theta=topline_theta, phi=topline_phi, color='r');

leftline_theta = np.linspace(theta0 - np.pi/180*5, theta0 + np.pi/180*5)
leftline_phi = np.ones_like(leftline_theta)*(phi0 + np.pi/180*5)

hp.newprojplot(theta=leftline_theta, phi=leftline_phi, color='r');
leftline_phi = np.ones_like(leftline_theta)*(phi0 - np.pi/180*5)
hp.newprojplot(theta=leftline_theta, phi=leftline_phi, color='r');

plt.savefig('../figures/polint_CG.pdf', bbox_inches='tight', dpi=dpi)
#cg.plot(np.hypot(synch_CG[3], synch_CG[4]), min=0, max=5, cmap='binary_r',
#    llabel=r'\textsc{Cosmoglobe}', rlabel='\sigma_P', unit=r'\mathrm{\mu K}')
#plt.savefig('../figures/polint_CG_sigma.pdf', bbox_inches='tight')

cg.plot(np.hypot(synch_CG[3], synch_CG[4])/np.hypot(synch_BP[3], synch_BP[4]),
        width=w,
    cmap='bone_r',
    min=0.4, max=1,
    llabel=r'\mathrm{Ratio}',
    rlabel='\sigma_P^\mathrm{CG}/\sigma_P^\mathrm{BP}', fontsize=fontsize)
plt.savefig('../figures/polint_sigma_ratio.pdf', bbox_inches='tight', dpi=dpi)

plt.close('all')
plt.hist(np.hypot(synch_CG[3], synch_CG[4])/np.hypot(synch_BP[3], synch_BP[4]),
    100)
plt.yscale('log')
#plt.show()
