import numpy as np
import matplotlib.pyplot as plt
import cosmoglobe as cg
import healpy as hp
from glob import glob



WMAP9 = '/mn/stornext/d16/cmbco/ola/wmap/foreground_products'
#for mcmc in ['c', 'e', 'f', 'g']:
for mcmc in ['c']:
    fnames = glob(f'{WMAP9}/mcmc_{mcmc}/*synch_stk*')
    fnames.sort()
    synch_q = hp.read_map(fnames[0], field=(0,1,2))
    synch_u = hp.read_map(fnames[1], field=(0,1,2))
    cg.plot(1e3*np.hypot(synch_q[0], synch_u[0])*0.42, min=0, max=50, comp='synch',
        sig=0, llabel=r'\mathit{WMAP9}', rlabel='P', unit=r'\mathrm{\mu K}')
    plt.savefig('polint_WMAP9.png', bbox_inches='tight')
    cg.plot(1e3*synch_q**0.5*0.42, sig=2, llabel=r'\sigma_P', min=0, max=5, cmap='binary_r')
    #plt.savefig('polint_WMAP9_sigma.png', bbox_inches='tight')


DIR_Planck18 = '/mn/stornext/d16/cmbco/ola/planck_products/commander'

synch, header = hp.read_map(f'{DIR_Planck18}/COM_CompMap_QU-synchrotron-commander_2048_R3.00_full.fits',
    field=(0,1), h=True)

cg.plot(np.hypot(synch[0], synch[1]), min=0, max=50, comp='synch', 
    llabel=r'\mathrm{PR3}', rlabel='P')
plt.savefig('polint_PR3.png', bbox_inches='tight')

DIR_PR4 = '/mn/stornext/d16/cmbco/ola/npipe'

synch, header = hp.read_map(f'{DIR_PR4}/npipe6v20_comm_synch_n2048_40arc_QU_rc1.fits',
    field=(1,2), h=True)

cg.plot(np.hypot(synch[0], synch[1]), min=0, max=50, comp='synch', 
    llabel=r'\mathrm{PR4}', rlabel='P')
plt.savefig('polint_PR4.png', bbox_inches='tight')

DIR_BP = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2'

synch_BP, h = hp.read_map(f'{DIR_BP}/BP_synch_IQU_n1024_v2.fits',
    field=(0,1,2,7,8), h=True)

cg.plot(np.hypot(synch_BP[1], synch_BP[2]), min=0, max=50, comp='synch', sig=0,
    llabel=r'\textsc{BeyondPlanck}', rlabel='P', unit=r'\mathrm{\mu K}')
plt.savefig('polint_BP.png', bbox_inches='tight')
cg.plot(np.hypot(synch_BP[3], synch_BP[4]), min=0, max=10, cmap='binary_r',
    llabel=r'\textsc{BeyondPlanck}', rlabel='\sigma_P')
#plt.savefig('polint_BP_sigma.png', bbox_inches='tight')

DIR_CG = '/mn/stornext/d16/cmbco/cg/v1'

synch_CG, h = hp.read_map(f'{DIR_CG}/CG_synch_IQU_n1024_v1.fits',
    field=(0,1,2,7,8), h=True)

cg.plot(np.hypot(synch_CG[1], synch_CG[2]), min=0, max=50, comp='synch', sig=0,
    llabel=r'\textsc{Cosmoglobe}', rlabel='P', unit=r'\mathrm{\mu K}')
plt.savefig('polint_CG.png', bbox_inches='tight')
cg.plot(np.hypot(synch_CG[3], synch_CG[4]), min=0, max=10, cmap='binary_r',
    llabel=r'\textsc{Cosmoglobe}', rlabel='\sigma_P')
#plt.savefig('polint_CG_sigma.png', bbox_inches='tight')
