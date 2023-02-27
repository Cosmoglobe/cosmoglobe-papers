import matplotlib.pyplot as plt
import cosmoglobe as cg
import numpy as np

DIR = '/mn/stornext/d5/data/duncanwa/WMAP/v1'

cg.plot(f'{DIR}/CG_synch_IQU_n1024_v1.fits', comp='synch', llabel='T',
    rlabel=r'\langle A_\mathrm s\rangle', scale=1e-6,
    unit=r'$\mathrm{K_{RJ}}$', min=10, max=200, ticks=[10,100,200], width=4)
plt.savefig('synch_I.pdf', bbox_inches='tight')

cg.plot(f'{DIR}/CG_freefree_I_n1024_v1.fits', comp='ff', width=4,
    rlabel=r'\langle A_\mathrm{ff}\rangle', llabel='T')
plt.savefig('ff_I.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_ame_I_n1024_v1.fits', comp='ame', width=4, llabel='T',
    rlabel=r'\langle A_\mathrm{ame}\rangle')
plt.savefig('ame_I.pdf', bbox_inches='tight')
cg.plot(f'{DIR}/CG_dust_IQU_n1024_v1.fits', comp='dust', min=3, max=300,
    ticks=[3, 30, 300], width=4,llabel='T',
    rlabel=r'\langle A_\mathrm{d}\rangle')
plt.savefig('dust_I.pdf', bbox_inches='tight')
