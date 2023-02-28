from setup_matplotlib import *

import cosmoglobe as cg
import matplotlib.pyplot as plt
import numpy as np

DIR = '/mn/stornext/d5/data/duncanwa/WMAP'
cg.plot(f'{DIR}/chains_grid_g0_1.177/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel=r'g_0=1.178', cbar=False, width=4)
plt.savefig('ame_g01.177.pdf', bbox_inches='tight', dpi=300)
cg.plot(f'{DIR}/chains_grid_g0_1.181/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel='g_0=1.181', cbar=False, width=4)
plt.savefig('ame_g01.181.pdf', bbox_inches='tight', dpi=300)
cg.plot(f'{DIR}/chains_grid_g0_1.185/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel='g_0=1.186', cbar=False, width=4)
plt.savefig('ame_g01.186.pdf', bbox_inches='tight', dpi=300)

cg.standalone_colorbar("coolwarm", ticks=[-100,0,100], extend='both',
    unit=r'$\mathrm{\mu K}$ @ 22\,GHz')
plt.savefig('cbar_100uK.pdf', bbox_inches='tight', dpi=300)
