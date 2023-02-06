import cosmoglobe as cg
import matplotlib.pyplot as plt
import numpy as np

DIR = '/mn/stornext/d5/data/duncanwa/WMAP'
cg.plot(f'{DIR}/chains_grid_g0_1.178/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel=r'1.178', cbar=False, width=4)
plt.savefig('ame_g01.178.png', bbox_inches='tight', dpi=300)
cg.plot(f'{DIR}/chains_grid_g0_1.181/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel='1.181', cbar=False, width=4)
plt.savefig('ame_g01.181.png', bbox_inches='tight', dpi=300)
cg.plot(f'{DIR}/chains_grid_g0_1.185/ame_c0001_k000005.fits', min=-1e2, max=1e2, 
    cmap='coolwarm', llabel='1.185', cbar=False, width=4)
plt.savefig('ame_g01.185.png', bbox_inches='tight', dpi=300)

cg.standalone_colorbar("coolwarm", ticks=[-100,0,100], extend='both',
    unit=r'$\mathrm{\mu K}$ @ 22\,GHz')
plt.savefig('cbar_100uK.png', bbox_inches='tight', dpi=300)
