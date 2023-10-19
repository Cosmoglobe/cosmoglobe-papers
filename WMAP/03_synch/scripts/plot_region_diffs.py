import cosmoglobe as cg
import matplotlib.pyplot as plt
import numpy as np

norm_dict = {'badcolor':'0.3'}
cg.plot('CG_023-WMAP_K_diff_wmap9_v1_t1.5_maskreg_all.fits', sig=1, min=-5e-3,
        max=5e-3, norm_dict=norm_dict, sub=(2,2,1), cbar=False, rlabel='Q',
        llabel='K')
cg.plot('CG_023-WMAP_K_diff_wmap9_v1_t1.5_maskreg_all.fits', sig=2, min=-5e-3,
        max=5e-3, norm_dict=norm_dict, sub=(2,2,2), cbar=False, rlabel='U')
cg.plot('CG_030_diff_dx12_v1_t1.5_maskreg_all.fits', sig=1, min=-5, max=5,
        norm_dict=norm_dict, sub=(2,2,3), cbar=False, llabel='30')
cg.plot('CG_030_diff_dx12_v1_t1.5_maskreg_all.fits', sig=2, min=-5, max=5,
        norm_dict=norm_dict, sub=(2,2,4), cbar=False)
plt.subplots_adjust(wspace=0.05, hspace=-0.15)

plt.savefig('../figures/CG_diff_regions.pdf', bbox_inches='tight', dpi=150)

