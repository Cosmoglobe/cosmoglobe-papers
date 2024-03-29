import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
n_dof = 300

f = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_b_230203/chisq_c0001_k000100.fits'

chisq = hp.read_map(f, field=(0,1,2))
chisq_IQU = chisq.sum(axis=0)
chisq_red = (chisq_IQU - 3*n_dof)/(2*3*n_dof)**0.5
cg.plot(chisq_red, min=-3, max=3, cmap='RdBu_r', cbar=False, width=4)
plt.savefig('chisq_IQU.pdf', bbox_inches='tight')

chisq_red = (chisq - n_dof)/(2*n_dof)**0.5


cg.plot(chisq_red, width=4, min=-3, max=3, cbar=False, cmap='RdBu_r', rlabel='I')
plt.savefig('chisq_I.pdf', bbox_inches='tight')
cg.plot(chisq_red, sig=1, width=4, min=-3, max=3, cbar=False, cmap='RdBu_r', rlabel='Q')
plt.savefig('chisq_Q.pdf', bbox_inches='tight')
cg.plot(chisq_red, sig=1, width=4, min=-3, max=3, cbar=False, cmap='RdBu_r', rlabel='U')
plt.savefig('chisq_U.pdf', bbox_inches='tight')



cg.standalone_colorbar("RdBu_r", ticks=[-3,0,3], extend='both',
        unit=r"$\chi^2\ [\sigma]$",width=3)

plt.savefig('cbar_3sigma.pdf', bbox_inches='tight')

