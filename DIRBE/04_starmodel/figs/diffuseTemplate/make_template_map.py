import cosmoglobe as cg
import healpy as hp

import matplotlib.pyplot as plt


template = '/mn/stornext/d16/cmbco/bp/mathew/commander3/wiseDiffuse/diffuse_star_template_DIRBE01_smoothe.fits'

map_in = hp.read_map(template)

cg.plot(map_in, cmap='Oranges', unit='mJy', rlabel='DIRBE 01', llabel='S_{diffuse}')

plt.savefig('diffuse_stars.pdf')
