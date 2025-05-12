import healpy as hp
import h5py
import numpy as np
import cosmoglobe as cg

import matplotlib.pyplot as plt


out_map = np.zeros(512*512*12)

#fir catalogue

fir_cats = ['/mn/stornext/d23/cmbco/cg/dirbe/data/crossmatch_gaia_allwise_feb25_mag8_5arcsec_unique_7magcutoff_v3_nonparametrized_template_1.dat', '/mn/stornext/d23/cmbco/cg/dirbe/data/crossmatch_gaia_allwise_feb25_mag8_5arcsec_unique_7magcutoff_v3_nonparametrized_template_2.dat']

for fir_cat in fir_cats:
    coords = np.loadtxt(fir_cat, usecols=(0,1))

    glon = coords[:,0]
    glat = coords[:,1] 

    print(np.shape(glat))

    pix = hp.ang2pix(512, glon, glat, lonlat=True)

    for p in pix:
        out_map[p] +=1

#gaia catalogue
gaia_cats = ['/mn/stornext/d23/cmbco/cg/dirbe/data/commander_star_model_v3_feb25_mag8_5arcsec_7magcutoff_parametrized_template_1.h5', '/mn/stornext/d23/cmbco/cg/dirbe/data/commander_star_model_v3_feb25_mag8_5arcsec_7magcutoff_parametrized_template_2.h5']

for gaia_cat in gaia_cats:
    cat = h5py.File(gaia_cat)

    coords = cat['coordinates']

    print(max(coords[:,0]), min(coords[:,0]), max(coords[:,1]), min(coords[:,1]))
    print(np.shape(coords), np.shape(coords[:,0]))
    pix = hp.ang2pix(512, np.degrees(coords[:,0]), np.degrees(coords[:,1]), lonlat=True)

    for p in pix:
        out_map[p] += 1

print(max(out_map))

plt.figure()

cg.plot(out_map, min=0, max=4, cmap='gist_rainbow', rlabel='\mathrm{Star \ Count}')

plt.savefig('source_count.pdf')
