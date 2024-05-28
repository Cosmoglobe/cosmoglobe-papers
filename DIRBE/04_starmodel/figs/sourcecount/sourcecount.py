import healpy as hp
import h5py
import numpy as np
import cosmoglobe as cg

import matplotlib.pyplot as plt


out_map = np.zeros(512*512*12)

#fir catalogue

fir_cat = '/mn/stornext/d5/data/metins/dirbe/data/bsmith_table1_commander_v4.dat'

coords = np.loadtxt(fir_cat, usecols=(0,1))

glon = coords[:,0]
glat = coords[:,1] 

print(np.shape(glat))

pix = hp.ang2pix(512, glon, glat, lonlat=True)

for p in pix:
    out_map[p] +=1

#gaia catalogue

gaia_cat = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'

cat = h5py.File(gaia_cat)

coords = cat['coordinates']

print(max(coords[:,0]), min(coords[:,0]), max(coords[:,1]), min(coords[:,1]))
print(np.shape(coords), np.shape(coords[:,0]))
pix = hp.ang2pix(512, np.degrees(coords[:,0]), np.degrees(coords[:,1]), lonlat=True)

for p in pix:
    out_map[p] += 1

print(max(out_map))

plt.figure()

cg.plot(out_map, min=0, max=6, cmap='gist_rainbow', rlabel='\mathrm{Star \ Count}')

plt.savefig('source_count.pdf')
