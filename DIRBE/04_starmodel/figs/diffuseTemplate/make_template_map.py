import cosmoglobe as cg
import healpy as hp
import numpy as np
import h5py
from astropy.io import fits
import quadcube

import matplotlib.pyplot as plt


template = '/mn/stornext/d23/cmbco/dirbe/data/diffuse_star_template_DIRBE01_combined_x2000.fits'

residuals = '/mn/stornext/d5/data/duncanwa/DIRBE/release_tests/DR2_v2/goodness/'

map_in = hp.read_map(template)

#map_in *= 2000

chains_dir = '/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/'

chains = ['chains_prod4_c1', 'chains_prod4_c2']#, 'chains_v1.00_c03', 'chains_v1.00_c04','chains_v1.00_c05','chains_v1.00_c06']

burn_in = 6
num_samps = 206

count = 0
scale = [0*10]

for chain in chains:
    cf = h5py.File(chains_dir + chain + '/chain_c0001.h5')
    for samp in range(burn_in, num_samps, 2):
        count +=1
        scale += cf['/' + str(samp).zfill(6) + '/stars_diff/SED'][()][:,2]
        
scale = scale/count

print(np.flip(scale))

map_in *= scale[9]

print(max(map_in))

cg.plot(map_in, cmap='Oranges', unit='MJy/sr', min=1, max=10000, rlabel='DIRBE 01', llabel='S_{diffuse}', norm='log')

plt.savefig('diffuse_stars_log.pdf', bbox_inches='tight')

cg.plot(map_in, cmap='Oranges', unit='MJy/sr', min=1, max=7, rlabel='DIRBE 01', llabel='S_{diffuse}')

plt.savefig('diffuse_stars.pdf', bbox_inches='tight')



#4096*pix9 + 85*16 + 5

dirbe_in = fits.open('DIRBE_BAND1A_FSM.FITS')

vecs = quadcube.pix2vec(4096*dirbe_in[1].data['Pixel_no'])#+85*16+5)
pix = hp.vec2pix(128, *vecs)

binned_map = np.zeros(12*128*128)
total_hits = np.zeros_like(binned_map)

bin_count = np.bincount(pix, weights=dirbe_in[1].data['Photomet'], minlength=len(binned_map))
binned_map[: len(bin_count)] += bin_count
unique_pix, counts = np.unique(pix, return_counts=True)
total_hits[unique_pix] += counts
non_zero_inds = total_hits > 0
binned_map[non_zero_inds] /= total_hits[non_zero_inds]

plt.figure()

map_128 = hp.ud_grade(map_in, 128)

r = hp.Rotator(coord=['E', 'G'])

binned_map = r.rotate_map_alms(binned_map)
binned_map[binned_map<0] = 0
cg.plot(binned_map, cmap='Oranges', unit='MJy/sr', min=0, max=2, rlabel='DIRBE 01', llabel='S_{fsm}')

diff = map_128 - binned_map

plt.savefig('dirbe_template.pdf', bbox_inches='tight')
plt.figure()

cg.plot(diff, cmap='Oranges', unit='MJy/sr', min=-1, max=1, rlabel='DIRBE 01', llabel='S_{diffuse} - S_{fsm}')

plt.savefig('diffuse_diff.pdf', bbox_inches='tight')

plt.figure()

res1 = hp.read_map(residuals + 'CG_res_01_I_n512_0arcmin_DR2_v2.fits')

cg.plot(res1, unit='MJy/sr', min=-0.1, max=0.1, rlabel='DIRBE 01', llabel='Res')
plt.savefig('band01_res.pdf', bbox_inches='tight')
