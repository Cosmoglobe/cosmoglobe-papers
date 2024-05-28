import healpy as hp

import matplotlib.pyplot as plt

map1 = hp.read_map('quadcube_centroid.fits')
map2 = hp.read_map('quadcube_corner.fits')

diff = map1-map2

hp.mollview(diff, min=-0.1, max=0.1)

plt.savefig('quadcube_diff.pdf')
