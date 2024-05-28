import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import h5py

dirs = ['/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c01/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c02/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c03/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c04/']


pix = hp.ang2pix(512, 0, -5, lonlat=True)

print(pix, 3145728/2)

rand_pixels = [hp.ang2pix(512, 0, 0, lonlat=True), hp.ang2pix(512, 0, -80, lonlat=True), hp.ang2pix(512, 0, 8, lonlat=True)]

colors = ['blue', 'red', 'green', 'orange']

for sdir, colour in zip(dirs, colors):

    emissions = [[],[],[]]

    curr_ind = 1
    while True:
        try:
            stars = hp.read_map(sdir + 'stars_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
            stars2 = hp.read_map(sdir + 'stars2_c0001_k' + str(curr_ind).zfill(6) + '.fits')
            gaia = stars = hp.read_map(sdir + 'gaia_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
        except FileNotFoundError as e:
            print('File not found', curr_ind)
            break

        for i in range(0,3):

            tot_emission = stars[rand_pixels[i]] + gaia[rand_pixels[i]] + stars2[rand_pixels[i]]

            #tot_emission = stars2[rand_pixels[i]]

            emissions[i].append(tot_emission)
        curr_ind +=2

        plt.plot(emissions[0], color=colour)
        plt.plot(emissions[1], color=colour)
        plt.plot(emissions[2], color=colour)

plt.xlabel('Sample')
plt.ylabel('Combined Star Emission')
plt.savefig('total_star_trace.pdf')    
