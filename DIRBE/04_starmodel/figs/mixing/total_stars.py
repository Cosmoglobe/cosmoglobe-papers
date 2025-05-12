import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import h5py

dirs = ['/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/chains_prod4_c1/',
        '/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/chains_prod4_c2/']
#        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c03/',
#        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/chains_v1.00_c04/']


pix = [1352704, 2429499, 3121308] #(0,8), LMC (280.46, -32.89), (0, -80)
locs = ['(0,8)', '(280.46, -32.89)', '(0, -80)']

print(pix, 3145728/2)

colors = ['blue', 'red']#, 'green', 'orange']

emissions = np.empty((len(colors), len(pix), 1000))
star_arr = np.empty((len(colors), len(pix), 1000))
exgal = np.empty((len(colors), len(pix), 1000))
diffuse = np.empty((len(colors), len(pix), 1000))

for d, (sdir, colour) in enumerate(zip(dirs, colors)):

    curr_ind = 1
    while True:
        i = int((curr_ind-1)/2)
        try:
            stars = hp.read_map(sdir + 'stars_mbb1_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
            stars += hp.read_map(sdir + 'stars_mbb2_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
            stars2 = hp.ud_grade(hp.read_map(sdir + 'stars_diff_c0001_k' + str(curr_ind).zfill(6) + '.fits')*1.82047631, 512)
            gaia = hp.read_map(sdir + 'gaia1_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
            gaia += hp.read_map(sdir + 'gaia2_01a_c0001_k' + str(curr_ind).zfill(6) + '.fits')
        except FileNotFoundError as e:
            print('File not found', curr_ind)
            break

        for p, px in enumerate(pix):

            tot_emission = stars[px] + gaia[px] + stars2[px]

            emissions[d,p,i] = tot_emission
            star_arr[d,p,i] = gaia[px]
            exgal[d,p,i] = stars[px]
            diffuse[d,p,i] = stars2[px]

        print(d, p, i, curr_ind)
        curr_ind +=2


for p, (location, px) in enumerate(zip(locs, pix)):
    plt.figure()
    fig, ax = plt.subplots(4, 1, sharex=True, figsize=(3, 10))

    fig.suptitle('(glon, glat)=' + location, y=0.9)

    for d, (sdir, colour) in enumerate(zip(dirs, colors)):

        ax[0].plot(emissions[d,p,:], color=colour)
        ax[1].plot(star_arr[d,p,:], color=colour)
        ax[2].plot(exgal[d,p,:], color=colour)
        ax[3].plot(diffuse[d,p,:], color=colour)

    #ax[0].axvline(10, linestyle='dashed', color='black')
    #ax[1].axvline(10, linestyle='dashed', color='black')
    #ax[2].axvline(10, linestyle='dashed', color='black')
    #ax[3].axvline(10, linestyle='dashed', color='black')

    if(p == 0):
        ax[0].set_ylabel('Combined Emission (MJy/sr)')
        ax[1].set_ylabel('Gaia Stars (MJy/sr)')
        ax[2].set_ylabel('MBB Sources (MJy/sr)')
        ax[3].set_ylabel('Diffuse Sources (MJy/sr)')

    plt.xlim((0,100))
    plt.xlabel('Sample')
    plt.savefig('trace_' + str(px)  + '.pdf', bbox_inches='tight')    
