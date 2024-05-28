import healpy
import numpy as np
import matplotlib.pyplot as plt
import h5py

dirs = ['/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.23/chains_v0.23_c01/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.23/chains_v0.23_c02/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.23/chains_v0.23_c03/',
        '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.23/chains_v0.23_c04/']

colors = ['blue', 'red', 'green', 'orange']


rand_stars = [191352, 199550, 212190] #rand-ish bright

#[305694, 560770, 528111] brightest
#[536892, 17000, 123456] #rand
rand_fir = [5639, 45123, 60488]

for sdir, colour in zip(dirs, colors):

    chains = h5py.File(sdir + 'chain_c0001.h5')

    samples = []
    monos = []
    albedos = []
    stars = [[],[],[]]
    firs = [[],[],[]]


    first_sample = np.mean(healpy.read_map(sdir + 'stars2_c0001_k000000.fits'))
    curr_ind = 1
    while True:
        try:
            curr_sample = np.mean(healpy.read_map(sdir+f'stars2_c0001_k{curr_ind:06d}.fits'))
            print('opened ' + sdir+f'stars2_c0001_k{curr_ind:06d}.fits')
            samples.append(curr_sample/first_sample)
        except FileNotFoundError as e:
            break

        mono = np.mean(np.genfromtxt(sdir + 'md_c0001_k'+ str(curr_ind).zfill(6)+ '.dat', usecols=(1))[0:2])
        monos.append(mono)

        albedo = np.mean(chains['/' + str(curr_ind).zfill(6) + '/zodi/comps/cloud/albedo'][0:2])
        albedos.append(albedo)


        for i in range(0,3):
            star = chains['/' + str(curr_ind).zfill(6) + '/gaia/amp'][0][rand_stars[i]]
            stars[i].append(star)

            fir = chains['/' + str(curr_ind).zfill(6) + '/stars/amp'][0][rand_fir[i]]
            firs[i].append(fir)

        curr_ind+=2

    chains.close()

    plt.plot(stars[0], color=colour)
    plt.plot(stars[1], color=colour)
    plt.plot(stars[2], color=colour)



plt.xlabel('Sample')
plt.ylabel('Gaia Amplitude')
plt.savefig('gaia_trace.pdf', bbox_inches='tight')
