import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.00/'
catalog = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'


n_samps = 16
n_burnin = 0

bands = ['01', '02', '03', '04', '05', '06']
scalings = [1, 7.133751886520310892e-01, 3.387313865841199978e-01, 1.778017878718340106e-01, 3.973778737936959488e-02, 1.574453456726415679e-02]
chains = ['chains_v1.00_c01', 'chains_v1.00_c02', 'chains_v1.00_c03', 'chains_v1.00_c04','chains_v1.00_c05','chains_v1.00_c06']

nmaps = int((n_samps - n_burnin)/2)


for band, scaling in zip(bands, scalings):
    for chain in chains:
    
        map_mean = hp.read_map(chains_dir + chain + '/gaia_' + band + 'a_c0001_k'+str(n_burnin).zfill(6) + '.fits')

        map_mean += hp.read_map(chains_dir + chain + '/stars_' + band + 'a_c0001_k' + str(n_burnin).zfill(6) + '.fits')

        map_mean += hp.read_map(chains_dir + chain + '/stars2_c0001_k' + str(n_burnin).zfill(6) + '.fits') * scaling

        map_square = np.power(map_mean, 2)/nmaps

        map_mean = map_mean/ nmaps

        #one pass mean and stddev

        for i in range(n_burnin+1, n_samps, 2):
            map_in = hp.read_map(chains_dir + chain +'/gaia_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            map_mean += map_in/nmaps

            ptsrc_in = hp.read_map(chains_dir + chain +'/stars_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            map_mean += ptsrc_in/nmaps

            diffuse_in = hp.read_map(chains_dir + chain +'/stars2_c0001_k'+str(i).zfill(6) + '.fits') * scaling

            map_mean += diffuse_in/nmaps

            map_square += np.power(map_in,  2) /nmaps

        map_square -= np.power(map_mean, 2)
        map_square = np.sqrt(map_square)

        rlabel = 'DIRBE' + band

        plt.figure()
        map_mean[map_mean<0] = 0
        print(np.sum(map_mean/np.size(map_mean)))
        cg.plot(map_mean, min=0, max=1, cmap='Oranges', rlabel=rlabel, llabel='S_{stars}', unit='MJy/sr')

        plt.savefig('all_stars_mean_'+band+'.pdf')
        plt.close()

        plt.figure()
        cg.plot(map_square, min=0, max=0.000001, cmap='grey', rlabel=rlabel, llabel='RMS_{stars}', unit='MJy/sr')
        plt.savefig('all_stars_std_'+band+'.pdf')
        plt.close()

