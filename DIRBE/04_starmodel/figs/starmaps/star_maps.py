import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.15/'
catalog = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'


n_samps = 43
n_burnin = 20

bands = ['01', '02', '03', '04', '05', '06']
chains = ['chains_v0.15_c01', 'chains_v0.15_c02']

nmaps = int((n_samps - n_burnin)/2)


for band in bands:
    for chain in chains:
    
        map_mean = hp.read_map(chains_dir + chain + '/gaia_' + band + 'a_c0001_k'+str(n_burnin).zfill(6) + '.fits')

        ptsrc_mean = hp.read_map(chains_dir + chain + '/stars_' + band + 'a_c0001_k' + str(n_burnin).zfill(6) + '.fits')

        map_square = np.power(map_mean, 2)/nmaps

        map_mean = map_mean/ nmaps

    
        ptsrc_square = np.power(ptsrc_mean, 2)/nmaps
        ptsrc_mean = ptsrc_mean/nmaps

        #one pass mean and stddev

        for i in range(n_burnin+1, n_samps, 2):
            map_in = hp.read_map(chains_dir + chain +'/gaia_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            map_mean += map_in/nmaps
            map_square += np.power(map_in,  2) /nmaps
        
            ptsrc_in = hp.read_map(chains_dir + chain +'/stars_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            ptsrc_mean += ptsrc_in/nmaps
            ptsrc_square += np.power(ptsrc_in, 2)/nmaps


        map_square -= np.power(map_mean, 2)
        map_square = np.sqrt(map_square)

        ptsrc_square -= np.power(ptsrc_mean, 2)
        ptsrc_square = np.sqrt(ptsrc_square)

        rlabel = 'DIRBE' + band

        plt.figure()
        map_mean[map_mean<0] = 0
        print(np.sum(map_mean/np.size(map_mean)))
        cg.plot(map_mean, min=0, max=2, cmap='Oranges', rlabel=rlabel, llabel='S_{stars}', unit='mJy')

        plt.savefig('stars_mean_'+band+'.pdf')
        plt.close()

        plt.figure()
        cg.plot(map_square, min=0, cmap='grey', rlabel=rlabel, llabel='RMS_{stars}', unit='mJy')
        plt.savefig('stars_std_'+band+'.pdf')
        plt.close()

        plt.figure()
        ptsrc_mean[ptsrc_mean<0] = 0
        cg.plot(ptsrc_mean, min=0, max=2, cmap='Reds', rlabel=rlabel, llabel='S_{ptsrc}', unit='mJy')
        
        plt.savefig('ptsrc_mean_'+band+'.pdf')

        plt.figure()
        cg.plot(ptsrc_square, min=0, cmap='grey', rlabel=rlabel, llabel='RMS_{ptsrc}', unit='mJy')
        plt.savefig('ptsrc_std_'+band+'.pdf')

 
