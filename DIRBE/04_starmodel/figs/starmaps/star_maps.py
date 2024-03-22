import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.10/chains_v0.10'
catalog = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'


n_samps = 22
n_burnin = 10

bands = ['01', '02', '03', '04', '05', '06']


nmaps = n_samps - n_burnin

for band in bands:
    
    map_mean = hp.read_map(chains_dir + '/gaia_' + band + 'a_c0001_k'+str(n_burnin).zfill(6) + '.fits')

    map_square = np.power(map_mean, 2)/nmaps

    map_mean = map_mean/ nmaps

    #one pass mean and stddev

    for i in range(n_burnin+1, n_samps):
        map_in = hp.read_map(chains_dir + '/gaia_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

        map_mean += map_in/nmaps

        map_square += np.power(map_in,  2) /nmaps
        
    map_square -= np.power(map_mean, 2)
    map_square = np.sqrt(map_square)

    rlabel = 'DIRBE' + band

    plt.figure()
    map_mean[map_mean<0] = 0
    print(min(map_mean), max(map_mean), np.percentile(map_mean, 97.5))
    cg.plot(map_mean, min=0, max=2, cmap='Oranges', rlabel=rlabel, llabel='S_{stars}', unit='mJy')

    plt.savefig('stars_mean_'+band+'.pdf')
    plt.close()

    plt.figure()
    cg.plot(map_mean, min=0, max=1, cmap='grey', rlabel=rlabel, llabel='RMS_{stars}', unit='mJy')
    plt.savefig('stars_std_'+band+'.pdf')
    plt.close()

 
