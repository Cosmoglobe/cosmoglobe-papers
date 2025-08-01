import cosmoglobe as cg
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp


chains_dir = '/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/'

n_samps = 206
n_burnin = 6

bands = ['01', '02', '03', '04', '05', '06']
#output by ../diffuseTemplate/make_template_map.py
scalings = [1.82047631, 1.14461531, 0.56037995, 0.4172607, 0.04969875, 0.02000054]
chains = ['chains_prod4_c1', 'chains_prod4_c2']#, 'chains_v1.00_c03', 'chains_v1.00_c04','chains_v1.00_c05','chains_v1.00_c06']

nmaps = int((n_samps - n_burnin)/2)


for scaling, band in zip(scalings, bands):
    for chain in chains:
    
        map_mean = hp.read_map(chains_dir + chain + '/gaia1_' + band + 'a_c0001_k'+str(n_burnin).zfill(6) + '.fits')

        map_mean += hp.read_map(chains_dir + chain + '/gaia2_' + band + 'a_c0001_k'+str(n_burnin).zfill(6) + '.fits')

        map_gaia = map_mean/nmaps

        map_mean += hp.read_map(chains_dir + chain + '/stars_mbb1_' + band + 'a_c0001_k' + str(n_burnin).zfill(6) + '.fits')

        map_mean += hp.read_map(chains_dir + chain + '/stars_mbb2_' + band + 'a_c0001_k' + str(n_burnin).zfill(6) + '.fits')

        map_extra = (map_mean - map_gaia)/nmaps

        map_mean += hp.ud_grade(hp.read_map(chains_dir + chain + '/stars_diff_c0001_k' + str(n_burnin).zfill(6) + '.fits') * scaling, 512)

        map_diffuse = hp.ud_grade(hp.read_map(chains_dir + chain + '/stars_diff_c0001_k' + str(n_burnin).zfill(6) + '.fits') * scaling, 512)/nmaps

        map_square = np.power(map_mean, 2)/nmaps

        map_mean = map_mean/ nmaps

        #one pass mean and stddev

        for i in range(n_burnin+1, n_samps, 2):
            print('sample', i)
            map_in = hp.read_map(chains_dir + chain +'/gaia1_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')
            map_in += hp.read_map(chains_dir + chain +'/gaia2_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            map_mean += map_in/nmaps
            map_gaia += map_in/nmaps

            ptsrc_in = hp.read_map(chains_dir + chain +'/stars_mbb1_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')
            ptsrc_in += hp.read_map(chains_dir + chain +'/stars_mbb2_' + band + 'a_c0001_k'+str(i).zfill(6) + '.fits')

            map_mean += ptsrc_in/nmaps
            map_extra += ptsrc_in/nmaps

            diffuse_in = hp.ud_grade(hp.read_map(chains_dir + chain +'/stars_diff_c0001_k'+str(i).zfill(6) + '.fits') * scaling, 512)

            map_mean += diffuse_in/nmaps
            map_diffuse += diffuse_in/nmaps
    
            map_square += np.power(map_in,  2) /nmaps

        map_square -= np.power(map_mean, 2)
        map_square = np.sqrt(map_square)

        rlabel = 'DIRBE' + band

        plt.figure()
        map_mean[map_mean<0] = 0
        print(np.sum(map_mean/np.size(map_mean)))
        cg.plot(map_mean, min=0, max=1, cmap='Oranges', rlabel=rlabel, llabel='S_{stars}', unit='MJy/sr')

        plt.savefig('all_stars_mean_'+band+'.pdf', bbox_inches='tight')
        plt.close()

        plt.figure()
        cg.plot(map_square, min=0, max=0.000001, cmap='grey', rlabel=rlabel, llabel='RMS_{stars}', unit='MJy/sr')
        plt.savefig('all_stars_std_'+band+'.pdf', bbox_inches='tight')
        plt.close()

        hp.write_map('all_stars_mean_' + band + '.fits', map_mean, overwrite=True)
        hp.write_map('all_stars_std_' + band + '.fits', map_square, overwrite=True)

        hp.write_map('stars_mean_' + band + '.fits', map_gaia, overwrite=True)
        hp.write_map('diffuse_mean_' + band + '.fits', map_diffuse, overwrite=True)
        hp.write_map('exgal_mean_' + band + '.fits', map_extra, overwrite=True)
