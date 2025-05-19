import cosmoglobe as cg
import numpy as np
import h5py 

import matplotlib.pyplot as plt

mean_emission = [0.5228979990674149, 0.483617344561016, 0.21522293944861112, 0.12503070382008852, 0.02770944012167716, 0.01172312982552367]

dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25] #microns
scalings = [1.82047631, 1.14461531, 0.56037995, 0.4172607, 0.04969875, 0.02000054]


chains_dir = '/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/'
catalog1 = '/mn/stornext/d23/cmbco/dirbe/data/commander_star_model_v3_feb25_mag8_5arcsec_7magcutoff_parametrized_template_1.h5'
catalog2 = '/mn/stornext/d23/cmbco/dirbe/data/commander_star_model_v3_feb25_mag8_5arcsec_7magcutoff_parametrized_template_2.h5' 

n_samps = 200
n_burnin = 10

chains = ['chains_prod4_c1', 'chains_prod4_c2']#, 'chains_v1.02_c03', 'chains_v1.02_c04','chains_v1.02_c05','chains_v1.02_c06']

count = 0

data1 = np.zeros(135486)
data2 = np.zeros(46623)

for chain in chains:

    chainfile = h5py.File(chains_dir + chain + '/chain_c0001.h5')

    for i in range(n_burnin, n_samps):
        try:
            data1 += chainfile['/' + str(i).zfill(6) + '/gaia1/amp'][:][0]
            data2 += chainfile['/' + str(i).zfill(6) + '/gaia2/amp'][:][0]
            count+= 1
        except KeyError:
            break

print(count)
data1 = data1 / count
data2 = data2 / count

cat1 = h5py.File(catalog1)
cat2 = h5py.File(catalog2)

plt.figure()

#plt.plot(dirbe_bands, data[0]/cat['reported_values'][0][0]*cat['reported_values'][0][0:6], alpha=0.1, color='blue', label='Star SEDs')

#for i in range(1, len(data)):
#    if i % 1000 == 0:
#        print(i)
#    plt.plot(dirbe_bands, data[i]/cat['reported_values'][i][0]*cat['reported_values'][i][0:6], alpha=0.1, color='blue')

for i, band in enumerate(dirbe_bands):
    #bin each set of SEDs vertically in 10000 bins
   
   
    catdata1 = cat1['reported_values'][()]
    catdata2 = cat2['reported_values'][()]
 
    bin_data1 = data1[:]/catdata1[:,0]*catdata1[:,i]
    bin_data2 = data2[:]/catdata2[:,0]*catdata2[:,i]

    bin_data = np.append(bin_data1, bin_data2)

    bin_data[bin_data < 0] = np.nan

    print(np.shape(bin_data), bin_data[0], bin_data[3000])

    plot = plt.violinplot(bin_data, positions=[band], side='high', widths=1, showextrema=False)

    for bd in plot['bodies']:
        bd.set_color('blue')


#plt.plot(dirbe_bands, scalings, color='red', label='Full Sky Mean')


plt.yscale('log')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Amplitude (mJy/Sr)')
#plt.legend(loc='best')


plt.title('Star SEDs')
plt.savefig('star_seds.pdf')
