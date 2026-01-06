import cosmoglobe as cg
import numpy as np
import h5py 

import matplotlib.pyplot as plt

dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25] #microns
scalings = np.array([1.9566207, 1.31449094, 0.43703894, 0.24464909, 0.13159147, 0.07952643])
stddevs = np.array([0.00333873, 0.00777074, 0.00488747, 0.00427976, 0.00988609, 0.01207999])

chains_dir = '/mn/stornext/d23/cmbco/dirbe/DR2_v3.02/'

gaia_cats = ['/mn/stornext/d23/cmbco/dirbe/data/mar25_mag12_5arcsec_withmcmc_mag7cutoff_parametrized_template_1.h5', '/mn/stornext/d23/cmbco/dirbe/data/mar25_mag12_5arcsec_withmcmc_mag7cutoff_parametrized_template_2.h5', '/mn/stornext/d23/cmbco/dirbe/data/mar25_mag12_5arcsec_withmcmc_mag7cutoff_parametrized_template_3.h5']


n_samps = 485
n_burnin = 50

chains = ['chains_dirbe_prod6_c1', 'chains_dirbe_prod6_c2', 'chains_dirbe_prod6_c3', 'chains_dirbe_prod6_c4','chains_dirbe_prod6_c5']#,'chains_dirbe_prod6_c6']

count = 1

data1 = []
medians = [[],[],[],[],[],[]]

cat1 = h5py.File(gaia_cats[0])
cat2 = h5py.File(gaia_cats[1])
cat3 = h5py.File(gaia_cats[2])


catdata1 = cat1['reported_values'][()]
catdata2 = cat2['reported_values'][()]
catdata3 = cat3['reported_values'][()]


for chain in chains:

    chainfile = h5py.File(chains_dir + chain + '/chain_c0001.h5')

    for i in range(n_burnin, n_samps, 2):
        if(len(data1) == 0) :
            data1 = chainfile['/' + str(i).zfill(6) + '/gaia1/amp'][:][0]
            data2 = chainfile['/' + str(i).zfill(6) + '/gaia2/amp'][:][0]
            data3 = chainfile['/' + str(i).zfill(6) + '/gaia3/amp'][:][0]
            for j,band in enumerate(dirbe_bands):
                new_data1 = data1[:]/catdata1[:,0]*catdata1[:,j]
                new_data2 = data2[:]/catdata2[:,0]*catdata2[:,j]
                new_data3 = data3[:]/catdata3[:,0]*catdata3[:,j]

                new_data = np.append(new_data1, new_data2)
                new_data = np.append(new_data, new_data3)

                new_data = new_data[new_data > 0]

                medians[j].append(np.median(new_data))

        else:

            try:
                data1 += chainfile['/' + str(i).zfill(6) + '/gaia1/amp'][:][0]
                data2 += chainfile['/' + str(i).zfill(6) + '/gaia2/amp'][:][0]
                data3 += chainfile['/' + str(i).zfill(6) + '/gaia3/amp'][:][0]
                
                for j,band in enumerate(dirbe_bands):
                    new_data1 = chainfile['/' + str(i).zfill(6) + '/gaia1/amp'][:][0]/catdata1[:,0]*catdata1[:,j]
                    new_data2 = chainfile['/' + str(i).zfill(6) + '/gaia2/amp'][:][0]/catdata2[:,0]*catdata2[:,j]
                    new_data3 = chainfile['/' + str(i).zfill(6) + '/gaia3/amp'][:][0]/catdata3[:,0]*catdata3[:,j]

                    new_data = np.append(new_data1, new_data2)
                    new_data = np.append(new_data, new_data3)
                    new_data = new_data[new_data > 0]

                    medians[j].append(np.median(new_data))

                count+= 1
            except KeyError:
                break

print(count)
data1 = data1 / count
data2 = data2 / count
data3 = data3 / count

cat1 = h5py.File(gaia_cats[0])
cat2 = h5py.File(gaia_cats[1])
cat3 = h5py.File(gaia_cats[2])


#plt.plot(dirbe_bands, data[0]/cat['reported_values'][0][0]*cat['reported_values'][0][0:6], alpha=0.1, color='blue', label='Star SEDs')

#for i in range(1, len(data)):
#    if i % 1000 == 0:
#        print(i)
#    plt.plot(dirbe_bands, data[i]/cat['reported_values'][i][0]*cat['reported_values'][i][0:6], alpha=0.1, color='blue')

means = np.zeros(6)

plt.figure()
for i, band in enumerate(dirbe_bands):

    bin_data1 = data1[:]/catdata1[:,0]*catdata1[:,i]
    bin_data2 = data2[:]/catdata2[:,0]*catdata2[:,i]
    bin_data3 = data3[:]/catdata3[:,0]*catdata3[:,i]

    bin_data = np.append(bin_data1, bin_data2)
    bin_data = np.append(bin_data, bin_data3)

    bins = np.logspace(1e-12, 10, num=50)
    plt.hist(bin_data, label=str(band), log=True, bins=bins)

plt.xscale('log')
plt.legend(loc='best')
plt.savefig('testhists.pdf')


plt.figure()
for i, band in enumerate(dirbe_bands):
    #bin each set of SEDs vertically in 10000 bins
    
    bin_data1 = data1[:]/catdata1[:,0]*catdata1[:,i]
    bin_data2 = data2[:]/catdata2[:,0]*catdata2[:,i]
    bin_data3 = data3[:]/catdata3[:,0]*catdata3[:,i]

    bin_data = np.append(bin_data1, bin_data2)
    bin_data = np.append(bin_data, bin_data3)

    bin_data = bin_data[bin_data > 0]

    print(np.shape(bin_data), bin_data[0], bin_data[3000])
    print(np.mean(bin_data), np.max(bin_data), np.min(bin_data))

    #plot = plt.violinplot(bin_data, positions=[band], side='high', widths=1, showextrema=False)

    means[i] = np.median(bin_data)

    #for bd in plot['bodies']:
    #    bd.set_color('blue')


#plt.plot(dirbe_bands, scalings, color='red', label='Full Sky Mean')

#plt.plot(dirbe_bands, means, color='orange', label='Median SED')
#plt.plot(dirbe_bands, scalings * (means[0] / scalings[0]), color='red', label='Diffuse SED')

print(np.std(medians, axis=1))

plt.errorbar(dirbe_bands, means, yerr=np.std(medians, axis=1), color='orange', label='Stars Median')

plt.errorbar(dirbe_bands, scalings * (means[0] / scalings[0]), yerr=stddevs*(means[0] / scalings[0]), color='red', label='Diffuse SED')

plt.yscale('log')
#plt.ylim(1e-12, 10)

#plt.xscale('log')

plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Amplitude (MJy/sr)')
plt.legend(loc='best')

#plt.title('Star SEDs')
plt.savefig('star_seds.pdf')
