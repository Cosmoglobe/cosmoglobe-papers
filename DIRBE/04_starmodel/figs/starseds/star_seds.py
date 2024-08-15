import cosmoglobe as cg
import numpy as np
import h5py 

import matplotlib.pyplot as plt

mean_emission = [0.5228979990674149, 0.483617344561016, 0.21522293944861112, 0.12503070382008852, 0.02770944012167716, 0.01172312982552367]

dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25] #microns
scalings = [1, 7.133751886520310892e-01, 3.387313865841199978e-01, 1.778017878718340106e-01, 3.973778737936959488e-02, 1.574453456726415679e-02]


chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v1.02/'
catalog = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'


n_samps = 200
n_burnin = 49

chains = ['chains_v1.02_c01', 'chains_v1.02_c02', 'chains_v1.02_c03', 'chains_v1.02_c04','chains_v1.02_c05','chains_v1.02_c06']

count = 0

data = np.zeros(717454)

for chain in chains:

    chainfile = h5py.File(chains_dir + chain + '/chain_c0001.h5')

    for i in range(n_burnin, n_samps):
        try:
            data += chainfile['/' + str(i).zfill(6) + '/gaia/amp'][:][0]
            count+= 1
        except KeyError:
            break

print(count)
data = data / count

cat = h5py.File(catalog)

plt.figure()

#plt.plot(dirbe_bands, data[0]/cat['reported_values'][0][0]*cat['reported_values'][0][0:6], alpha=0.1, color='blue', label='Star SEDs')

#for i in range(1, len(data)):
#    if i % 1000 == 0:
#        print(i)
#    plt.plot(dirbe_bands, data[i]/cat['reported_values'][i][0]*cat['reported_values'][i][0:6], alpha=0.1, color='blue')

for i, band in enumerate(dirbe_bands):
    #bin each set of SEDs vertically in 10000 bins
   
   
    catdata = cat['reported_values'][()]

    print(np.shape(data), np.shape(catdata), np.shape(catdata[:,0]))
 
    bin_data = data[:]/catdata[:,0]*catdata[:,i]

    print(np.shape(bin_data), bin_data[0], bin_data[3000])

    plt.violinplot(bin_data, positions=[band], points=500, side='high', widths=1)


#plt.plot(dirbe_bands, scalings, color='red', label='Full Sky Mean')

#plt.yscale('log')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Amplitude (mJy/Sr)')
#plt.legend(loc='best')


plt.title('Star SEDs')
plt.savefig('star_seds.pdf')
