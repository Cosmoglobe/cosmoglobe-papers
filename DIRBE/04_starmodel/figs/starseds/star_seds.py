import cosmoglobe as cg
import h5py 

import matplotlib.pyplot as plt

mean_emission = [0.5228979990674149, 0.483617344561016, 0.21522293944861112, 0.12503070382008852, 0.02770944012167716, 0.01172312982552367]

dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25] #microns

chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.10/chains_v0.10'
catalog = '/mn/stornext/d5/data/metins/dirbe/data/commander_star_model.h5'


n_samps = 22
n_burnin = 10

chainfile = h5py.File(chains_dir + '/chain_c0001.h5')

data = chainfile['/' + str(n_burnin).zfill(6) + '/gaia/amp'][:]

for i in range(n_burnin, n_samps):
    data += chainfile['/' + str(i).zfill(6) + '/gaia/amp'][:]

data = data / (n_samps - n_burnin)
data = data[0]

cat = h5py.File(catalog)

plt.figure()

plt.plot(dirbe_bands, data[0]/cat['reported_values'][0][0]*cat['reported_values'][0][0:6], alpha=0.1, color='blue', label='Star SEDs')

for i in range(1, 10000): #len(data)):
    if i % 1000 == 0:
        print(i)
    plt.plot(dirbe_bands, data[i]/cat['reported_values'][i][0]*cat['reported_values'][i][0:6], alpha=0.1, color='blue')

plt.plot(dirbe_bands, mean_emission, color='red', label='Full Sky Mean')

plt.yscale('log')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Amplitude (mJy)')
plt.legend(loc='best')


plt.title('Star SEDs')
plt.savefig('star_seds.pdf')
