import matplotlib.pyplot as plt

import numpy as np
import cosmoglobe as cg

import astropy.units as u
import astropy.constants as const

dirbe_bands = [1.25, 2.2, 3.5, 4.9, 12, 25] #microns

chainsdir = '/mn/stornext/d23/cmbco/dirbe/DR2_v2.00/'

chains = ['chains_prod4_c1', 'chains_prod4_c2']#, 'chains_v1.05_c03', 'chains_v1.05_c04', 'chains_v1.05_c05', 'chains_v1.05_c06',]

amps = []
theta1 = []
theta2 = []

kb = 1.3806503e-23
h = 1.0545726691251021e-34 * 2 *np.pi 
nu_ref = 239833.996
c = 299792458

GHz_to_um = 299792.458

def evalStar(nu, theta1, theta2):
    x = h/(kb*theta2)
    power_uk = ((np.exp(x*nu_ref) - 1)/(np.exp(x*nu) -1) * np.power(nu/nu_ref, theta1+1))

    C = (2*kb*(1e9*nu)**2/c**2)

    #print(C)

    return power_uk * C


fir_cat1 = '/mn/stornext/d23/cmbco/dirbe/data/crossmatch_gaia_allwise_feb25_mag8_5arcsec_unique_7magcutoff_v3_nonparametrized_template_1.dat'
fir_cat2 = '/mn/stornext/d23/cmbco/dirbe/data/crossmatch_gaia_allwise_feb25_mag8_5arcsec_unique_7magcutoff_v3_nonparametrized_template_2.dat'

coords1 = np.loadtxt(fir_cat1, usecols=(0,1))
coords2 = np.loadtxt(fir_cat2, usecols=(0,1))

glon = np.append(coords1[:,0], coords2[:,0])
glat = np.append(coords1[:,1], coords2[:,1])


for chain in chains:
    print(chain)
    chain = cg.Chain(chainsdir+ chain + '/chain_c0001.h5', burn_in=10)

    if len(amps) == 0:

        amps_1 = chain.mean('stars_mbb1/amp')/len(chains)
        thetas_1 = chain.mean('stars_mbb1/specind')
        theta1_1 = thetas_1[0,0,:]/len(chains)
        theta2_1 = thetas_1[1,0,:]/len(chains)

        amps_2 = chain.mean('stars_mbb2/amp')/len(chains)
        thetas_2 = chain.mean('stars_mbb2/specind')
        theta1_2 = thetas_2[0,0,:]/len(chains)
        theta2_2 = thetas_2[1,0,:]/len(chains)

    
    else:
        amps_1 += chain.mean('stars_mbb1/amp')/len(chains)

        thetas_1 = chain.mean('stars_mbb1/specind')
        theta1_1 += thetas[0,0,:]/len(chains)
        theta2_1 += thetas[1,0,:]/len(chains)

        amps_2 += chain.mean('stars_mbb2/amp')/len(chains)

        thetas_2 = chain.mean('stars_mbb2/specind')
        theta1_2 += thetas[0,0,:]/len(chains)
        theta2_2 += thetas[1,0,:]/len(chains)


microns = np.logspace(0,1.5)

amps = np.append(amps_1, amps_2, axis=1)
theta1 = np.append(theta1_1, theta1_2)
theta2 = np.append(theta2_1, theta2_2)
 

plt.figure()

data = np.zeros((len(theta1), len(dirbe_bands)))


for i in range(0, len(theta1)):
    data[i] = evalStar(np.divide(GHz_to_um,dirbe_bands), theta1[i], theta2[i]) * amps[0,i]*1e6
    #if(data[i,0] < 0):
    #    print(i, data[i], theta1[i], theta2[i], amps[0,i])

ax = plt.gca()

print(np.shape(data), np.shape(data[:,0]))

for i, band in enumerate(dirbe_bands):
    print(i)
    bin_data = data[:,i]
    bin_data = bin_data[bin_data > 0]

    print(np.shape(bin_data))

    if(band == 1.25):
        inds = bin_data.argsort()    
        
        #for i in range(1, 11):
            #print(glon[inds][-i], glat[inds][-i], bin_data[inds][-i])
        

        for i in range(1,11):
            index = np.random.randint(0,len(bin_data))
            print(glon[inds][index], glat[inds][index], bin_data[inds][index])

    bins = np.logspace(np.log10(min(bin_data)), np.log10(max(bin_data)), num=500)

    hist, edges = np.histogram(bin_data, bins=bins, density=False)

    stats ={}
    stats['coords'] = 0.5*(edges[:-1] + edges[1:])
    stats['vals'] = hist
    stats['mean'] = np.mean(bin_data)

    stats['median'] = np.median(bin_data)
    stats['min'] = np.min(bin_data)
    stats['max'] = np.max(bin_data)


    label=''
    if(band == 1.25):
        label='Extra Source SEDs'

    plot= ax.violin([stats], positions=[band], side='high', widths=1.8, showextrema=False)

    for bd in plot['bodies']:
        bd.set_color('blue')

#median SED

amps_pos = amps[0][amps[0] > 0]

median_curve = np.median(amps_pos)*evalStar(np.divide(GHz_to_um,microns), np.median(theta1), np.median(theta2))*1e6

plt.plot(microns, median_curve, color='red', label='Median SED')


plt.yscale('log')
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Emission (MJy/sr)')
plt.legend(loc='best')
plt.xlim(0.01, 30)
#plt.ylim(1e-7, 10)

#plt.title('Extra Source SEDs')
#
plt.savefig('exgal_spectra.pdf')
    

