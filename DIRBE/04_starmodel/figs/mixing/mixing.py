import cosmoglobe as cg
import numpy as np
import h5py
import matplotlib.pyplot as plt

chains_dir = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.15/'

chain1 = 'chains_v0.22_c01'
chain2 = 'chains_v0.22_c02'
chain3 = 'chains_v0.22_c03'
chain4 = 'chains_v0.22_c04'
chain5 = 'chains_v0.22_c05'
chain6 = 'chains_v0.22_c06'


num_samps = 43
burn_in = 20

plt.figure()

for chain in [chain1, chain2, chain3, chain4, chain5. chain6]:
    chainfile = h5py.File(chains_dir + chain + '/chain_c0001.h5')
    data = []
    
    eyes = []

    for i in range(burn_in, num_samps, 2):
        amp = chainfile['/' + str(i).zfill(6) + '/stars2/amp_alm'][0][0]

        print(i, amp)

        data.append(amp)
        eyes.append(i)

    plt.plot(eyes, data, label=chain)


plt.savefig('template_mixing.pdf')
