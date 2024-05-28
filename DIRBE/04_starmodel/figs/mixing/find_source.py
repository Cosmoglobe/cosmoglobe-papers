import healpy
import numpy as np
import matplotlib.pyplot as plt
import h5py

file = '/mn/stornext/d16/cmbco/cg/dirbe/DR2_v0.23/chains_v0.23_c01/chain_c0001.h5'

chains = h5py.File(file)

stars = chains['/' + str(83).zfill(6) + '/gaia/amp'][0]

for i in range(len(stars)):
    if stars[i] > 1:
        print(i)
    

#ind = np.argpartition(stars, -3)[-3:]
#print(ind)
