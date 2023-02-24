import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft
import h5py 

n_data = 7


bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
bands2 = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

n_pids = [1066, 1069, 1330, 1328, 897, 895, 1346, 1327, 1329, 1328]

Nobs_list = [12, 12, 15, 15, 20, 20, 30, 30, 30, 30]

freqs = [23, 30, 40, 40, 60, 60, 90, 90, 90, 90]
n_dets = 4



data_sampavg = {}

data_pidavg = {}

n_chains = 2

chain_id = ['a_230206', 'b_230203']

n_samp = [108, 160]  # extend this if the chains get longer
n_samples = n_samp[1]  # the longest
for n_pid, freq, band, band2 in zip(n_pids, freqs, bands, bands2):
    print(band)
    data_band = np.zeros((n_chains, n_samples, n_dets, n_data, n_pid))
    data_band[:] = np.nan
    for j in range(n_chains):
        #filename = '/mn/stornext/d16/cmbco/bp/delivery/v8.00/chains/chains_BP8_c51/chain_c000%i.h5' % (j+1)
        #filename = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/v2/BP_c000%i_v2.h5' % (j+1)
        # filename = '/mn/stornext/d16/cmbco/bp/delivery/v10.00/chains/chains_BP10_c%i/chain_c0001.h5' % (j + 22)
        filename = '/mn/stornext/d5/data/duncanwa/WMAP/chains_CG_' + chain_id[j] + '/chain_c0001.h5'
        print(filename)
        with h5py.File(filename, mode="r") as my_file:
            for i in range(n_samp[j]):
                mjd_array = np.zeros((n_dets, n_pid))
                istr = "{0:0=3d}".format(i)
                # print(istr)
                folder = '/000' + str(istr) + '/tod/%03i-WMAP_' % freq + band +'/'
                mjd_array[:, :] = np.array(my_file[folder + 'MJD'])[None, :]
                data_band[j, i, :, 0, :] = my_file[folder + 'accept']
                wh = np.where(data_band[j, i, :, 0, :] == 1.0)
                data_band[j, i, :, 1, :][wh] = np.array(my_file[folder + 'xi_n'])[2][wh]  # alpha
                data_band[j, i, :, 2, :][wh] = np.array(my_file[folder + 'chisq'])[wh]
                data_band[j, i, :, 3, :][wh] = np.array(my_file[folder + 'xi_n'])[1][wh]  # fknee
                data_band[j, i, :, 4, :][wh] = np.array(my_file[folder + 'gain'])[wh]
                data_band[j, i, :, 5, :][wh] = np.array(my_file[folder + 'xi_n'])[0][wh]  # sigma
                data_band[j, i, :, 6, :][wh] = mjd_array[wh]
    burnin = 10

    data_sampavg[band2] = np.nanmean(data_band[:,burnin:], (0, 1))
    data_pidavg[band2] = np.nanmean(data_band, 4)




np.save('noise_data_sampavg_WMAP_2chain.npy', data_sampavg)

np.save('noise_data_pidavg_WMAP_2chain.npy', data_pidavg)
