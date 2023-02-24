import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft
from scipy import interpolate
from matplotlib.pyplot import cm
import sys
import h5py

def one_over_f(freq, sigma0, fknee):
    return sigma0 ** 2 * (1.0 + fknee * freq)

def bin_data_maxbin(arr, slope=1.0, cut=20, maxbin=2000):
    n = len(arr)
    n_coadd = []
    cumsum = 0
    while np.sum(n_coadd) < n:
        n_add = min(1 + int((cumsum//cut) ** slope), maxbin)
        n_coadd.append(n_add)
        cumsum += n_add
    sum = np.sum(n_coadd[:-1])
    n_coadd[-2] += n - sum
    n_coadd = np.array(n_coadd[:-1]).astype(int)
    newarr = np.zeros(len(n_coadd))
    cumsum = 0
    for i in range(len(n_coadd)):
        newarr[i] = np.sum(arr[cumsum:cumsum+n_coadd[i]]) / n_coadd[i]
        cumsum += n_coadd[i]
    return newarr


# path = '/mn/stornext/d16/cmbco/bp/dwatts/WMAP/faz_chains_DAs/'
path = '/mn/stornext/d5/data/duncanwa/WMAP/chains_writeW4_230221/'
# path = '/mn/stornext/d5/data/duncanwa/WMAP/chains_writeW2_230217/'
filename = path + 'tod_000050_samp000002.h5'
gain = np.zeros(4)
xi = np.zeros((4, 5))
with h5py.File(filename, 'r') as f:
    gain[0] = f['gain_000001'][()]
    gain[1] = f['gain_000002'][()]
    gain[2] = f['gain_000003'][()]
    gain[3] = f['gain_000004'][()]
    r = f['tod'][()] - gain[:, None] * f['s_tot'][()]
    n_corr = f['n_corr'][()]
    mask = f['mask'][()]
    xi[0] = f['xi_n_000001'][()]
    xi[1] = f['xi_n_000002'][()]
    xi[2] = f['xi_n_000003'][()]
    xi[3] = f['xi_n_000004'][()]

    # detector = 3
    for detector in range(4):

        Nobs = 30  #for W-band

        dt = 1.536 / Nobs
        samprate = 1 / dt

        n_samples = len(r[0])

        # dt = 1 / samprate  # seconds
        sigma0 = xi[detector, 0] #/ gain[detector]

        r = r * mask + (n_corr + np.random.randn(*n_corr.shape) * sigma0) * (1 - mask)
        
        fknee = xi[detector, 1]
        alpha = xi[detector, 2]
        f = np.abs(fft.rfftfreq(n_samples) * samprate)


        psd = np.abs(fft.rfft(r[detector])) ** 2 / (n_samples * samprate)
        psd_n = np.abs(fft.rfft(n_corr[detector])) ** 2 / (n_samples * samprate)
        psd_r = np.abs(fft.rfft(r[detector] - n_corr[detector])) ** 2 / (n_samples * samprate)
        cut = 3
        slope = 0.7
        maxbin = 4000
        psd3 = bin_data_maxbin(psd, slope=slope, cut=cut, maxbin=maxbin)
        psd3_n = bin_data_maxbin(psd_n, slope=slope, cut=cut, maxbin=maxbin)
        psd3_r = bin_data_maxbin(psd_r, slope=slope, cut=cut, maxbin=maxbin)
        f3 = bin_data_maxbin(f, slope=slope, cut=cut, maxbin=maxbin)

        wn = sigma0 ** 2 / samprate # sigma0_def ** 2

        S = wn * (f3 / fknee) ** (alpha)
        full = wn + S
        plt.figure(figsize=(7, 5))
        plt.loglog(f3, psd3, 'C0', alpha=0.5, label=r'$d - g s_\mathrm{tot}$')
        plt.loglog(f3, psd3_r, 'k', alpha=0.5, label=r'$d - g s_\mathrm{tot} - n_\mathrm{corr}$')
        plt.loglog(f3, psd3_n, 'C1', alpha=0.5, label=r'$n_\mathrm{corr}$')

        plt.loglog(f3, full, '--', color='C0', label='current model (sample)')
        plt.loglog(f3, S, '--', color='C1') 
        plt.loglog(f3, full * 0 + wn, '--', color='k')
        plt.xlim(4e-6, 2e1)
        plt.legend(frameon=False)
        
        plt.ylabel('PSD [du${}^2$ Hz${}^{-1}$]')
        plt.xlabel('Frequency [Hz]')
        plt.savefig('ps_test_W4_det%i.png' % (detector+1), bbox_inches='tight', dpi=150)
        # plt.savefig('ps_test_W4_det%i.pdf' % (detector+1), bbox_inches='tight')
        # plt.show()

# bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

# Nobs_list = [12, 12, 15, 15, 20, 20, 30, 30, 30, 30]

# freqs = [23, 30, 40, 40, 60, 60, 90, 90, 90, 90]


# for freqband, band, Nobs in zip(freqs, bands, Nobs_list):
#     for pid in x:
#         # filename = path + 'tod_090-WMAP_W1_000500_samp000001.h5'
#         filename = path + 'tod_%03i-WMAP_%s_%06i_samp000001.h5' % (freqband, band, pid)