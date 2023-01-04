import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import h5py

from scipy.optimize import curve_fit

def one_over_f(freq, sigma0, fknee):
    return sigma0 ** 2 * (1.0 + fknee * freq)

path = '/mn/stornext/d16/cmbco/bp/dwatts/WMAP/faz_chains_DAs/'

x = np.arange(0, 1500, 25, dtype=int)

bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

Nobs_list = [12, 12, 15, 15, 20, 20, 30, 30, 30, 30]

freqs = [23, 30, 40, 40, 60, 60, 90, 90, 90, 90]


for freqband, band, Nobs in zip(freqs, bands, Nobs_list):

    chi2_diff_all = []
    chi2_corr_all = []
    param_all = []
    pids = []
    chi2_binned_all = []

    for pid in x:
        # filename = path + 'tod_090-WMAP_W1_000500_samp000001.h5'
        filename = path + 'tod_%03i-WMAP_%s_%06i_samp000001.h5' % (freqband, band, pid)
        print(filename)
        try:
            gain = np.zeros(4)
            with h5py.File(filename, 'r') as f:
                gain[0] = f['gain_000001'][()]
                gain[1] = f['gain_000002'][()]
                gain[2] = f['gain_000003'][()]
                gain[3] = f['gain_000004'][()]
                r = f['tod'][()] - gain[:, None] * f['s_tot'][()] - f['n_corr'][()]
                mask = f['mask'][()]
                # s_tot = f['s_tot'][0]
                # n_corr = f['n_corr'][0]
        except:
            print('Pid: %i not found' % pid)
            continue

        # tod = tod - s_tot - n_corr

        print(filename)


        # Nobs = 30

        dt = 1.536/Nobs
        samprate = 1/ dt

        dt = 1 / samprate  # seconds

        chi2_diff_list = []
        chi2_corr_list = []
        param_list = []
        chi2_binned_list = []

        for i in range(4):
            print(i)
            tod = r[i]
            
            n = len(tod)
            print(n)

            freq = fft.rfftfreq(n, dt)
            p = np.abs(fft.rfft(tod)) ** 2 / (n)

            high_lims = [2.0 * Nobs / 30.0, max(freq)+1e-5]

            bins_low = np.logspace(-5, np.log10(high_lims[0]), 11)
            bins_high = np.logspace(np.log10(high_lims[0]), np.log10(high_lims[1]), 10)
            nmodes_low = np.histogram(freq, bins=bins_low)[0]

            bin_freqs_low = np.histogram(freq, bins=bins_low, weights=freq)[0] / nmodes_low
            ps_low = np.histogram(freq, bins=bins_low, weights=p)[0] / nmodes_low
            nmodes_high = np.histogram(freq, bins=bins_high)[0]

            bin_freqs_high = np.histogram(freq, bins=bins_high, weights=freq)[0] / nmodes_high



            sigma0_diff = np.std(tod[1:] - tod[:-1]) / np.sqrt(2)


            ps_high = np.histogram(freq, bins=bins_high, weights=p)[0] / nmodes_high
            #p0 = (3, -1e-2)
            p0 = (3, 1e-2)
            popt, pcov = curve_fit(one_over_f, bin_freqs_high, ps_high, p0=p0, sigma=ps_high/np.sqrt(nmodes_high))

            sigma0 = popt[0]
            fknee = popt[1]
            
            # use estimated sigma0 to fill gaps, then calculate spectrum again
            tod = tod * mask[i] + (1 - mask[i]) * np.random.randn(n) * sigma0

            p = np.abs(fft.rfft(tod)) ** 2 / (n)
            ps_low = np.histogram(freq, bins=bins_low, weights=p)[0] / nmodes_low
            ps_high = np.histogram(freq, bins=bins_high, weights=p)[0] / nmodes_high

            p0 = (sigma0, fknee)
            popt, pcov = curve_fit(one_over_f, bin_freqs_high, ps_high, p0=p0, sigma=ps_high/np.sqrt(nmodes_high))

            sigma0 = popt[0]
            fknee = popt[1]

            chi2_binned = (ps_low / one_over_f(bin_freqs_low, sigma0, fknee) * nmodes_low - nmodes_low) / (np.sqrt(2 * nmodes_low))
            chi2_binned_diff = (ps_low / sigma0_diff ** 2 * nmodes_low - nmodes_low) / (np.sqrt(2 * nmodes_low))

            # print(sigma0, fknee)

            # plt.loglog(bin_freqs_low, ps_low)
            # plt.loglog(bin_freqs_high, ps_high)
            # plt.loglog(freq, one_over_f(freq, sigma0, fknee))
            # # plt.loglog(freq, one_over_f(freq, sigma0_diff, 1000000))
            # plt.grid()
            # plt.show()

            chi2 = (np.sum((tod / sigma0) ** 2) - n) / np.sqrt(2.0 * n)
            chi2_diff = (np.sum((tod / sigma0_diff) ** 2) - n) / np.sqrt(2.0 * n)
            n_f = len(freq)
            chi2_fourier = (np.sum(p / one_over_f(freq, sigma0, fknee)) - n_f) / (np.sqrt(2.0 * n_f))
            chi2_fourier_diff = (np.sum(p / one_over_f(freq, sigma0_diff, 1e8)) - n_f) / (np.sqrt(2.0 * n_f))

            tod_corrected = fft.irfft(fft.rfft(tod) * sigma0 / np.sqrt(one_over_f(freq, sigma0, fknee)))
            chi2_corr = (np.sum((tod_corrected / sigma0) ** 2) - n) / np.sqrt(2.0 * n)

            print(chi2_corr, chi2_diff)
            chi2_corr_list.append(chi2_corr)
            chi2_diff_list.append(chi2_diff)
            param_list.append([sigma0, fknee, sigma0_diff])
            chi2_binned_list.append([chi2_binned, chi2_binned_diff])
            print([sigma0, fknee, sigma0_diff])

        chi2_corr_all.append(chi2_corr_list)
        chi2_diff_all.append(chi2_diff_list)
        param_all.append(param_list)
        pids.append(pid)
        chi2_binned_all.append(chi2_binned_list)

    param = np.array(param_all)
    chi2_corr = np.array(chi2_corr_all)
    chi2_diff = np.array(chi2_diff_all)
    pid = np.array(pids).astype(int)
    chi2_binned = np.array(chi2_binned_all)


    np.save('param_%s' % band, param)
    np.savetxt('chi2_corr_%s.txt' % band, chi2_corr)
    np.savetxt('chi2_diff_%s.txt' % band, chi2_diff)
    np.save('chi2_binned_%s' % band, chi2_binned)
    np.savetxt('pids_%s.txt' % band, pid)
