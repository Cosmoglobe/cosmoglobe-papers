import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import matplotlib.patches as mpatches


bands = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

Nobs_list = [12, 12, 15, 15, 20, 20, 30, 30, 30, 30]

for band, Nobs in zip(bands, Nobs_list):

    chi2_corr = np.loadtxt('chi2_corr_%s.txt' % band)
    chi2_diff = np.loadtxt('chi2_diff_%s.txt' % band)
    chi2_binned = np.load('chi2_binned_%s.npy' % band)

    n_pid, n_det, n_chi2, n_bin = chi2_binned.shape

    chi2_binned_corr = chi2_binned[:, :, 0, :].reshape(n_pid * n_det, n_bin)
    chi2_binned_diff = chi2_binned[:, :, 1, :].reshape(n_pid * n_det, n_bin)

    s = 10

    dt = 1.536/Nobs
    samprate = 1/ dt
    n = 4194304 #7508448 # 4194304  # this is sometimes a different number, but that is not very important

    dt = 1 / samprate  # seconds

    freq = fft.rfftfreq(n, dt)

    bins_low = np.logspace(-5, np.log10(2.0 * Nobs / 30.0), 11)
    nmodes_low = np.histogram(freq, bins=bins_low)[0]
    bin_freqs_low = np.histogram(freq, bins=bins_low, weights=freq)[0] / nmodes_low


    labels = []
    def add_label(violin, label):
        color = violin["bodies"][0].get_facecolor().flatten()
        labels.append((mpatches.Patch(color=color), label))

    bin_widths = bins_low[1:] - bins_low[:-1]
    plt.figure()
    viol1 = plt.violinplot(chi2_binned_corr, positions=bin_freqs_low,
        widths=bin_widths/2, showextrema=False, showmedians=True)
    viol2 = plt.violinplot(chi2_binned_diff, positions=bin_freqs_low,
        widths=bin_widths/2, showextrema=False, showmedians=True)
    add_label(viol1, 'corrected')  
    add_label(viol2, 'uncorrected')  
    
    plt.xscale('log')
    plt.title('%s' % band)
    plt.ylabel(r'$\chi^2$')
    plt.xlabel('Frequency [Hz]')
    plt.ylim((-15, 15))
    plt.grid()
    plt.legend(*zip(*labels))
    plt.savefig('chi_per_frequency_%s.png' % band, bbox_inches='tight')
