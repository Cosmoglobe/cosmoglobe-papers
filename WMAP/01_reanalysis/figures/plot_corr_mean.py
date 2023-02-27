import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import corner
import pandas as pd
import seaborn as sns

bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']

n_pids = [1066, 1069, 1330, 1328, 897, 895, 1346, 1327, 1329, 1328]

n_dets = 4
n_stats = 7
folder = 'plots/'


data = np.load('noise_data_sampavg_WMAP_2chain.npy', allow_pickle=True)

det_str = ['13', '14', '23', '24']
stats = (1, 3, 5)  # this is alpha, fknee and sigma, in that order


n_stats_used = len(stats)


# column and row labels of corr-plot
full_list = [[band + det_str[j] for stat in stats for j in range(n_dets)] for band in bands]

# binning to a common mjd-grid
n_mjd_bins = 500
mjd_bins = np.linspace(52100, 55500, n_mjd_bins+1)
all_a = []

mask = np.zeros(n_mjd_bins) + 1.0

for band in bands:
    band_data = data.item()[band]
    band_a = np.zeros((n_dets, n_stats, n_mjd_bins))

    mjd = band_data[0, 6, :]
    n_hit = np.histogram(mjd, bins=mjd_bins)[0]
    for i in range(n_stats):
        for j in range(n_dets):
            band_a[j, i, :] = np.histogram(mjd, bins=mjd_bins, weights=band_data[j, i, :])[0] / n_hit
    all_a.append(band_a)
    mask[(n_hit == 0)] = 0.0


i = 0

wh = np.where(mask == 1.0)[0]

n_pid = len(wh)
dfs = []
for a in all_a:

    a[:, 5] /=  a[:, 4]  # calibrate 
    a = a[:, stats, :].transpose((1, 0, 2))
    a = a[:, :, wh]
    a = a.reshape((n_dets * n_stats_used, n_pid))

    df = pd.DataFrame(a.T, columns=full_list[i])
    dfs.append(df)

    i += 1


df = pd.concat(dfs, axis=1, ignore_index=False, keys=None,
          levels=None, names=None, verify_integrity=False, copy=True)

corr = df.corr()



f, ax = plt.subplots(figsize=(18, 15))
#plt.title('correlation of statistics')
sns.set(font_scale=1.4)
sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool), cmap=sns.diverging_palette(220, 10, as_cmap=True),
            square=False, ax=ax, vmin=-1, vmax=1, cbar_kws={'label': 'Correlation of noise parameters'})

ax = plt.gca()
pos = ax.get_position()

dx = pos.width/10
x0 = pos.x0 + pos.width/30
dy = pos.height/10
y0 = pos.y0 + pos.height/25

fontsize = 18

labels = [r'$\alpha$', r'$f_\mathrm{knee}$', r'$\sigma_0$']*10
ax.set_xticks(np.arange(2,120, 4))
ax.set_yticks(np.arange(2,120, 4))
ax.set_xticklabels(labels, rotation=0)
ax.set_yticklabels(labels, rotation=90)

for i in range(len(bands)):
    plt.text(x0 + i*dx, 0.07, r'\textit{' + bands[i] + '}', transform=f.transFigure, fontsize=fontsize, rotation=0)
    plt.text(0.085, y0 + i*dy ,   r'\textit{' + bands[len(bands) - i - 1] + '}', transform=f.transFigure, fontsize=fontsize, rotation=90)


xl = np.array([0.0 - 0.01, 12.0 - 0.05, 30.0, 66.03]) * 0.999
xl = []
for i in range(11):
    xl.append(i*12*0.999)
plt.vlines(xl, ymin=0.0, ymax=120.0, linestyle='-', color='k', lw=1.05, alpha=1.0)
plt.hlines(xl, xmin=0.0, xmax=120.0, linestyle='-', color='k', lw=1.05, alpha=1.0)
xl = []
for i in range(41):
    #xl.append(i*12*0.999)
    xl += (np.array([4.0 - 0.01 + i*12, 8.0- 0.01 + i*12]) * 0.999).tolist()
#xl = np.array([4.0 - 0.01, 8.0- 0.01, 18.0, 24.0, 42.0, 54.0]) * 0.999
plt.vlines(xl, ymin=0.0, ymax=120.0, linestyle='-', color='k', lw=0.5, alpha=0.6)
plt.hlines(xl, xmin=0.0, xmax=120.0, linestyle='-', color='k', lw=0.5, alpha=0.6)
#plt.tight_layout()
plt.savefig('noise_parameter_correlation.pdf', bbox_inches='tight')
# plt.show()
