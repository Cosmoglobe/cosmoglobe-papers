# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
import plotly.colors as pcol
import matplotlib as mpl
import cosmoglobe as cg

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

DIR = "/mn/stornext/d5/data/duncanwa/WMAP/v1"

c1 = cg.Chain(f'{DIR}/CG_c0001_v1.h5')
c2 = cg.Chain(f'{DIR}/CG_c0002_v1.h5')

width = 8.8

# Load data
trans_1_K = c1.get('tod/023-WMAP_K/x_im')
trans_1_Ka = c1.get('tod/030-WMAP_Ka/x_im')
trans_1_Q1 = c1.get('tod/040-WMAP_Q1/x_im')
trans_1_Q2 = c1.get('tod/040-WMAP_Q2/x_im')
trans_1_V1 = c1.get('tod/060-WMAP_V1/x_im')
trans_1_V2 = c1.get('tod/060-WMAP_V2/x_im')
trans_1_W1 = c1.get('tod/090-WMAP_W1/x_im')
trans_1_W2 = c1.get('tod/090-WMAP_W2/x_im')
trans_1_W3 = c1.get('tod/090-WMAP_W3/x_im')
trans_1_W4 = c1.get('tod/090-WMAP_W4/x_im')

trans_2_K =  c2.get('tod/023-WMAP_K/x_im')
trans_2_Ka = c2.get('tod/030-WMAP_Ka/x_im')
trans_2_Q1 = c2.get('tod/040-WMAP_Q1/x_im')
trans_2_Q2 = c2.get('tod/040-WMAP_Q2/x_im')
trans_2_V1 = c2.get('tod/060-WMAP_V1/x_im')
trans_2_V2 = c2.get('tod/060-WMAP_V2/x_im')
trans_2_W1 = c2.get('tod/090-WMAP_W1/x_im')
trans_2_W2 = c2.get('tod/090-WMAP_W2/x_im')
trans_2_W3 = c2.get('tod/090-WMAP_W3/x_im')
trans_2_W4 = c2.get('tod/090-WMAP_W4/x_im')

transK = np.concatenate((trans_1_K, trans_2_K))
transKa = np.concatenate((trans_1_Ka, trans_2_Ka))
transQ1 = np.concatenate((trans_1_Q1, trans_2_Q1))
transQ2 = np.concatenate((trans_1_Q2, trans_2_Q2))
transV1 = np.concatenate((trans_1_V1, trans_2_V1))
transV2 = np.concatenate((trans_1_V2, trans_2_V2))
transW1 = np.concatenate((trans_1_W1, trans_2_W1))
transW2 = np.concatenate((trans_1_W2, trans_2_W2))
transW3 = np.concatenate((trans_1_W3, trans_2_W3))
transW4 = np.concatenate((trans_1_W4, trans_2_W4))


x = np.arange(0,2)

muK = np.mean(transK, 0)
muKa = np.mean(transKa, 0)
muQ1 = np.mean(transQ1, 0)
muQ2 = np.mean(transQ2, 0)
muV1 = np.mean(transV1, 0)
muV2 = np.mean(transV2, 0)
muW1 = np.mean(transW1, 0)
muW2 = np.mean(transW2, 0)
muW3 = np.mean(transW3, 0)
muW4 = np.mean(transW4, 0)

scale = 3.
rmsK  = scale*np.std(transK, 0)
rmsKa = scale*np.std(transKa, 0)
rmsQ1 = scale*np.std(transQ1, 0)
rmsQ2 = scale*np.std(transQ2, 0)
rmsV1 = scale*np.std(transV1, 0)
rmsV2 = scale*np.std(transV2, 0)
rmsW1 = scale*np.std(transW1, 0)
rmsW2 = scale*np.std(transW2, 0)
rmsW3 = scale*np.std(transW3, 0)
rmsW4 = scale*np.std(transW4, 0)

wmapKmu = [-0.00067, 0.00536]
wmapKrms = [scale*0.00017, scale*0.00014]

wmapKamu = [0.00353, 0.00154]
wmapKarms = [scale*0.00014, scale*0.00008]

wmapQ1mu = [-0.00013, 0.00414]
wmapQ1rms = [scale*0.00046, scale*0.00025]

wmapQ2mu = [0.00756, 0.00986]
wmapQ2rms = [scale*0.00052, scale*0.00115]

wmapV1mu = [0.00053, 0.00250]
wmapV1rms = [scale*0.00020, scale*0.00057]

wmapV2mu = [0.00352, 0.00245]
wmapV2rms = [scale*0.00033, scale*0.00098]

wmapW1mu = [0.01134, 0.00173]
wmapW1rms = [scale*0.00199, scale*0.00036]

wmapW2mu = [0.01017, 0.01142]
wmapW2rms = [scale*0.00216, scale*0.00121]

wmapW3mu = [-0.00122, 0.00463]
wmapW3rms = [scale*0.00062, scale*0.00041]

wmapW4mu = [0.02311, 0.02054]
wmapW4rms = [scale*0.00380, scale*0.00202]

print("K11 & ", '% 8.5f' % muK[0], "\pm", '% 8.5f' % (1/scale * rmsK[0]), " & ", '% 8.5f' % wmapKmu[0], "\pm", '% 8.5f' % (1/scale * wmapKrms[0]), "\cr")
print("K12 & ", '% 8.5f' % muK[1], "\pm", '% 8.5f' % (1/scale * rmsK[1]), " & ", '% 8.5f' % wmapKmu[1], "\pm", '% 8.5f' % (1/scale * wmapKrms[1]), "\cr")

print("Ka11 & ", '% 8.5f' % muKa[0], "\pm", '% 8.5f' % (1/scale * rmsKa[0]), " & ", '% 8.5f' % wmapKamu[0], "\pm", '% 8.5f' % (1/scale * wmapKrms[0]), "\cr")
print("Ka12 & ", '% 8.5f' % muKa[1], "\pm", '% 8.5f' % (1/scale * rmsKa[1]), " & ", '% 8.5f' % wmapKamu[1], "\pm", '% 8.5f' % (1/scale * wmapKarms[1]), "\cr")

print("Q11 & ", '% 8.5f' % muQ1[0], "\pm", '% 8.5f' % (1/scale * rmsQ1[0]), " & ", '% 8.5f' % wmapQ1mu[0], "\pm", '% 8.5f' % (1/scale * wmapQ1rms[0]), "\cr")
print("Q12 & ", '% 8.5f' % muQ1[1], "\pm", '% 8.5f' % (1/scale * rmsQ1[1]), " & ", '% 8.5f' % wmapQ1mu[1], "\pm", '% 8.5f' % (1/scale * wmapQ1rms[1]), "\cr")
print("Q21 & ", '% 8.5f' % muQ2[0], "\pm", '% 8.5f' % (1/scale * rmsQ2[0]), " & ", '% 8.5f' % wmapQ2mu[0], "\pm", '% 8.5f' % (1/scale * wmapQ2rms[0]), "\cr")
print("Q22 & ", '% 8.5f' % muQ2[1], "\pm", '% 8.5f' % (1/scale * rmsQ2[1]), " & ", '% 8.5f' % wmapQ2mu[1], "\pm", '% 8.5f' % (1/scale * wmapQ2rms[1]), "\cr")

print("V11 & ", '% 8.5f' % muV1[0], "\pm", '% 8.5f' % (1/scale * rmsV1[0]), " & ", '% 8.5f' % wmapV1mu[0], "\pm", '% 8.5f' % (1/scale * wmapV1rms[0]), "\cr")
print("V12 & ", '% 8.5f' % muV1[1], "\pm", '% 8.5f' % (1/scale * rmsV1[1]), " & ", '% 8.5f' % wmapV1mu[1], "\pm", '% 8.5f' % (1/scale * wmapV1rms[1]), "\cr")
print("V21 & ", '% 8.5f' % muV2[0], "\pm", '% 8.5f' % (1/scale * rmsV2[0]), " & ", '% 8.5f' % wmapV2mu[0], "\pm", '% 8.5f' % (1/scale * wmapV2rms[0]), "\cr")
print("V22 & ", '% 8.5f' % muV2[1], "\pm", '% 8.5f' % (1/scale * rmsV2[1]), " & ", '% 8.5f' % wmapV2mu[1], "\pm", '% 8.5f' % (1/scale * wmapV2rms[1]), "\cr")

print("W11 & ", '% 8.5f' % muW1[0], "\pm", '% 8.5f' % (1/scale * rmsW1[0]), " & ", '% 8.5f' % wmapW1mu[0], "\pm", '% 8.5f' % (1/scale * wmapW1rms[0]), "\cr")
print("W12 & ", '% 8.5f' % muW1[1], "\pm", '% 8.5f' % (1/scale * rmsW1[1]), " & ", '% 8.5f' % wmapW1mu[1], "\pm", '% 8.5f' % (1/scale * wmapW1rms[1]), "\cr")
print("W21 & ", '% 8.5f' % muW2[0], "\pm", '% 8.5f' % (1/scale * rmsW2[0]), " & ", '% 8.5f' % wmapW2mu[0], "\pm", '% 8.5f' % (1/scale * wmapW2rms[0]), "\cr")
print("W22 & ", '% 8.5f' % muW2[1], "\pm", '% 8.5f' % (1/scale * rmsW2[1]), " & ", '% 8.5f' % wmapW2mu[1], "\pm", '% 8.5f' % (1/scale * wmapW2rms[1]), "\cr")
print("W31 & ", '% 8.5f' % muW3[0], "\pm", '% 8.5f' % (1/scale * rmsW3[0]), " & ", '% 8.5f' % wmapW3mu[0], "\pm", '% 8.5f' % (1/scale * wmapW3rms[0]), "\cr")
print("W32 & ", '% 8.5f' % muW3[1], "\pm", '% 8.5f' % (1/scale * rmsW3[1]), " & ", '% 8.5f' % wmapW3mu[1], "\pm", '% 8.5f' % (1/scale * wmapW3rms[1]), "\cr")
print("W41 & ", '% 8.5f' % muW4[0], "\pm", '% 8.5f' % (1/scale * rmsW4[0]), " & ", '% 8.5f' % wmapW4mu[0], "\pm", '% 8.5f' % (1/scale * wmapW4rms[0]), "\cr")
print("W42 & ", '% 8.5f' % muW4[1], "\pm", '% 8.5f' % (1/scale * rmsW4[1]), " & ", '% 8.5f' % wmapW4mu[1], "\pm", '% 8.5f' % (1/scale * wmapW4rms[1]), "\cr")

#df 

vmin = -110
vmax =  160

# Create the plot
fig = plt.figure(figsize=(1.4*cm2inch(width), cm2inch(width)))


fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0)



ax1 = plt.subplot2grid((1, 1), (0, 0))

plt.locator_params(nbins=5)


# x axis
#plt.hlines(0, 0, 3300)

# labels
plt.ylabel(r"Transmission imbalance, $x_{\mathrm{im}}$");
plt.xlabel(r"Radiometer"); 
ax1.yaxis.labelpad = 10*width/17.; ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

ax1.errorbar(x+0., muK, yerr=rmsK, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='red', label='Cosmoglobe')
ax1.errorbar(x+0.3, wmapKmu, yerr=wmapKrms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='red', label='WMAP')
ax1.errorbar(x+4., muKa, yerr=rmsKa, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='orange')
ax1.errorbar(x+4.3, wmapKamu, yerr=wmapKarms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='orange')
ax1.errorbar(x+8., muQ1, yerr=rmsQ1, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='green')
ax1.errorbar(x+8.3, wmapQ1mu, yerr=wmapQ1rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='green')
ax1.errorbar(x+10., muQ2, yerr=rmsQ2, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='green')
ax1.errorbar(x+10.3, wmapQ2mu, yerr=wmapQ2rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3,color='green')
ax1.errorbar(x+14., muV1, yerr=rmsV1, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='blue')
ax1.errorbar(x+14.3, wmapV1mu, yerr=wmapV1rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='blue')
ax1.errorbar(x+16., muV2, yerr=rmsV2, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='blue')
ax1.errorbar(x+16.3, wmapV2mu, yerr=wmapV2rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='blue')
ax1.errorbar(x+20., muW1, yerr=rmsW1, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='purple')
ax1.errorbar(x+20.3, wmapW1mu, yerr=wmapW1rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='purple')
ax1.errorbar(x+22., muW2, yerr=rmsW2, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='purple')
ax1.errorbar(x+22.3, wmapW2mu, yerr=wmapW2rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='purple')
ax1.errorbar(x+24., muW3, yerr=rmsW3, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='purple')
ax1.errorbar(x+24.3, wmapW3mu, yerr=wmapW3rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3, color='purple')
ax1.errorbar(x+26., muW4, yerr=rmsW4, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, color='purple')
ax1.errorbar(x+26.3, wmapW4mu, yerr=wmapW4rms, fmt='.', ms=3, capsize=1, capthick=1, elinewidth=1, alpha=0.3,  color='purple')

plt.plot([0,28], [0,0], "k", color='black', linewidth=1, linestyle='--')

# reduce ticks for small figures
#if width < 10:
#    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
plt.grid(False, which="major", axis="both")

# axes limits
#ax1.set_xscale("log")
#ax1.set_yscale("log")
#plt.ylim([-1.5, 1.5]);
#plt.xlim([0, 28]);

plt.xticks([0,1,4,5,8,9,10,11,14,15,16,17,20,21,22,23,24,25,26,27], [r"K11", r"K12", r"Ka11", r"Ka12", r"Q11", r"Q12", r"Q21", r"Q22", r"V11", r"V12", r"V21", r"V22", r"W11", r"W12", r"W21", r"W22", r"W31", r"W32", r"W41", r"W42"])
#plt.yticks([-1,-0,1], [r"$-1.0$", r"$0.0$", r"$1.0$"])
plt.setp( ax1.get_xticklabels(), visible=True)

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

#plt.text(80,2000,r"$TT$, 30 GHz", fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for ticklabel in ax1.xaxis.get_ticklabels():
    ticklabel.set_rotation(45.)    

# legend
leg = plt.legend(frameon=True, loc=2, fontsize=10)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(0)







# save to pdf with right bounding box
#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#plt.savefig("test.pdf", bbox_inches=extent, pad_inches=0.02)
plt.savefig("../../figures/x_im_CG_v1.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#plt.savefig("test.pdf", bbox_inches=[0,1,0,1], pad_inches=0.02)

# Make table

