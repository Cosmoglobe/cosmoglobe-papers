# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import numpy as N
import scipy.stats
import healpy as hp
import plotly.colors as pcol
import matplotlib as mpl

cmap = "Plotly"
colors = getattr(pcol.qualitative, cmap)
colors.insert(3, colors.pop(-1))
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

width = 8.8

# Load data

nrad = 40
npar = 4
vals = N.zeros((nrad,npar,2))

c = 0
for r in ["023-WMAP_K", "030-WMAP_Ka", "040-WMAP_Q1", "040-WMAP_Q2", "060-WMAP_V1", "060-WMAP_V2", "090-WMAP_W1", "090-WMAP_W2", "090-WMAP_W3", "090-WMAP_W4"]:
    for i in range(4):
        data = N.loadtxt('noise_CG_'+r+str(i)+'_v1.dat')
        data[:,3] = data[:,0]/data[:,3] # Convert to mK
        ind = np.nonzero(data[:,0])
        for p in range(4):
            vals[c,p,0] = N.mean(data[ind,p])
            vals[c,p,1] = N.std(data[ind,p])
        #print c, vals[c,:,0]
        #print c, vals[c,:,1]
        c = c+1


vals[0:40,1,:] = vals[0:40,1,:] * 1000. # Hz to mHz

vals[0:4,3,:] = vals[0:4,3,:]   * N.sqrt(1.536/12) / N.sqrt(2.)# K
vals[4:8,3,:] = vals[4:8,3,:]   * N.sqrt(1.536/12) / N.sqrt(2.)# Ka
vals[8:16,3,:] = vals[8:16,3,:]  * N.sqrt(1.536/15) / N.sqrt(2.) # Q
vals[16:24,3,:] = vals[16:24,3,:] * N.sqrt(1.536/20) / N.sqrt(2.) # V
vals[24:40,3,:] = vals[24:40,3,:] * N.sqrt(1.536/30) / N.sqrt(2.) # W

wmap = N.zeros((20,2))
wmap[0,:] = ["0.66", "0.40"]
wmap[1,:] = ["0.75", "0.51"]
wmap[2,:] = ["0.71", "0.71"]
wmap[3,:] = ["0.72", "0.32"]
wmap[4,:] = ["0.92", "1.09"]
wmap[5,:] = ["1.02", "0.35"]
wmap[6,:] = ["0.85", "5.76"]
wmap[7,:] = ["0.99", "8.62"]
wmap[8,:] = ["1.22", "0.09"]
wmap[9,:] = ["1.11", "1.41"]
wmap[10,:] = ["0.97", "0.88"]
wmap[11,:] = ["1.10", "8.35"]
wmap[12,:] = ["1.35", "7.88"]
wmap[13,:] = ["1.61", "0.66"]
wmap[14,:] = ["1.61", "9.02"]
wmap[15,:] = ["1.72", "7.47"]
wmap[16,:] = ["1.65", "0.93"]
wmap[17,:] = ["1.86", "0.28"]
wmap[18,:] = ["1.71", "46.5"]
wmap[19,:] = ["1.65", "26.0"]

gsfc = N.zeros((20,2))
gsfc[0,:] = ["0.72", "6.13"]
gsfc[1,:] = ["0.87", "5.37"]
gsfc[2,:] = ["0.75", "1.66"]
gsfc[3,:] = ["0.77", "1.29"]
gsfc[4,:] = ["0.99", "3.21"]
gsfc[5,:] = ["0.95", "3.13"]
gsfc[6,:] = ["0.89", "1.92"]
gsfc[7,:] = ["1.04", "4.61"]
gsfc[8,:] = ["1.25", "2.56"]
gsfc[9,:] = ["1.07", "4.49"]
gsfc[10,:] = ["1.01", "2.43"]
gsfc[11,:] = ["1.13", "3.06"]
gsfc[12,:] = ["1.18", "16.2"]
gsfc[13,:] = ["1.41", "15.1"]
gsfc[14,:] = ["1.38", "1.76"]
gsfc[15,:] = ["1.44", "0.77"]
gsfc[16,:] = ["1.47", "1.84"]
gsfc[17,:] = ["1.69", "2.39"]
gsfc[18,:] = ["1.60", "8.46"]
gsfc[19,:] = ["1.43", "5.31"]

labels = ["K11", "", "K12", "", "Ka11", "", "Ka12", "", "Q11",  "","Q12",  "","Q21",  "","Q22",  "","V11",  "","V12",  "","V21",  "","V22", "","W11",  "","W12",  "","W21",  "","W22",  "","W31",  "","W32",  "","W41",  "","W42", ""]

for i in range(40):
    if N.mod(i,2) == 0 :
        print  labels[i], " & ", str(N.mod(i,2)+1), " & ", gsfc[i/2,0], " & ", wmap[i/2,0], " & ", '% 6.3f' % vals[i,3,0], "\pm", '% 6.3f' % vals[i,3,1], " & ", gsfc[i/2,1], " & ", wmap[i/2,1], " & ", '% 6.2f' % vals[i,1,0], "\pm", '% 6.2f' % vals[i,1,1], " & ", '% 6.2f' % vals[i,2,0], "\pm", '% 6.2f' % vals[i,2,1], "\cr"
    else:
        print  "\omit & ", str(N.mod(i,2)+1), " &  &  & ", '% 6.3f' % vals[i,3,0], "\pm", '% 6.3f' % vals[i,3,1], " &  &  & ", '% 6.2f' % vals[i,1,0], "\pm", '% 6.2f' % vals[i,1,1], " & ", '% 6.2f' % vals[i,2,0], "\pm", '% 6.2f' % vals[i,2,1], "\cr"
        


