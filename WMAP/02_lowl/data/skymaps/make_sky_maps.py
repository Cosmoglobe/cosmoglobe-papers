from setup_matplotlib import *
import healpy as hp
import matplotlib.cm as cm
import numpy as N
import matplotlib.colors as col

directory="./"

for filename, outfile, coltype, vmin, vmax, dataset, freq, unit, scale, offset, comp, labelpos, labels, bar in [
    # K1
    #(directory+"CG_K_n16.fits", "CG_K_n16_Q.pdf", 0, -30, 30., r"K1", r"$Q^{\mathrm{CG}}$", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_K_n16.fits", "CG_K_n16_U.pdf", 0, -30, 30., r"", r"$U^{\mathrm{CG}}$", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_K_map.fits", "WMAP_K_n16_Q.pdf", 0, -30, 30., r"", r"$Q^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_Ka_map.fits", "WMAP_K_n16_U.pdf", 0, -30, 30., r"", r"$U^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_K_10deg.fits", "diff_K_n16_10deg_Q.pdf", 0, -3, 3., r"", r"$Q^{\mathrm{CG}}-Q^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_K_10deg.fits", "diff_K_n16_10deg_U.pdf", 0, -3, 3., r"", r"$U^{\mathrm{CG}}-U^{\mathrm{WMAP}}$", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    # Ka1
    #(directory+"CG_Ka_n16.fits", "CG_Ka_n16_Q.pdf", 0, -30, 30., r"Ka1", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_Ka_n16.fits", "CG_Ka_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_Ka_map.fits", "WMAP_Ka_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_Ka_map.fits", "WMAP_Ka_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_Ka_10deg.fits", "diff_Ka_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_Ka_10deg.fits", "diff_Ka_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    # Q1
    #(directory+"CG_Q1_n16.fits", "CG_Q1_n16_Q.pdf", 0, -30, 30., r"Q1", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_Q1_n16.fits", "CG_Q1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_Q1_map.fits", "WMAP_Q1_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_Q1_map.fits", "WMAP_Q1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_Q1_10deg.fits", "diff_Q1_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_Q1_10deg.fits", "diff_Q1_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # Q2
    #(directory+"CG_Q2_n16.fits", "CG_Q2_n16_Q.pdf", 0, -30, 30., r"Q2", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_Q2_n16.fits", "CG_Q2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_Q2_map.fits", "WMAP_Q2_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_Q2_map.fits", "WMAP_Q2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_Q2_10deg.fits", "diff_Q2_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_Q2_10deg.fits", "diff_Q2_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                               
    # V1
    #(directory+"CG_V1_n16.fits", "CG_V1_n16_Q.pdf", 0, -30, 30., r"V1", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_V1_n16.fits", "CG_V1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_V1_map.fits", "WMAP_V1_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_V1_map.fits", "WMAP_V1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_V1_10deg.fits", "diff_V1_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_V1_10deg.fits", "diff_V1_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # V2
    #(directory+"CG_V2_n16.fits", "CG_V2_n16_Q.pdf", 0, -30, 30., r"V2", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_V2_n16.fits", "CG_V2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_V2_map.fits", "WMAP_V2_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_V2_map.fits", "WMAP_V2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_V2_10deg.fits", "diff_V2_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_V2_10deg.fits", "diff_V2_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # W1
    #(directory+"CG_W1_n16.fits", "CG_W1_n16_Q.pdf", 0, -30, 30., r"W1", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_W1_n16.fits", "CG_W1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_W1_map.fits", "WMAP_W1_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_W1_map.fits", "WMAP_W1_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_W1_10deg.fits", "diff_W1_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_W1_10deg.fits", "diff_W1_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # W2
    #(directory+"CG_W2_n16.fits", "CG_W2_n16_Q.pdf", 0, -30, 30., r"W2", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_W2_n16.fits", "CG_W2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_W2_map.fits", "WMAP_W2_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_W2_map.fits", "WMAP_W2_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_W2_10deg.fits", "diff_W2_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_W2_10deg.fits", "diff_W2_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # W3
    #(directory+"CG_W3_n16.fits", "CG_W3_n16_Q.pdf", 0, -30, 30., r"W3", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"CG_W3_n16.fits", "CG_W3_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"wmap_9yr_coadd_W3_map.fits", "WMAP_W3_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),        
    #(directory+"wmap_9yr_coadd_W3_map.fits", "WMAP_W3_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"diff_W3_10deg.fits", "diff_W3_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),        
    #(directory+"diff_W3_10deg.fits", "diff_W3_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),                                
    # W4
    #(directory+"CG_W4_n16.fits", "CG_W4_n16_Q.pdf", 0, -30, 30., r"W4", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,30], [r"$-30$", r"$30$"], 1),        
    #(directory+"CG_W4_n16.fits", "CG_W4_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,30], [r"$-30$", r"$30$"], 1),
    #(directory+"wmap_9yr_coadd_W4_map.fits", "WMAP_W4_n16_Q.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-30,30], [r"$-30$", r"$30$"], 1),        
    #(directory+"wmap_9yr_coadd_W4_map.fits", "WMAP_W4_n16_U.pdf", 0, -30, 30., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-30,30], [r"$-30$", r"$30$"], 1),
    #(directory+"diff_W4_10deg.fits", "diff_W4_n16_10deg_Q.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 1, [-3,3], [r"$-3$", r"$3$"], 1),        
    #(directory+"diff_W4_10deg.fits", "diff_W4_n16_10deg_U.pdf", 0, -3, 3., r"", r"", r"$\mathrm{\mu K}$", 1000, 0, 2, [-3,3], [r"$-3$", r"$3$"], 1),
    # Likelihood inputs
    (directory+"wt_r3_9yr.KaQV.map_q", "wt_r3_9yr.KaQV.map_q.pdf", 0, -3000, 3000., r"", r"", r"$\mathrm{\mu K}$", 1, 0, 0, [-3000,3000], [r"$-3000$", r"$3000$"], 1),                                        
]:
    m = hp.ma(hp.read_map(filename,comp))*scale - offset
    if coltype == 2:
        ind = N.where(m > 0)
        m[ind] = N.log10(0.5*(m[ind]+N.sqrt(4.+m[ind]*m[ind])))
        ind = N.where(m < 0)
        m[ind] = -N.log10(0.5*(N.abs(m[ind])+N.sqrt(4.+m[ind]*m[ind])))
    nside = hp.npix2nside(len(m))

    # setup colormap
    from matplotlib.colors import ListedColormap
    colombi1_cmap = ListedColormap(np.loadtxt("parchment1.dat")/255.)

    #colombi1_cmap = ListedColormap(turbo_cmap)
    #colombi1_cmap = plt.get_cmap("bwr")
    
    #startcolor = 'black'  # a dark olive 
    #midcolor = color    # a bright yellow
    #endcolor = 'white'    # medium dark red
    #colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","DeepSkyBlue","Blue","black",col1,(1.,156/256.,57/256.,1),"white"])
    #colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["gray"])
    #colombi1_cmap = plt.get_cmap("gray")
    
    #if coltype == 0:
    #    colombi1_cmap = plt.get_cmap("gray")
    #elif coltype == 1:
    #    colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","orange","black","LawnGreen","white"])
    #elif coltype == 2:
    #    colombi1_cmap = col.LinearSegmentedColormap.from_list('own2',["white","DeepSkyBlue","Blue","black",(195/256.,5/256.,0,1),(1.,156/256.,57/256.,1),"white"])


    
    colombi1_cmap.set_bad("gray") # color of missing pixels
    #colombi1_cmap.set_under("white") # color of background, necessary if you want to use
    # this colormap directly with hp.mollview(m, cmap=colombi1_cmap)

    use_mask = False

    # using directly matplotlib instead of mollview has higher
    # quality output, I plan to merge this into healpy

    # ratio is always 1/2
    xsize = 2000
    ysize = 1000

    #unit = r"$\mathrm{\mu K}$"

    # this is the mollview min and max
    #vmin = 0; vmax = 10

    theta = np.linspace(np.pi, 0, ysize)
    phi = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)

    if use_mask:
        # mask
        m.mask = np.logical_not(hp.read_map("wmap_polarization_analysis_mask_r4_9yr_v5.fits",0))
        grid_mask = m.mask[grid_pix]
        grid_map = np.ma.MaskedArray(m[grid_pix], grid_mask)
    else:
        grid_map = m[grid_pix]

    from matplotlib.projections.geo import GeoAxes

    class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
        """Shifts labelling by pi

        Shifts labelling from -180,180 to 0-360"""
        def __call__(self, x, pos=None):
            if x != 0:
                x *= -1
            if x < 0:
                x += 2*np.pi
            return GeoAxes.ThetaFormatter.__call__(self, x, pos)

    #for width in [4.5]:
    for width in [5.5]:        
        for cmap, colormaptag in [(None, ''), (colombi1_cmap, "colombi1_")]:

            fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/2.)))
            # matplotlib is doing the mollveide projection
            ax = fig.add_subplot(111,projection='mollweide')

            # remove white space around the image
            plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)

            # rasterized makes the map bitmap while the labels remain vectorial
            # flip longitude to the astro convention
            image = plt.pcolormesh(longitude[::-1], latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap)

            # graticule
            ax.set_longitude_grid(60)
            ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
            if width < 10:
                ax.set_latitude_grid(45)
                ax.set_longitude_grid_ends(90)


            # colorbar
            if bar :
                cb = fig.colorbar(image, orientation='horizontal', shrink=.4, pad=0.05, ticks=labelpos)
                cb.ax.xaxis.set_label_text(unit)
                cb.ax.xaxis.labelpad = -8
                cb.set_ticklabels(labels)
                # workaround for issue with viewers, see colorbar docstring
                cb.solids.set_edgecolor("face")
                cb.set_ticklabels(labels)

            #ax.tick_params(axis='x', labelsize=10)
            #ax.tick_params(axis='y', labelsize=10)

            # remove longitude tick labels
            ax.xaxis.set_ticklabels([])
            # remove horizontal grid
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([])

            plt.grid(True)
            plt.title(r"%s" % freq, fontsize=16)
            #plt.text(4.65,  1.2, r"%s" % freq, ha='center', va='center')
            plt.text(-4.65,  1.2, r"%s" % dataset, ha='center', va='center')
            plt.savefig(outfile, bbox_inches='tight', pad_inches=0.02)
            #plt.savefig(outfile, pad_inches=0.02)

