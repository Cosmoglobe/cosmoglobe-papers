from setup_matplotlib import *
import healpy as hp
import matplotlib.cm as cm
import numpy as N
import matplotlib.colors as col

directory="./"

for filename, outfile, coltype, vmin, vmax, dataset, freq, unit, scale, offset, comp, labelpos, labels, bar in [
    #(directory+"BP_030_IQU_full_n0512_BP7.fits", "bp_030_map_full_n0512_T.pdf", 0, -3400., 3400., r"30 GHz", r"$\left<T\right>$", r"$\mathrm{\mu K}$", 1, 0, [-3400,0,3400], [r"$-3400$", r"0", r"$3400$"], 0),
    #(directory+"BP_030_IQU_full_n0512_BP7_1deg.fits", "bp_030_map_full_n0512_1deg_Q.pdf", 0, -30.,     30., r"30 GHz", r"$\left<Q\right>$", r"$\mathrm{\mu K}$", 1, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"BP_030_IQU_full_n0512_BP7_1deg.fits", "bp_030_map_full_n0512_1deg_U.pdf", 0, -30.,     30., r"30 GHz", r"$\left<U\right>$", r"$\mathrm{\mu K}$", 1, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),

    #(directory+"BP_044_IQU_full_n0512_BP7.fits", "bp_044_map_full_n0512_T.pdf", 0, -3400., 3400., r"44 GHz", r"$\left<T\right>$", r"$\mathrm{\mu K}$", 1, 0, [-3400,0,3400], [r"$-3400$", r"0", r"$3400$"], 0),
    #(directory+"BP_044_IQU_full_n0512_BP7_1deg.fits", "bp_044_map_full_n0512_1deg_Q.pdf", 0, -30.,     30., r"44 GHz", r"$\left<Q\right>$", r"$\mathrm{\mu K}$", 1, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),
    #(directory+"BP_044_IQU_full_n0512_BP7_1deg.fits", "bp_044_map_full_n0512_1deg_U.pdf", 0, -30.,     30., r"44 GHz", r"$\left<U\right>$", r"$\mathrm{\mu K}$", 1, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 0),

    #(directory+"BP_070_IQU_full_n1024_BP7.fits", "bp_070_map_full_n1024_T.pdf", 0, -3400., 3400., r"70 GHz", r"$\left<T\right>$", r"$\mathrm{\mu K}$", 1, 0, [-3400,0,3400], [r"$-3400$", r"0", r"$3400$"], 1),
    #(directory+"BP_070_IQU_full_n1024_BP7_1deg.fits", "bp_070_map_full_n1024_1deg_Q.pdf", 0, -30.,     30., r"70 GHz", r"$\left<Q\right>$", r"$\mathrm{\mu K}$", 1, 1, [-30,0,30], [r"$-30$", r"0", r"$30$"], 1),
    #(directory+"BP_070_IQU_full_n1024_BP7_1deg.fits", "bp_070_map_full_n1024_1deg_U.pdf", 0, -30.,     30., r"70 GHz", r"$\left<U\right>$", r"$\mathrm{\mu K}$", 1, 2, [-30,0,30], [r"$-30$", r"0", r"$30$"], 1),        
        
    #(directory+"BP_030_IQU_full_n0512_BP7.fits", "bp_030_stddev_full_n0512_T.pdf", 0, 0., 20., r"30 GHz", r"$\sigma_T$", r"$\mathrm{\mu K}$", 1, 4, [0,20], [r"$0$", r"20"], 0),
    #(directory+"BP_030_IQU_full_n0512_BP7.fits", "bp_030_stddev_full_n0512_Q.pdf", 0, 0., 20., r"30 GHz", r"$\sigma_{Q}$", r"$\mathrm{\mu K}$", 1, 5, [0,20], [r"$0$", r"20"], 0),
    #(directory+"BP_030_IQU_full_n0512_BP7.fits", "bp_030_stddev_full_n0512_U.pdf", 0, 0., 20., r"30 GHz", r"$\sigma_{U}$", r"$\mathrm{\mu K}$", 1, 6, [0,20], [r"$0$", r"20"], 0),

    #(directory+"BP_044_IQU_full_n0512_BP7.fits", "bp_044_stddev_full_n0512_T.pdf", 0, 0., 20., r"44 GHz", r"$\sigma_T$", r"$\mathrm{\mu K}$", 1, 4, [0,20], [r"$0$", r"20"], 0),
    #(directory+"BP_044_IQU_full_n0512_BP7.fits", "bp_044_stddev_full_n0512_Q.pdf", 0, 0., 20., r"44 GHz", r"$\sigma_{Q}$", r"$\mathrm{\mu K}$", 1, 5, [0,20], [r"$0$", r"20"], 0),
    #(directory+"BP_044_IQU_full_n0512_BP7.fits", "bp_044_stddev_full_n0512_U.pdf", 0, 0., 20., r"44 GHz", r"$\sigma_{U}$", r"$\mathrm{\mu K}$", 1, 6, [0,20], [r"$0$", r"20"], 0),

    #(directory+"BP_070_IQU_full_n1024_BP7.fits", "bp_070_stddev_full_n1024_T.pdf", 0, 0., 20., r"70 GHz", r"$\sigma_T$", r"$\mathrm{\mu K}$", 1, 4, [0,20], [r"$0$", r"20"], 1),
    #(directory+"BP_070_IQU_full_n1024_BP7.fits", "bp_070_stddev_full_n1024_Q.pdf", 0, 0., 20., r"70 GHz", r"$\sigma_{Q}$", r"$\mathrm{\mu K}$", 1, 5, [0,20], [r"$0$", r"20"], 1),
    #(directory+"BP_070_IQU_full_n1024_BP7.fits", "bp_070_stddev_full_n1024_U.pdf", 0, 0., 20., r"70 GHz", r"$\sigma_{U}$", r"$\mathrm{\mu K}$", 1, 6, [0,20], [r"$0$", r"20"], 1),        

    (directory+"samplediff_tod_map_030_BP8.fits", "bp_030_sampdiff_full_T_BP8.pdf", 0, -3., 3., r"30 GHz", r"$\Delta T$", r"$\mathrm{\mu K}$", 1, 2, 0, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    (directory+"samplediff_tod_map_030_BP8.fits", "bp_030_sampdiff_full_Q_BP8.pdf", 0, -3., 3., r"30 GHz", r"$\Delta Q$", r"$\mathrm{\mu K}$", 1, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    (directory+"samplediff_tod_map_030_BP8.fits", "bp_030_sampdiff_full_U_BP8.pdf", 0, -3., 3., r"30 GHz", r"$\Delta U$", r"$\mathrm{\mu K}$", 1, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),

    #(directory+"samplediff_tod_map_044_BP8.fits", "bp_044_sampdiff_full_T_BP8.pdf", 0, -3., 3., r"44 GHz", r"$\Delta T$", r"$\mathrm{\mu K}$", 1, 1.0, 0, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    #(directory+"samplediff_tod_map_044_BP8.fits", "bp_044_sampdiff_full_Q_BP8.pdf", 0, -3., 3., r"44 GHz", r"$\Delta Q$", r"$\mathrm{\mu K}$", 1, 0, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),
    #(directory+"samplediff_tod_map_044_BP8.fits", "bp_044_sampdiff_full_U_BP8.pdf", 0, -3., 3., r"44 GHz", r"$\Delta U$", r"$\mathrm{\mu K}$", 1, 0, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 0),

    #(directory+"samplediff_tod_map_070_BP8.fits", "bp_070_sampdiff_full_T_BP8.pdf", 0, -3., 3., r"70 GHz", r"$\Delta T$", r"$\mathrm{\mu K}$", 1, 2.1, 0, [-3,3], [r"$-3$", r"$3$"], 1),
    #(directory+"samplediff_tod_map_070_BP8.fits", "bp_070_sampdiff_full_Q_BP8.pdf", 0, -3., 3., r"70 GHz", r"$\Delta Q$", r"$\mathrm{\mu K}$", 1, 0, 1, [-3,3], [r"$-3$", r"$3$"], 1),
    #(directory+"samplediff_tod_map_070_BP8.fits", "bp_070_sampdiff_full_U_BP8.pdf", 0, -3., 3., r"70 GHz", r"$\Delta U$", r"$\mathrm{\mu K}$", 1, 0, 2, [-3,3], [r"$-3$", r"$3$"], 1),        

    #(directory+"diff_tod_map_044_7deg.fits", "bp_044_sampdiff_full_n0512_T.pdf", 0, -3., 3., r"", r"$T^{(1)}-T^{(2)}$", r"$\mathrm{\mu K}$", 1, 0, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),
    #(directory+"diff_tod_map_044_7deg.fits", "bp_044_sampdiff_full_n0512_Q.pdf", 0, -3., 3., r"", r"$Q^{(1)}-Q^{(2)}$", r"$\mathrm{\mu K}$", 1, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),
    #(directory+"diff_tod_map_044_7deg.fits", "bp_044_sampdiff_full_n0512_U.pdf", 0, -3., 3., r"", r"$U^{(1)}-U^{(2)}$", r"$\mathrm{\mu K}$", 1, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),

    #(directory+"diff_tod_map_070_7deg.fits", "bp_070_sampdiff_full_n0512_T.pdf", 0, -3., 3., r"", r"$T^{(1)}-T^{(2)}$", r"$\mathrm{\mu K}$", 1, 0, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),
    #(directory+"diff_tod_map_070_7deg.fits", "bp_070_sampdiff_full_n0512_Q.pdf", 0, -3., 3., r"", r"$Q^{(1)}-Q^{(2)}$", r"$\mathrm{\mu K}$", 1, 1, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),
    #(directory+"diff_tod_map_070_7deg.fits", "bp_070_sampdiff_full_n0512_U.pdf", 0, -3., 3., r"", r"$U^{(1)}-U^{(2)}$", r"$\mathrm{\mu K}$", 1, 2, [-3,0,3], [r"$-3$", r"0", r"$3$"], 1),        


    #(directory+"tod_030_rms_c0001_k000061.fits", "bp_030_rmswn_full_n0512_T_BP8.pdf", 0, 0., 70., r"30 GHz", r"$\sigma_T^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 0, [0,70], [r"$0$", r"70"], 0),        
    #(directory+"tod_030_rms_c0001_k000061.fits", "bp_030_rmswn_full_n0512_Q_BP8.pdf", 0, 0., 70., r"30 GHz", r"$\sigma_Q^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 1, [0,70], [r"$0$", r"70"], 0),        
    #(directory+"tod_030_rms_c0001_k000061.fits", "bp_030_rmswn_full_n0512_U_BP8.pdf", 0, 0., 70., r"30 GHz", r"$\sigma_U^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 2, [0,70], [r"$0$", r"70"], 0),

    #(directory+"tod_044_rms_c0001_k000061.fits", "bp_044_rmswn_full_n0512_T_BP8.pdf", 0, 0., 70., r"44 GHz", r"$\sigma_T^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 0, [0,70], [r"$0$", r"70"], 0),        
    #(directory+"tod_044_rms_c0001_k000061.fits", "bp_044_rmswn_full_n0512_Q_BP8.pdf", 0, 0., 70., r"44 GHz", r"$\sigma_Q^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 1, [0,70], [r"$0$", r"70"], 0),        
    #(directory+"tod_044_rms_c0001_k000061.fits", "bp_044_rmswn_full_n0512_U_BP8.pdf", 0, 0., 70., r"44 GHz", r"$\sigma_U^{\mathrm{wn}}$", r"$\mathrm{\mu K}$", 1, 0, 2, [0,70], [r"$0$", r"70"], 0),

    #(directory+"tod_070_rms_c0001_k000061.fits", "bp_070_rmswn_full_n1024_T_BP8.pdf", 0, 0., 70., r"70 GHz", r"$\sigma_T^{\mathrm{wn}}/2$", r"$\mathrm{\mu K}$", 0.5, 0, 0, [0,70], [r"$0$", r"70"], 1),        
    #(directory+"tod_070_rms_c0001_k000061.fits", "bp_070_rmswn_full_n1024_Q_BP8.pdf", 0, 0., 70., r"70 GHz", r"$\sigma_Q^{\mathrm{wn}}/2$", r"$\mathrm{\mu K}$", 0.5, 0, 1, [0,70], [r"$0$", r"70"], 1),        
    #(directory+"tod_070_rms_c0001_k000061.fits", "bp_070_rmswn_full_n1024_U_BP8.pdf", 0, 0., 70., r"70 GHz", r"$\sigma_U^{\mathrm{wn}}/2$", r"$\mathrm{\mu K}$", 0.5, 0, 2, [0,70], [r"$0$", r"70"], 1),                        
        
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


    
    #colombi1_cmap.set_bad("gray") # color of missing pixels
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
        m.mask = np.logical_not(hp.read_map("dx12_v3_common_mask_pol_005a_0032.fits",0))
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

    for width in [7.]:
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
            plt.text(4.65,  1.2, r"%s" % freq, ha='center', va='center')
            plt.text(-4.65,  1.2, r"%s" % dataset, ha='center', va='center')
            plt.savefig(outfile, bbox_inches='tight', pad_inches=0.02)
            #plt.savefig(outfile, pad_inches=0.02)

