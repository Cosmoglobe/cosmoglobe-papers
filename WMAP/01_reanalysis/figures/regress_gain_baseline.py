import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
from scipy.interpolate import interp1d

from tqdm import tqdm

from glob import glob

import time

from fits_to_h5 import quat_to_sky_coords, ang2pix_multiprocessing, get_psi
import ducc0.totalconvolve as totalconvolve
Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])


def get_ephem(time, planet):
    # Obtained these from https://ssd.jpl.nasa.gov/horizons.cgi
    t_p, ra_p, dec_p = np.loadtxt(f"/mn/stornext/u3/duncanwa/Commander/commander3/todscripts/wmap/{planet}_ephem.txt").T

    f_ra = interp1d(t_p, ra_p, fill_value="extrapolate")
    f_dec = interp1d(t_p, dec_p, fill_value="extrapolate")

    return f_ra(time), f_dec(time)


def get_flags(data, test=False, center=True, bands=np.arange(10)):

    Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])

    t2jd = data[1].header["TIME2JD"]

    quat = data[1].data["QUATERN"]
    Nobs_arr = Nobs_array[bands]

    ll_A, ll_B, p_A, p_B, t_list = quat_to_sky_coords(
        data,
        lonlat=True,
        center=center,
        ret_times=True,
        coord_out="C",
        Nobs_array=Nobs_arr,
        n_ind=bands,
    )

    time_majorframe = data[2].data["TIME"] + t2jd
    # time_majorframe = data[1].data['TIME'] + t2jd

    daflags = data[2].data["daflags"]

    planets = ["mars", "jupiter", "saturn", "uranus", "neptune"]
    radii = np.array(
        [
            [2.0, 3.0, 2.0, 2.0, 2.0],  # K (yr!=2)
            [1.5, 2.5, 1.5, 1.5, 1.5],  # Ka
            [1.5, 2.5, 1.5, 1.5, 1.5],  # Q1
            [1.5, 2.5, 1.5, 1.5, 1.5],  # Q2
            [1.5, 2.2, 1.5, 1.5, 1.5],  # V1
            [1.5, 2.2, 1.5, 1.5, 1.5],  # V2
            [1.5, 2.0, 1.5, 1.5, 1.5],  # W1
            [1.5, 2.0, 1.5, 1.5, 1.5],  # W2
            [1.5, 2.0, 1.5, 1.5, 1.5],  # W3
            [1.5, 2.0, 1.5, 1.5, 1.5],  # W4
        ]
    )
    # radii   = np.array([
    #          [7,   7,      7,      7,       7], #K (yr!=2)
    #          [7,   7,      7,      7,       7], #Ka
    #          [7,   7,      7,      7,       7], #Q1
    #          [7,   7,      7,      7,       7], #Q2
    #          [7,   7,      7,      7,       7], #V1
    #          [7,   7,      7,      7,       7], #V2
    #          [7,   7,      7,      7,       7], #W1
    #          [7,   7,      7,      7,       7], #W2
    #          [7,   7,      7,      7,       7], #W3
    #          [7,   7,      7,      7,       7], #W4
    #          ])

    radii = radii * (np.pi / 180)
    # I should be calculating the distances using each individual observation,
    # not the daflags shape, since that's going to be the same size as the major
    # frames, not the individual ones. Additionally, I am using the time for
    # each major frame, not the actual time for each individual observation.
    myflags = []
    daflags_copy = []
    for band in bands:
        myflags.append([])
        daflags_copy.append([])

    for b, band in enumerate(bands):
        myflags[b] = np.zeros(len(t_list[b]))
        daflags_copy[b] = np.zeros(len(t_list[b]))
        for i, p in enumerate(planets):
            # t = t_list[b] + t2jd
            t = t_list[b] + time_majorframe[0]
            ra_p, dec_p = get_ephem(t, p)
            ll_p = np.array([ra_p, dec_p])

            d_A = hp.rotator.angdist(ll_A[b].T, ll_p, lonlat=True)
            inds = d_A <= radii[band][i]
            myflags[b][inds] += 2 ** (2 * i + 1)

            d_B = hp.rotator.angdist(ll_B[b].T, ll_p, lonlat=True)
            inds = d_B <= radii[band][i]
            myflags[b][inds] += 2 ** (2 * i + 1 + 1)

    for b, band in enumerate(bands):
        Nobs = Nobs_array[band]
        for i in range(Nobs):
            daflags_copy[b][i::Nobs] = daflags[:, band]
        ind1 = daflags_copy[b] % 2 == 1
        myflags[b] = np.where(ind1, daflags_copy[b], myflags[b])

    if test:
        return daflags_copy, myflags
    else:
        return myflags


def get_sidelobe_alms(band="Q1", lmax=128, kmax=100, theta_c=0, psi=0):
    # LOS geometry extracted from program.pars
    dir_A_los = np.array(
        [
            [0.03993743194318, 0.92448267167832, -0.37912635267982],
            [-0.03836350153280, 0.92543717887494, -0.37695393578810],
            [-0.03157188095163, 0.95219265474988, -0.30386241059657],
            [0.03193385161530, 0.95220162163922, -0.30379647935526],
            [-0.03317333754910, 0.94156429439011, -0.33519577742792],
            [0.03337676771235, 0.94149468374332, -0.33537106592570],
            [-0.00918939185649, 0.93943847522010, -0.34259437583453],
            [-0.00950701394255, 0.94586439605663, -0.32442281201900],
            [0.00980040822398, 0.94576779947882, -0.32469558276581],
            [0.00980808738477, 0.93934799994236, -0.34282522723123],
        ]
    )
    dir_B_los = np.array(
        [
            [0.03794083653062, -0.92391755783762, -0.38070571212253],
            [-0.04002167684949, -0.92463440201100, -0.37874726137612],
            [-0.03340297596219, -0.95176877819247, -0.30499251475222],
            [0.03014337784306, -0.95192770480751, -0.30483605690947],
            [-0.03503633693827, -0.94094544143324, -0.33674045100040],
            [0.03144454385558, -0.94113854675448, -0.33655530968115],
            [-0.01147317267740, -0.93883247845653, -0.34418300902847],
            [-0.01159000320270, -0.94535005109668, -0.32585112047876],
            [0.00768184749607, -0.94540702221088, -0.32580139897397],
            [0.00751408106677, -0.93889226303920, -0.34412912836731],
        ]
    )

    bands = np.array(["K1", "Ka1", "Q1", "Q2", "V1", "V2", "W1", "W2", "W3", "W4"])
    SIDELOBE_DIR = "/mn/stornext/d16/cmbco/ola/wmap/ancillary_data/far_sidelobe_maps"
    # Construct sidelobe model
    sidelobe = hp.read_map(f"{SIDELOBE_DIR}/wmap_sidelobe_map_{band}_3yr_v2.fits")
    sidelobe = hp.reorder(sidelobe, n2r=True)
    # sidelobe = hp.read_map(f'{SIDELOBE_DIR}/map_{band.lower()}_sidelobes_yr1_v1.fits')
    # Higher resolution makes the interpolation from rotation less of a mess.
    sidelobe = hp.ud_grade(sidelobe, 1024)

    # Normalized such that \int B_A d\Omega = 1, converting from
    # sum(abs(sidelobe)) = 2*N_pix normalization
    beam_A = sidelobe / (4 * np.pi)
    beam_A[beam_A < 0] = 0
    beam_B = sidelobe / (4 * np.pi)
    beam_B[beam_B > 0] = 0
    beam_B = -beam_B

    # This is only possible for the 4pi beam
    # beam_A = beam_A/(sum(beam_A)*hp.nside2pixarea(2048))
    # beam_B = beam_B/(sum(beam_B)*hp.nside2pixarea(2048))


    # Angle psi is roughly the right value based on some tests

    dir_A = dir_A_los[band == bands][0]
    theta = np.arccos(dir_A[2])
    phi = np.arctan2(dir_A[1], dir_A[0])

    # Rotate so that main beam A is pointing in the z-direction
    r = hp.rotator.Rotator(rot=(phi, -theta, psi), deg=False, eulertype="Y")
    beam_A = r.rotate_map_pixel(beam_A)

    dir_B = dir_B_los[band == bands][0]
    theta = np.arccos(dir_B[2])
    phi = np.arctan2(dir_B[1], dir_B[0])

    # Rotate so that main beam B is pointing in the z-direction
    r = hp.rotator.Rotator(rot=(phi, -theta, -psi), deg=False, eulertype="Y")
    beam_B = r.rotate_map_pixel(beam_B)

    if theta_c > 0:
        pix = np.arange(len(beam_A))
        thetaphi = hp.pix2ang(hp.npix2nside(len(beam_A)), pix)
        r = hp.rotator.angdist(thetaphi, np.array([0, 0]))
        beam_A[r < theta_c * np.pi / 180] = 0
        beam_B[r < theta_c * np.pi / 180] = 0
        # hp.mollview(beam_A, rot=(0,90,0), min=0, max=0.5)
        # plt.show()

    # beam_A = hp.ud_grade(beam_A, 128)
    # beam_B = hp.ud_grade(beam_B, 128)

    blm_A = hp.map2alm(beam_A, lmax=lmax, mmax=kmax)
    blm_B = hp.map2alm(beam_B, lmax=lmax, mmax=kmax)

    # blm_A = blm_A[np.newaxis,:].astype('complex128')
    # blm_B = blm_B[np.newaxis,:].astype('complex128')
    blm_A = np.array([blm_A, blm_A * 0, blm_A * 0])
    blm_B = np.array([blm_B, blm_B * 0, blm_B * 0])

    return blm_A, blm_B


def make_dipole_alms(amp=3355, l=263.99, b=48.26, lmax=128, band="K1"):
    ipix = np.arange(12 * 512**2)
    x, y, z = hp.pix2vec(512, ipix)

    theta, phi = np.pi / 180 * l, np.pi / 2 - np.pi / 180 * b
    amps = amp * np.array(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]
    )

    dipole = x * amps[0] + y * amps[1] + z * amps[2]
    dipole = np.array([dipole, 0 * dipole, 0 * dipole])

    m = (
        hp.read_map(
            f"/mn/stornext/d16/cmbco/ola/wmap/freq_maps/wmap_iqusmap_r9_9yr_{band}_v5.fits",
            field=(0, 1, 2),
        )
    )

    slm = hp.map2alm(m + dipole*1e-3, lmax=lmax)

    return slm



bands = ['K1', 'Ka1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
#labels = ["13", "14", "23", "24"]
#
#band = bands[0]
#band_labels = [f"{band}{labels[i]}" for i in range(4)]

band = 'K1'
theta_c = 2.8
psi = 135
theta_cs = np.array([2.8, 2.5, 2.2, 2.2, 1.8, 1.8, 1.5, 1.5, 1.5, 1.5])
psis = np.array([135, 45, 135, 45, 45, 135, 135, 45, 135, 45])

x, y, z = hp.pix2vec(512, np.arange(12 * 512**2))
dip_W = -0.233 * x - 2.222 * y + 2.504 * z
lmax = 128
kmax = 100

lmax = 32
kmax = 16


W_DIR = '/mn/stornext/d16/cmbco/ola/wmap/tods'
W_DIR_CAL = f'{W_DIR}/calibrated'
W_DIR_UCAL = f'{W_DIR}/uncalibrated'
fnames_cal = glob(f'{W_DIR_CAL}/*.fits')
fnames_ucal = glob(f'{W_DIR_UCAL}/*.fits')

fnames_cal.sort()
fnames_ucal.sort()

slms = []
blm_As = []
blm_Bs = []

for bi in range(len(bands)):
    slm = make_dipole_alms(lmax=lmax, band = bands[bi])
    blm_A, blm_B = get_sidelobe_alms(
        band=bands[bi], lmax=lmax, kmax=kmax, theta_c=theta_cs[bi], psi=np.pi / 180 * psis[bi]
    )
    slms.append(slm)
    blm_As.append(blm_A)
    blm_Bs.append(blm_B)


for n in tqdm(range(len(fnames_cal))):

    t0 = time.time()

    data_cal = fits.open(fnames_cal[n])
    data_ucal = fits.open(fnames_ucal[n])
    
    band_labels = data_cal[2].columns.names[1:-6]
    
    t_cal = data_cal[2].data['time']
    t_ucal = data_ucal[2].data['time']
    
    assert np.allclose(t_cal, t_ucal), 'The time arrays are different'
    
    flags_all = [[], [], [], [], [], [], [], [], [], []]
    psi_A_all = [[], [], [], [], [], [], [], [], [], []]
    psi_B_all = [[], [], [], [], [], [], [], [], [], []]
    pix_A_all = [[], [], [], [], [], [], [], [], [], []]
    pix_B_all = [[], [], [], [], [], [], [], [], [], []]
    genflags = data_ucal[2].data["genflags"] * 2**11
    daflags0 = data_ucal[2].data["daflags"]
    #daflags = get_flags(data_ucal)

    for i in range(10):
        Nobs = Nobs_array[i]
        flags_all[i] = np.repeat(daflags0[:,i], Nobs)
    gal_A, gal_B, pol_A, pol_B = quat_to_sky_coords(data_ucal, center=True)
    psi_A = get_psi(gal_A, pol_A, band_labels[::4])
    psi_B = get_psi(gal_B, pol_B, band_labels[1::4])
    
    
    t2jd = 2.45e6
    
    
    for bi in range(10):
       
        # totalconvolver interpolator, grid in theta,phi,psi
        interp_A = totalconvolve.Interpolator(
            slms[bi], blm_As[bi], separate=False, lmax=lmax, kmax=kmax, epsilon=1e-4, nthreads=0
        )
        interp_B = totalconvolve.Interpolator(
            slms[bi], blm_Bs[bi], separate=False, lmax=lmax, kmax=kmax, epsilon=1e-4, nthreads=0
        )
        npnt = len(psi_A[bi])
        ptg = np.zeros((npnt, 3))
        ptg[:, 0] = gal_A[bi][:,0]  # theta
        ptg[:, 1] = gal_A[bi][:,1]  # phi
        ptg[:, 2] = psi_A[bi]  # psi
        res_A = interp_A.interpol(ptg)[0]
        
        ptg[:, 0] = gal_B[bi][:,0]  # theta
        ptg[:, 1] = gal_B[bi][:,1]  # phi
        ptg[:, 2] = psi_B[bi]  # psi
        res_B = interp_B.interpol(ptg)[0]
        
        sl = res_A - res_B
        
        for j in range(4):
            d_ucal = data_ucal[2].data[band_labels[4*bi+j]].flatten()
            d_cal = data_cal[2].data[band_labels[4*bi+j]].flatten()
            t = np.arange(len(d_ucal))
            #inds = ~np.logical_and(flags_all[0], 2047)
            inds = (flags_all[bi] % 2) != 1
            if not any(inds):
                continue
            
            if (n % 100 == 0):
               fig, axes = plt.subplots(sharex=True, nrows=3)
               fig.suptitle(band_labels[4*bi+j])
            for i in range(24):
                d_ucali = np.array_split(d_ucal, 24)[i]
                d_cali = np.array_split(d_cal, 24)[i]
                ti = np.array_split(t, 24)[i]
                indsi = np.array_split(inds, 24)[i]
                sli = np.array_split(sl, 24)[i]
                t_lab = np.array_split(t_ucal,24)[i] + t2jd - 2_400_000.5
   
                y = d_cali[indsi] + sli[indsi]
                A = np.vstack([y, np.ones_like(ti[indsi]), (ti[indsi]-ti[0]),
                  (ti[indsi]-ti[0])**2,
                  (ti[indsi]-ti[0])**3]).T
                try:
                    X = np.linalg.inv(A.T.dot(A)).dot(A.T.dot(d_ucali[indsi]))
                except np.linalg.LinAlgError:
                    continue
                g0, baseline, slope, b2, b3 = X
                y = d_cali + sli
                A = np.vstack([y, np.ones_like(ti), (ti-ti[0]),
                  (ti-ti[0])**2,
                  (ti-ti[0])**3]).T
               
                if (n % 100 == 0):
                    axes[0].plot(ti[indsi] - t[0], d_ucali[indsi])
                    axes[1].plot(ti[indsi] - t[0], d_cali[indsi])
                    
                    
                    axes[2].plot(ti - t[0],  (d_ucali - A.dot(X)))
                    axes[2].set_ylim([-1, 1])
               

                with open(f'{band_labels[4*bi+j]}_g0.txt', 'a') as f:
                  f.write(f'{t_lab[0]}\t{g0}\n')
                with open(f'{band_labels[4*bi+j]}_b0.txt', 'a') as f:
                  f.write(f'{t_lab[0]}\t{baseline}\n')
                with open(f'{band_labels[4*bi+j]}_b1.txt', 'a') as f:
                  f.write(f'{t_lab[0]}\t{slope}\n')
        
            if (n % 100 == 0):
                axes[0].set_ylabel('Raw [du]')
                axes[1].set_ylabel('Cal [mK]')
                axes[2].set_ylabel('Raw - g Cal [du]')
                plt.savefig(f'{band_labels[4*bi+j]}_timestreams_{n:06}.png',
                    bbox_inches='tight')
                plt.close()
