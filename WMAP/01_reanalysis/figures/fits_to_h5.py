# ================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
# ================================================================================

import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
import h5py
import healpy as hp

from glob import glob
import multiprocessing as mp
from multiprocessing import Pool
from joblib import Parallel, delayed

from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import os, sys

from tqdm import tqdm

from astroquery.jplhorizons import Horizons
from astropy.time import Time
from datetime import datetime

all_band_labels = np.array(
    [
        "K113",
        "K114",
        "K123",
        "K124",
        "Ka113",
        "Ka114",
        "Ka123",
        "Ka124",
        "Q113",
        "Q114",
        "Q123",
        "Q124",
        "Q213",
        "Q214",
        "Q223",
        "Q224",
        "V113",
        "V114",
        "V123",
        "V124",
        "V213",
        "V214",
        "V223",
        "V224",
        "W113",
        "W114",
        "W123",
        "W124",
        "W213",
        "W214",
        "W223",
        "W224",
        "W313",
        "W314",
        "W323",
        "W324",
        "W413",
        "W414",
        "W423",
        "W424",
    ]
)


Nobs_array = np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30])
fknees = np.array(
    [
        0.7,
        0.5,
        1.2,
        0.6,
        1.1,
        0.6,
        3.0,
        4.0,
        1.3,
        2.2,
        1.5,
        4.0,
        16.17,
        15.05,
        9.02,
        7.47,
        1.84,
        2.39,
        46.5,
        26.0,
    ]
)  # mHz
fknees *= 1e-3

alphas = np.array(
    [
        -0.8,
        -0.7,
        -1.0,
        -0.9,
        -1.0,
        -1.0,
        -1.1,
        -1.1,
        -1.0,
        -0.9,
        -1.1,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
    ]
)









def coord_trans(pos_in, coord_in, coord_out, lonlat=False):
    if coord_in == coord_out:
        return pos_in
    r = hp.rotator.Rotator(coord=[coord_in, coord_out])
    pos_out = r(pos_in.T).T

    if lonlat:
        if pos_out.shape[1] == 2:
            return pos_out
        elif pos_out.shape[1] == 3:
            return hp.vec2dir(pos_out.T, lonlat=lonlat).T
    else:
        return pos_out


def Q2M(Q):
    """
    PURPOSE:
        Converts quaternions to rotation matrices.

    CALLING SEQUENCE:
        M =Q2M(Q)

    INPUTS:
        Q - Quaternions.  May be a 2-D array dimensioned 4xN or
            simply a vector dimensioned 4.

    OUTPUTS:
        M - Cube of attitude rotation matrices, 3x3xN (or 3x3
            if only one input quaternion).
    """
    q1 = -Q[0, :]
    q2 = -Q[1, :]
    q3 = -Q[2, :]
    q4 = Q[3, :]

    q11 = q1**2
    q22 = q2**2
    q33 = q3**2
    q44 = q4**2
    s = q11 + q22 + q33 + q44
    w = abs(s - 1.0) > 1e-5
    if sum(w) > 0:
        s = np.sqrt(s)
        q1 = q1 / s
        q2 = q2 / s
        q3 = q3 / s
        q4 = q4 / s

    q12 = q1 * q2
    q13 = q1 * q3
    q14 = q1 * q4
    q23 = q2 * q3
    q24 = q2 * q4
    q34 = q3 * q4

    M = np.zeros((len(q1), 3, 3))

    M[:, 0, 0] = q11 - q22 - q33 + q44
    M[:, 0, 1] = 2.0 * (q12 + q34)
    M[:, 0, 2] = 2.0 * (q13 - q24)
    M[:, 1, 0] = 2.0 * (q12 - q34)
    M[:, 1, 1] = -q11 + q22 - q33 + q44
    M[:, 1, 2] = 2.0 * (q23 + q14)
    M[:, 2, 0] = 2.0 * (q13 + q24)
    M[:, 2, 1] = 2.0 * (q23 - q14)
    M[:, 2, 2] = -q11 - q22 + q33 + q44

    M = np.transpose(M, [1, 2, 0])

    return M


def gamma_from_pol(gal, pol):
    # gal and pol are galactic lonlat vectors
    dir_gal = hp.ang2vec(gal[:, 0] % np.pi, gal[:, 1] % (2 * np.pi), lonlat=False)
    dir_pol = hp.ang2vec(pol[:, 0] % np.pi, pol[:, 1] % (2 * np.pi), lonlat=False)

    dir_Z = np.array([0, 0, 1])

    sin_theta = np.sqrt(dir_gal[:, 0] ** 2 + dir_gal[:, 1] ** 2)

    dir_west_x = dir_gal[:, 1] / sin_theta
    dir_west_y = -dir_gal[:, 0] / sin_theta
    dir_west_z = dir_gal[:, 1] * 0
    dir_west = np.array([dir_west_x, dir_west_y, dir_west_z]).T
    dir_north = (dir_Z - dir_gal[:, 2][:, np.newaxis] * dir_gal) / sin_theta[
        :, np.newaxis
    ]

    # If north Galactic pole is observed
    ind = sin_theta == 0
    dir_west[ind] = np.array([1, 0, 0])
    dir_north[ind] = np.array([0, 1, 0])

    sin_gamma = (
        dir_pol[:, 0] * dir_west[:, 0]
        + dir_pol[:, 1] * dir_west[:, 1]
        + dir_pol[:, 2] * dir_west[:, 2]
    )
    cos_gamma = (
        dir_pol[:, 0] * dir_north[:, 0]
        + dir_pol[:, 1] * dir_north[:, 1]
        + dir_pol[:, 2] * dir_north[:, 2]
    )

    cos_2_gamma = 2 * cos_gamma**2 - 1
    sin_2_gamma = 2 * sin_gamma * cos_gamma

    return sin_gamma, cos_gamma


def quat_to_sky_coords(
    data,
    center=True,
    lonlat=False,
    nointerp=False,
    ret_times=False,
    coord_out="G",
    Nobs_array=np.array([12, 12, 15, 15, 20, 20, 30, 30, 30, 30]),
    n_ind=np.arange(10),
):
    """
    Quaternion is of form (N_frames, 30, 4), with one redundant frame at the
    beginning and two redundant ones at the end, that match the adjacent frames.
    """
    quat = data[1].data["QUATERN"]
    t_major = data[2].data["time"]
    # print(data[1].data['time'][0] - data[2].data['time'][0])
    nt = len(quat)
    Q = np.zeros((4, 33, nt))
    q0 = quat[:, 0::4]
    q1 = quat[:, 1::4]
    q2 = quat[:, 2::4]
    q3 = quat[:, 3::4]
    q0 = np.array([q0[0, 0]] + q0[:, 1:-2].flatten().tolist() + q0[-1, -2:].tolist())
    q1 = np.array([q1[0, 0]] + q1[:, 1:-2].flatten().tolist() + q1[-1, -2:].tolist())
    q2 = np.array([q2[0, 0]] + q2[:, 1:-2].flatten().tolist() + q2[-1, -2:].tolist())
    q3 = np.array([q3[0, 0]] + q3[:, 1:-2].flatten().tolist() + q3[-1, -2:].tolist())
    Q = np.zeros((4, 30 * nt + 3))
    Q[0] = q0
    Q[1] = q1
    Q[2] = q2
    Q[3] = q3
    t0 = np.arange(0, 30 * nt + 3)

    # inds_0 = np.arange(len(t_major))
    # inds_1 = np.arange(len(t_major)*30)
    # t_major = interp1d(inds_0, t_major, fill_value='extrapolate')(inds_1)

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

    dir_A_pol = np.array(
        [
            [
                0.69487757242271,
                -0.29835139515692,
                -0.65431766318192,
            ],
            [
                -0.69545992357813,
                -0.29560553030986,
                -0.65494493291187,
            ],
            [
                0.71383872060219,
                -0.19131247543171,
                -0.67367189173456,
            ],
            [
                -0.71390969181845,
                -0.19099503229669,
                -0.67368675923286,
            ],
            [
                -0.69832280289930,
                -0.26176968417604,
                -0.66619959126169,
            ],
            [
                0.69826122350352,
                -0.26204606404493,
                -0.66615548040223,
            ],
            [
                0.70944248806767,
                -0.23532277684296,
                -0.66431509603747,
            ],
            [
                -0.70476543555624,
                -0.23649685267332,
                -0.66886091193973,
            ],
            [
                0.70468980214241,
                -0.23690904054153,
                -0.66879472879665,
            ],
            [-0.70959923775957, -0.23501806310177, -0.66425554705017],
        ]
    )
    dir_B_pol = np.array(
        [
            [
                0.69546590081501,
                0.29798590641998,
                -0.65385899120425,
            ],
            [
                -0.69486414021667,
                0.29814186328140,
                -0.65442742607568,
            ],
            [
                0.71423586688235,
                0.19072845484161,
                -0.67341650037147,
            ],
            [
                -0.71357469183546,
                0.19306390125546,
                -0.67345192048426,
            ],
            [
                -0.69775710213559,
                0.26425762446771,
                -0.66580998365151,
            ],
            [
                0.69876566230957,
                0.26145991550208,
                -0.66585678772745,
            ],
            [
                0.71002796142313,
                0.23471528678222,
                -0.66390438178103,
            ],
            [
                -0.70422900931886,
                0.23906270891214,
                -0.66851366750529,
            ],
            [
                0.70521159225086,
                0.23611413753036,
                -0.66852578425466,
            ],
            [-0.70903152581832, 0.23766935833457, -0.66391834701609],
        ]
    )

    M = Q2M(Q)
    M = np.transpose(M, [2, 0, 1])

    gal_A = []
    pol_A = []
    gal_B = []
    pol_B = []
    t_list = []
    for n, Nobs in zip(n_ind, Nobs_array):
        # for each group from 0--4, the interpolation is valid between 1.5--2.5,
        # which is equivalent to cutting out the first 1.5 time units from the
        # beginning of the total array and the final set of quaternions does not
        # need the last half of the time interval.
        # for i in range(nt):
        #     for j in range(30):
        #          for k in range(Nobs):
        #               offset = 1 + (k+0.5)/Nobs
        #               or
        #               offset = k/Nobs + 1 + 0.5/Nobs
        #               interp(qt, offset, qout)
        if nointerp:
            M2 = M[1:-2]
            Npts = 30 * nt
        else:
            Npts = 30 * nt * Nobs
            k = np.arange(Npts)
            if center:
                t = 1 + (k + 0.5) / Nobs
                # t = 1+k/Nobs + 0.5
                # t = np.arange(0, 30*nt, 1/Nobs) + 0.5
            else:
                t = 1 + k / Nobs
                # t = np.arange(0, 30*nt, 1/Nobs)

            M2 = np.zeros((len(t), 3, 3))
            for i in range(3):
                for j in range(3):
                    inds = np.isfinite(M[:, i, j])
                    f = interp1d(
                        t0[inds],
                        M[:, i, j][inds],
                        kind="cubic",
                        fill_value="extrapolate",
                    )
                    M2[:, i, j] = f(t)
            # T is the sample number here, but I want the actual time elapsed
            t *= 1.536 / 3600 / 24
            t_list.append(t)

        dir_A_los_cel = []
        dir_B_los_cel = []
        dir_A_los_cel = np.sum(
            M2 * np.tile(dir_A_los[n, np.newaxis, np.newaxis, :], (Npts, 3, 1)), axis=2
        )
        dir_B_los_cel = np.sum(
            M2 * np.tile(dir_B_los[n, np.newaxis, np.newaxis, :], (Npts, 3, 1)), axis=2
        )

        dir_A_los_gal = coord_trans(dir_A_los_cel, "C", coord_out)
        Pll_A = np.array(hp.vec2ang(dir_A_los_gal, lonlat=lonlat))

        dir_B_los_gal = coord_trans(dir_B_los_cel, "C", coord_out)
        Pll_B = np.array(hp.vec2ang(dir_B_los_gal, lonlat=lonlat))
        gal_A.append(Pll_A.T)
        gal_B.append(Pll_B.T)

        dir_A_pol_cel = np.sum(
            M2 * np.tile(dir_A_pol[n, np.newaxis, np.newaxis, :], (Npts, 3, 1)), axis=2
        )
        dir_B_pol_cel = np.sum(
            M2 * np.tile(dir_B_pol[n, np.newaxis, np.newaxis, :], (Npts, 3, 1)), axis=2
        )

        dir_A_pol_gal = coord_trans(dir_A_pol_cel, "C", coord_out)
        Pll_A = np.array(hp.vec2ang(dir_A_pol_gal, lonlat=lonlat))

        dir_B_pol_gal = coord_trans(dir_B_pol_cel, "C", coord_out)
        Pll_B = np.array(hp.vec2ang(dir_B_pol_gal, lonlat=lonlat))
        pol_A.append(Pll_A.T)
        pol_B.append(Pll_B.T)

    if ret_times:
        return gal_A, gal_B, pol_A, pol_B, t_list
    else:
        return gal_A, gal_B, pol_A, pol_B


def get_psi_band(gal, pol):
    psi = []
    sing, cosg = gamma_from_pol(gal, pol)
    return np.arctan2(sing, cosg)


def get_psi(gal, pol, band_labels):
    psi = []
    for band in range(len(band_labels)):
        sing, cosg = gamma_from_pol(gal[band], pol[band])
        psi.append(np.arctan2(sing, cosg))
    return psi


def ang2pix_multiprocessing(nside, theta, phi):
    return hp.ang2pix(nside, theta, phi)
