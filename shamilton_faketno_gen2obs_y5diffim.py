import ephem
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree as KDTree
# from ccdBounds import *
import argparse


# This is a dictionary giving the approx corner positions of each CCD on
# the sky, in degrees from center of the focal plane.  Tuple is (xmin, xmax,
# ymin, ymax) with x to E and y to N.
ccdBounds = {'N1': (-1.0811, -0.782681, -0.157306, -0.00750506),
             'N2': (-0.771362, -0.472493, -0.157385, -0.00749848),
             'N3': (-0.461205, -0.161464, -0.157448, -0.00749265),
             'N4': (-0.150127, 0.149894, -0.15747, -0.00749085),
             'N5': (0.161033, 0.460796, -0.157638, -0.0074294),
             'N6': (0.472171, 0.771045, -0.157286, -0.00740563),
             'N7': (0.782398, 1.08083, -0.157141, -0.0074798),
             'N8': (-0.92615, -0.627492, -0.321782, -0.172004),
             'N9': (-0.616455, -0.317043, -0.322077, -0.172189),
             'N10': (-0.305679, -0.00571999, -0.322071, -0.17217),
             'N11': (0.00565427, 0.305554, -0.322243, -0.172254),
             'N12': (0.31684, 0.616183, -0.322099, -0.172063),
             'N13': (0.627264, 0.925858, -0.321792, -0.171887),
             'N14': (-0.926057, -0.62726, -0.485961, -0.336213),
             'N15': (-0.616498, -0.317089, -0.486444, -0.336606),
             'N16': (-0.30558, -0.00578257, -0.486753, -0.336864),
             'N17': (0.00532179, 0.305123, -0.486814, -0.33687),
             'N18': (0.316662, 0.616018, -0.486495, -0.336537),
             'N19': (0.62708, 0.92578, -0.485992, -0.336061),
             'N20': (-0.770814, -0.471826, -0.650617, -0.500679),
             'N21': (-0.460777, -0.161224, -0.650817, -0.501097),
             'N22': (-0.149847, 0.149886, -0.650816, -0.501308),
             'N23': (0.161001, 0.460566, -0.650946, -0.501263),
             'N24': (0.47163, 0.770632, -0.650495, -0.500592),
             'N25': (-0.615548, -0.316352, -0.814774, -0.665052),
             'N26': (-0.305399, -0.00591217, -0.814862, -0.665489),
             'N27': (0.00550714, 0.304979, -0.815022, -0.665418),
             'N28': (0.316126, 0.615276, -0.814707, -0.664908),
             'N29': (-0.46018, -0.16101, -0.97887, -0.829315),
             'N31': (0.160884, 0.460147, -0.978775, -0.829426),
             'S1': (-1.08096, -0.782554, 0.00715956, 0.15689),
             'S2': (-0.7713, -0.47242, 0.0074194, 0.157269),
             'S3': (-0.4611, -0.161377, 0.00723009, 0.157192),
             'S4': (-0.149836, 0.150222, 0.00737069, 0.157441),
             'S5': (0.161297, 0.461031, 0.0072399, 0.1572),
             'S6': (0.472537, 0.771441, 0.00728934, 0.157137),
             'S7': (0.782516, 1.08097, 0.00742809, 0.15709),
             'S8': (-0.92583, -0.627259, 0.171786, 0.32173),
             'S9': (-0.616329, -0.31694, 0.171889, 0.321823),
             'S10': (-0.305695, -0.00579187, 0.172216, 0.322179),
             'S11': (0.00556739, 0.305472, 0.172237, 0.322278),
             'S12': (0.316973, 0.61631, 0.172015, 0.322057),
             'S13': (0.627389, 0.925972, 0.171749, 0.321672),
             'S14': (-0.925847, -0.627123, 0.335898, 0.48578),
             'S15': (-0.616201, -0.316839, 0.336498, 0.486438),
             'S16': (-0.305558, -0.00574858, 0.336904, 0.486749),
             'S17': (0.00557115, 0.305423, 0.33675, 0.486491),
             'S18': (0.316635, 0.615931, 0.33649, 0.486573),
             'S19': (0.627207, 0.925969, 0.336118, 0.485923),
             'S20': (-0.770675, -0.471718, 0.500411, 0.65042),
             'S21': (-0.46072, -0.161101, 0.501198, 0.650786),
             'S22': (-0.149915, 0.14982, 0.501334, 0.650856),
             'S23': (0.160973, 0.460482, 0.501075, 0.650896),
             'S24': (0.47167, 0.770647, 0.50045, 0.650441),
             'S25': (-0.615564, -0.316325, 0.66501, 0.814674),
             'S26': (-0.30512, -0.0056517, 0.665531, 0.81505),
             'S27': (0.00560886, 0.305082, 0.665509, 0.815022),
             'S28': (0.316158, 0.615391, 0.665058, 0.814732),
             'S29': (-0.46021, -0.160988, 0.829248, 0.978699),
             'S30': (-0.150043, 0.149464, 0.829007, 0.978648),
             'S31': (0.160898, 0.460111, 0.82932, 0.978804)}

# The following dictionary maps CCD names to their id numbers in DESDM.
ccdNum = {'N1': 32,
          'N2': 33,
          'N3': 34,
          'N4': 35,
          'N5': 36,
          'N6': 37,
          'N7': 38,
          'N8': 39,
          'N9': 40,
          'N10': 41,
          'N11': 42,
          'N12': 43,
          'N13': 44,
          'N14': 45,
          'N15': 46,
          'N16': 47,
          'N17': 48,
          'N18': 49,
          'N19': 50,
          'N20': 51,
          'N21': 52,
          'N22': 53,
          'N23': 54,
          'N24': 55,
          'N25': 56,
          'N26': 57,
          'N27': 58,
          'N28': 59,
          'N29': 60,
          'N30': 61,
          'N31': 62,
          'S1': 25,
          'S2': 26,
          'S3': 27,
          'S4': 28,
          'S5': 29,
          'S6': 30,
          'S7': 31,
          'S8': 19,
          'S9': 20,
          'S10': 21,
          'S11': 22,
          'S12': 23,
          'S13': 24,
          'S14': 13,
          'S15': 14,
          'S16': 15,
          'S17': 16,
          'S18': 17,
          'S19': 18,
          'S20': 8,
          'S21': 9,
          'S22': 10,
          'S23': 11,
          'S24': 12,
          'S25': 4,
          'S26': 5,
          'S27': 6,
          'S28': 7,
          'S29': 1,
          'S30': 2,
          'S31': 3,
          'None': -99}

rad2deg = 180/np.pi
deg2rad = np.pi/180


def read_file(data_file):
    """
    Read in a data file. Can be any delimiter separated file that pandas can
    parse, or a pandas DataFrame itself

    :param data_file:  either the name of a file (str) or a pandas DataFrame
    :return: the pandas DataFrame of data
    """
    if isinstance(data_file, str):
        df = pd.read_table(data_file, sep=None, engine='python')
    elif isinstance(data_file, pd.DataFrame):
        df = data_file
    else:
        return None

    df.columns = df.columns.map(str.lower)

    return df


def create_tno(tno_df_row):
    """
    Create a pyEphem EllipticalBody object using a row of a pandas DataFrame
    that contains the orbital elements of a fake object:
        * a = semimajor axis in AU
        * e = eccentricity
        * i = inclination in degrees
        * aop = argument of perihelion in degrees
        * lan = longitude of ascending node in degrees
        * m = mean anomaly in degrees
        * _epoch_M = epochjd - 2415020
                   = epoch DJD  (pyEphem requires the DJD....)
                   = date at which the mean anomaly is calculated
        * epoch = coordinate system to use (set to J2000)
        * h = absolute magnitude of the object

    :param tno_df_row:  The pandas DataFrame row with the above parameters
                        at minimum
    :return:  the pyEphem EllipticalBody object
    """
    tno = ephem.EllipticalBody()
    offset_fakeID_snfake = 180000000
    tno.name = str(int(tno_df_row.fakeid) + offset_fakeID_snfake)
    tno._a = tno_df_row.a
    tno._e = tno_df_row.e
    tno._inc = tno_df_row.i
    tno._om = tno_df_row.aop
    tno._Om = tno_df_row.lan
    tno._M = tno_df_row.m
    tno._epoch_M = ephem.date(tno_df_row.epochjd - 2415020)
    tno._epoch = ephem.date('2000/01/01')
    tno._H = tno_df_row.h

    return tno


def tno_color(mag_r, band):
    """
    Toy model of TNO colors derived from solar absolute magniture and
    reflectivity slope S~25 observed for IOC bodies (Sheppard 2010)
    """
    col_index = {'g':0.9, 'r':0, 'i':-0.46, 'z':-0.84, 'Y':-1.01}
    return mag_r + col_index[band]


@np.vectorize
def compute_tno_position(tno, date):
    """
    Given a pyEphem EllipticalBody object and a date, compute the object's
    position and magnitude at the given date

    :param tno:  A pyEphem EllipticalBody object
    :param date:  A pyEphem date object
    :return:  The RA (wrapped from -pi to pi), Dec, and magnitude at the
              given date
    """
    tno.compute(date)

    return tno.a_ra.znorm, tno.a_dec, tno.mag


def get_objects_to_keep(exposure_df, tno_fakes):
    """
    Given a list of exposures and potential TNO fakes, determine which fakes
    will fall within 3 degrees of an exposure at the midpoint of the night.
    These objects are more likely to end up on a CCD at some point than the
    objects nowhere near any exposures

    :param exposure_df:  the pandas DataFrame of exposures
    :param tno_fakes:   the list of fake TNOs, given as
                        ephem.EllipticalBody() objects
    :return:  A 2x2 array containing the row indices of the nearby exposures
              for each TNO in tno_fakes
    """
    exp_centers = exposure_df[['telra', 'teldec']].values
    exp_tree = KDTree(exp_centers)
    midpoint = ephem.date(np.median(exposure_df.djd_obs.values))

    ras_midpoint, decs_midpoint, mags_midpoint = compute_tno_position(
            tno_fakes, midpoint)

    pos_at_midpoint = np.array((ras_midpoint, decs_midpoint)).T

    near_lists = np.array(exp_tree.query_ball_point(pos_at_midpoint,
                                                    3.0*np.pi/180))
    near_list_lengths = np.array([len(l) for l in near_lists])
    non_zeros = np.argwhere(near_list_lengths != 0)

    near_lists_keep = near_lists[non_zeros].flatten()
    tno_fakes_keep = tno_fakes[non_zeros].flatten()

    return near_lists_keep, tno_fakes_keep


@np.vectorize
def compute_chip(rockra, rockdec, expra, expdec):
    """
    Given the ra and dec of a point and of the center
    of an exposure, find the CCD containing that point.

    Returns a pair of the CCD name and number.
    """

    # compute difference in degrees (normalized between -180, +180)
    # the 180/pi is because ephem.Angle objects are natively in radians
    deltara = 180 / np.pi * ephem.degrees(rockra - expra).znorm
    deltadec = 180 / np.pi * ephem.degrees(rockdec - expdec).znorm

    ccdname = 'None'
    for k in ccdBounds:
        if (ccdBounds[k][0] < deltara < ccdBounds[k][1]) and \
                (ccdBounds[k][2] < deltadec < ccdBounds[k][3]):
            ccdname = k
    return ccdname, ccdNum[ccdname]


def get_good_observations(exposure_df, tno_fakes, nearby_exposure_list):
    """
    Given the exposure list, a list of TNO fakes, and a list of exposures
    those fakes are near, calculate the observations that fall onto a CCD chip

    :param exposure_df:  the pandas DataFrame of exposure info
    :param tno_fakes:  the list of TNO fakes to evaluate
    :param nearby_exposure_list:  the 2x2 array of exposures near each fake
                                  in tno_fakes
    :return:  a list of the good observations for each TNO fake, i.e. all
              observations that were determined to have fallen on a CCD chip
    """
    good_obs = []
    for idx_list, tno in zip(nearby_exposure_list, tno_fakes):
        exp_keep = exposure_df.iloc[idx_list]
        t_offset = exp_keep.exptime*ephem.second / 2
        ras, decs, mags = compute_tno_position(tno,
                                               exp_keep.djd_obs.values
                                               + t_offset)
        chip_names, chip_nums = compute_chip(ras, decs,
                                             exp_keep.telra.values,
                                             exp_keep.teldec.values)
        if (chip_nums != -99).any():
            good_chip_idxs = np.argwhere(chip_nums != -99).flatten()
            for ii in good_chip_idxs:
                if ras[ii] < 0:
                    ras[ii] += 2*np.pi
                obs = dict(fakeid=tno.name, date=exp_keep.date_obs.iloc[ii],
                           mjd_obs=exp_keep.mjd_obs.iloc[ii],
                           expnum=exp_keep.expnum.iloc[ii],
                           exptime=exp_keep.exptime.iloc[ii],
                           ra=ras[ii], dec=decs[ii],
                           ccdnum=chip_nums[ii],
                           band=exp_keep.band.iloc[ii],
                           nite=exp_keep.nite.iloc[ii],
                           mag=tno_color(mags[ii], exp_keep.band.iloc[ii]))
                good_obs.append(obs)

    return good_obs


def main():
    # Parse the command line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp-file", required=True,
                    help="The path to the file containing the exposure "
                         "information")
    ap.add_argument("--fakegen-file", required=True,
                    help="The path to the path containing the orbital "
                         "elements of the fakes to be generated.")
    ap.add_argument("--fakeobs-outfile", required=True,
                    help="The filename to write the generated fake "
                         "observations to")
    args = ap.parse_args()

    # Read in the exposures file and perform necessary alterations:
    #    Select only griz
    #    Compute the DJD (needed for pyEphem
    #    Convert telra and teldec into ephem.Angle() objects
    exposures = read_file(args.exp_file)
    exposures = exposures[exposures.band.isin(['g','r','i','z'])]
    exposures['djd_obs'] = exposures.mjd_obs - 15019.5
    exposures['telra'] = exposures.telra.apply(lambda x: ephem.hours(x).znorm)
    exposures['teldec'] = exposures.teldec.apply(ephem.degrees)

    # Read in the fake file containing orbital elements and create
    # ephem.EllipticalBody() objects
    tno_gen_df = read_file(args.fakegen_file)
    tno_fakes = tno_gen_df.apply(create_tno, axis=1).values

    # Determine the TNOs that are nearby an exposure in the list
    near_lists_keep, tno_fakes_keep = get_objects_to_keep(exposures, tno_fakes)

    # Determine which TNOs will actually fall onto a CCD and therefore be
    # embedded
    good_obs = get_good_observations(exposures, tno_fakes_keep, near_lists_keep)

    tno_observations = pd.DataFrame(good_obs)
    tno_observations.to_csv(args.fakeobs_outfile, index=False)


if __name__ == "__main__":
    main()
