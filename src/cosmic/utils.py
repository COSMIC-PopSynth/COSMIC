# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017 - 2021)
#
# This file is part of cosmic.
#
# cosmic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cosmic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cosmic.  If not, see <http://www.gnu.org/licenses/>.

"""`utils`
"""
import scipy
import numpy as np
import pandas as pd
import scipy.special as ss
import astropy.stats as astrostats
import warnings
import ast
import operator
import json
import itertools
import os.path
import h5py as h5
import re

import sys
if sys.version_info >= (3, 9):
    from importlib.resources import files as io_files
else:
    from importlib_resources import files as io_files

from configparser import ConfigParser
from .bse_utils.zcnsts import zcnsts

__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = [
    "Scott Coughlin <scott.coughlin@ligo.org>",
    "Michael Zevin <zevin@northwestern.edu>",
    "Tom Wagg <tomjwagg@gmail.com>",
]
__all__ = [
    "filter_bin_state",
    "conv_select",
    "mass_min_max_select",
    "idl_tabulate",
    "rndm",
    "param_transform",
    "dat_transform",
    "dat_un_transform",
    "knuth_bw_selector",
    "error_check",
    "check_initial_conditions",
    "convert_kstar_evol_type",
    "parse_inifile",
    "pop_write",
    "a_from_p",
    "p_from_a",
    "get_Z_from_FeH",
    "get_FeH_from_Z",
    "get_binfrac_of_Z",
    "get_porb_norm",
    "get_met_dep_binfrac",
    "explain_setting",
]


def filter_bin_state(bcm, bpp, method, kstar1_range, kstar2_range):
    """Filter the output of bpp and bcm, where the kstar ranges
    have already been selected by the conv_select module

    Parameters
    ----------
    bcm : `pandas.DataFrame`
        bcm dataframe

    bpp : `pandas.DataFrame`
        bpp dataframe

    method : `dict`,
        one or more methods by which to filter the
        bpp or bcm table, e.g. ``{'binary_state' : [0,1]}``;
        This means you do *not* want to select the final state of the binaries in the bcm array

    kstar1_range : `list`
        list containing all kstar1 values to retain

    kstar2_range : `list`
        list containing all kstar2 values to retain

    Returns
    -------
    bcm : `pandas.DataFrame`
        filtered bcm dataframe
    """
    _known_methods = ["binary_state", "timestep_conditions"]

    if not set(method.keys()).issubset(set(_known_methods)):
        raise ValueError(
            "You have supplied an "
            "unknown method to filter out "
            "the bpp or bcm array. Known methods are "
            "{0}".format(_known_methods)
        )

    for meth, use in method.items():
        if meth == "binary_state":
            bin_num_save = []

            # in order to filter on binary state we need the last entry of the bcm array for each binary
            bcm_last_entry = bcm.groupby("bin_num").last().reset_index()

            # in order to find the properities of disrupted or systems
            # that are alive today we can simply check the last entry in the bcm
            # array for the system and see what its properities are today
            bcm_0_2 = bcm_last_entry.loc[(bcm_last_entry.bin_state != 1)]
            bin_num_save.extend(bcm_0_2.bin_num.tolist())

            # in order to find the properities of merged systems
            # we actually need to search in the BPP array for the properities
            # of the objects right at merge because the bcm will report
            # the post merge object only
            bcm_1 = bcm_last_entry.loc[bcm_last_entry.bin_state == 1]

            # We now find the product of the kstar range lists so we can match the
            # merger_type column from the bcm array which tells us what objects
            # merged
            merger_objects_to_track = []
            merger_objects_to_track.extend(
                list(
                    map(
                        lambda x: "{0}{1}".format(
                            str(x[0]).zfill(2), str(x[1]).zfill(2)
                        ),
                        list(itertools.product(kstar1_range, kstar2_range)),
                    )
                )
            )
            merger_objects_to_track.extend(
                list(
                    map(
                        lambda x: "{0}{1}".format(
                            str(x[0]).zfill(2), str(x[1]).zfill(2)
                        ),
                        list(itertools.product(kstar2_range, kstar1_range)),
                    )
                )
            )
            bin_num_save.extend(
                bcm_1.loc[
                    bcm_1.merger_type.isin(merger_objects_to_track)
                ].bin_num.tolist()
            )

            bcm_last_entry = bcm_last_entry.loc[
                bcm_last_entry.bin_num.isin(bin_num_save)
            ]

            # this will tell use the binary state fraction of the systems with a certain final kstar type
            # before we throw out certain binary states if a user requested that.
            bin_state_fraction = bcm_last_entry.groupby("bin_state").tphys.count()
            bin_states = []
            for ii in range(3):
                try:
                    bin_states.append(bin_state_fraction.loc[ii])
                except Exception:
                    bin_states.append(0)
            bin_state_fraction = pd.DataFrame([bin_states], columns=[0, 1, 2])

            bcm = bcm.loc[
                bcm.bin_num.isin(
                    bcm_last_entry.loc[bcm_last_entry.bin_state.isin(use)].bin_num
                )
            ]

    return bcm, bin_state_fraction

def conv_select_singles(bcm_save, bpp_save, final_kstar_1): # fix
    """Select singles"""
    conv_save = bpp_save.loc[
            (bpp_save.kstar_1.isin(final_kstar_1))
        ]
    # select the formation parameters
    conv_save = conv_save.groupby("bin_num").first().reset_index()
    return conv_save


def conv_select(bcm_save, bpp_save, final_kstar_1, final_kstar_2, method, conv_lims):
    """Select bcm data for special convergence cases

    Parameters
    ----------
    bcm_save : `pandas.DataFrame`
        bcm dataframe containing all saved bcm data

    bpp_save : `pandas.DataFrame`
        bpp dataframe containing all saved bpp data

    final_kstar_1 : `list`
        contains list of final primary kstars specified by user

    final_kstar_2 : `list`
        contains list of final primary kstars specified by user

    method : `str`
        stage in binary evolution to check convergence for
        only one method may be supplied and they are specified
        in the inifile

    conv_lims : `dict`
        dictionary where keys are convergence params and the
        values are lists containing a [lo, hi] value to filter the
        convergence param between
        any non-specified convergence params will not be filtered

    Returns
    -------
    conv_save : `pandas.DataFrame`
        filtered dataframe containing binaries that fulfill
        user-specified convergence criteria

    """
    _known_methods = [
        "formation",
        "1_SN",
        "2_SN",
        "disruption",
        "final_state",
        "XRB_form",
    ]

    if method not in _known_methods:
        raise ValueError(
            "You have supplied an "
            "unknown method to filter the "
            "bcm array for convergence. Known methods are "
            "{0}".format(_known_methods)
        )

    if method == "formation":
        # filter the bpp array to find the systems that match the user-specified
        # final kstars
        conv_save = bpp_save.loc[
            ((bpp_save.kstar_1.isin(final_kstar_1)) 
             & (bpp_save.kstar_2.isin(final_kstar_2))
            )
            |
            ((bpp_save.kstar_1.isin(final_kstar_2))
             & (bpp_save.kstar_2.isin(final_kstar_1))
            )
        ]

        # select the formation parameters
        conv_save = conv_save.groupby("bin_num").first().reset_index()

    elif method == "1_SN":
        # select out the systems which will undergo a supernova
        conv_sn_ind = bpp_save.loc[bpp_save.evol_type.isin([15.0, 16.0])].bin_num

        # select out the systems which will produce the user specified final kstars
        # and undergo a supernova
        conv_sn_ind = bpp_save.loc[
            (bpp_save.bin_num.isin(conv_sn_ind))
            & (bpp_save.kstar_1.isin(final_kstar_1))
            & (bpp_save.kstar_2.isin(final_kstar_2))
            & (bpp_save.sep > 0)
        ].bin_num

        # select out the values just before the supernova(e)
        conv_sn = bpp_save.loc[
            (bpp_save.bin_num.isin(conv_sn_ind))
            & (bpp_save.evol_type.isin([15.0, 16.0]))
        ]

        # make sure to select out only the first supernova
        conv_save = conv_sn.groupby("bin_num").first().reset_index()

    elif method == "2_SN":
        # select out the systems which will undergo a supernova
        conv_sn_ind = bpp_save.loc[bpp_save.evol_type.isin([15.0, 16.0])].bin_num

        # select out the systems which will produce the user specified final kstars
        # and undergo a supernova
        conv_sn_ind = bpp_save.loc[
            (bpp_save.bin_num.isin(conv_sn_ind))
            & (bpp_save.kstar_1.isin(final_kstar_1))
            & (bpp_save.kstar_2.isin(final_kstar_2))
            & (bpp_save.sep > 0)
        ].bin_num
        # select out the values just before the supernova(e)
        conv_sn = bpp_save.loc[
            (bpp_save.bin_num.isin(conv_sn_ind))
            & (bpp_save.evol_type.isin([15.0, 16.0]))
        ]

        # select out only the systems that go through 2 supernovae
        conv_sn_2 = conv_sn.loc[conv_sn.groupby("bin_num").size() == 2]

        # make sure to select out only the second supernova
        conv_save = conv_sn_2.groupby("bin_num").nth(1).reset_index()

    elif method == "disruption":
        # filter the bpp array to find the systems that match the user-specified
        # final kstars
        conv_ind = bpp_save.loc[
            (bpp_save.kstar_1.isin(final_kstar_1))
            & (bpp_save.kstar_2.isin(final_kstar_2))
        ].bin_num.unique()

        conv_save = bpp_save.loc[(bpp_save.bin_num.isin(conv_ind))]

        # select out the parameters just before disruption
        # first reset the index:
        conv_save_reset = conv_save.reset_index()

        # next select out the index for the disrupted systems using evol_type == 11
        conv_save_reset_ind = conv_save_reset.loc[
            conv_save_reset.evol_type == 11.0
        ].index

        conv_save = conv_save_reset.iloc[conv_save_reset_ind]

    elif method == "final_state":
        # the bcm array is all that we need!
        conv_save = bcm_save

    elif method == "XRB_form":
        # select out the systems which undergo a SN
        conv_ind = bpp_save.loc[bpp_save.evol_type.isin([15.0, 16.0])].bin_num.unique()
        conv_sn = bpp_save.loc[bpp_save.bin_num.isin(conv_ind)]

        # select out systems when they first enter RLO after the 1st SN
        conv_xrb = conv_sn.loc[
            (conv_sn.kstar_1.isin(final_kstar_1))
            & (conv_sn.kstar_2.isin(final_kstar_2))
            & (conv_sn.RRLO_2 >= 1.0)
            & (conv_sn.sep > 0)
        ]
        conv_save = conv_xrb.groupby("bin_num").first().reset_index()

    if conv_lims:
        for key in conv_lims.keys():
            filter_lo = conv_lims[key][0]
            filter_hi = conv_lims[key][1]
            conv_save_lim = conv_save.loc[conv_save[key] < filter_hi]
            conv_lims_bin_num = conv_save_lim.loc[conv_save[key] > filter_lo].bin_num
    else:
        conv_lims_bin_num = conv_save.bin_num

    return conv_save, conv_lims_bin_num


def pop_write(
    dat_store,
    log_file,
    mass_list,
    number_list,
    bcm,
    bpp,
    initC,
    conv,
    kick_info,
    bin_state_nums,
    match,
    idx,
    **kwargs,
):
    """Writes all the good stuff that you want to save from runFixedPop in a
       single function

    Parameters
    ----------
    dat_store : `pandas HDFStore`
        H5 file to write to

    log_file : `file write`
        log file to write to
    mass_list : `list`
        list containing the mass of the singles, mass of the binaries,
        and mass of the stars

    n_list : `list`
        list containing the number of singles, number of binaries,
        and number of stars

    bcm : `pandas.DataFrame`
        bcm array to write

    bpp : `pandas.DataFrame`
        bpp array to write

    initCond : `pandas.DataFrame`
        initCond array to write

    conv : `pandas.DataFrame`
        conv array to write

    kick_info : `pandas.DataFrame`
        kick_info array to write

    bin_state_nums : `list`
        contains the count of binstates 0,1,2

    match : pandas.DataFrame
        contains the match values for each conv_param

    idx : `int`
        contains the index of the bcm so we can pick up where we left off
        if runFixedPop hits a wall time

    conv_singles : `pandas.DataFrame`
        kwargs conv_singles array to write

    bcm_singles : `pandas.DataFrame`
        kwargs bcm_singles array to write

    bpp_singles : `pandas.DataFrame`
        kwargs bpp_singles array to write

    initC_singles : `pandas.DataFrame`
        kwargs initC_singles array to write

    kick_info_singles : `pandas.DataFrame`
        kwargs kick_info_singles array to write

    Returns
    -------
    Nothing!
    """

    m_keys = ["mass_singles", "mass_binaries", "mass_stars"]
    n_keys = ["n_singles", "n_binaries", "n_stars"]
    for m_write, m_key, n_write, n_key in zip(mass_list, m_keys, number_list, n_keys):
        # save the total_sampled_mass so far
        dat_store.append(m_key, pd.DataFrame([m_write]))
        dat_store.append(n_key, pd.DataFrame([n_write]))
    log_file.write("The total mass sampled so far is: {0}\n".format(mass_list[2]))

    # Save the bcm dataframe
    dat_store.append("bcm", bcm)

    # Save the bpp dataframe
    dat_store.append("bpp", bpp)

    # Save the initial binaries
    # ensure that the index corresponds to bin_num
    dat_store.append("initC", initC.set_index("bin_num", drop=False))

    # Save the converging dataframe
    dat_store.append("conv", conv)

    # Save the converging dataframe
    dat_store.append("kick_info", kick_info)

    # Save number of systems in each bin state
    dat_store.append("bin_state_nums", bin_state_nums)

    # Save the matches
    dat_store.append("match", match)

    # Save the index
    dat_store.append("idx", pd.DataFrame([idx]))

    if "conv_singles" in kwargs.keys():

        # Save the singles conv dataframe
        dat_store.append("conv_singles", kwargs["conv_singles"])

        # Save the singles bcm dataframe
        dat_store.append("bcm_singles", kwargs["bcm_singles"])

        # Save the singles bpp dataframe
        dat_store.append("bpp_singles", kwargs["bpp_singles"])

        # save the singles initCond dataframe
        dat_store.append("initC_singles", kwargs["initC_singles"])

        # save the singles kick_info dataframe      
        dat_store.append("kick_info_singles", kwargs["kick_info_singles"])

    return


def a_from_p(p, m1, m2):
    """Computes the separation from orbital period with KEPLER III

    Parameters
    ----------
    p : float/array
        orbital period [day]
    m1 : float/array
        primary mass [msun]
    m2 : float/array
        secondary mass [msun]

    Returns
    -------
    sep : float/array
        separation [rsun]
    """

    p_yr = p / 365.25
    sep_3 = p_yr ** 2 * (m1 + m2)
    sep = sep_3 ** (1 / 3.0)
    sep_rsun = sep * 215.032
    return sep_rsun


def p_from_a(sep, m1, m2):
    """Computes separation from orbital period with kepler III

    Parameters
    ----------
    sep : float/array
        separation [rsun]
    m1 : float/array
        primary mass [msun]
    m2 : float/array
        secondary mass [msun]

    Returns
    -------
    p : float/array
        orbital period [day]
    """

    sep_au = sep / 215.032
    p_2 = sep_au ** 3 / (m1 + m2)
    p_day = (p_2 ** 0.5) * 365.25
    return p_day


def calc_Roche_radius(M1, M2, A):
    """Get Roche lobe radius (Eggleton 1983)

    Parameters
    ----------
    M1 : float
        Primary mass [any unit]
    M2 : float
        Secondary mass [any unit]
    A : float
        Orbital separation [any unit]

    Returns
    -------
    Roche radius : float
        in units of input 'A'
    """
    q = M1 / M2
    return (
        A
        * 0.49
        * q ** (2.0 / 3.0)
        / (0.6 * q ** (2.0 / 3.0) + np.log(1.0 + q ** (1.0 / 3.0)))
    )


def mass_min_max_select(kstar_1, kstar_2, **kwargs):
    """Select a minimum and maximum mass to filter out binaries in the initial
    parameter sample to reduce the number of unneccessary binaries evolved
    in BSE

    Parameters
    ----------
    kstar_1 : int, list
        BSE stellar type for the primary
        or minimum and maximum stellar types for the primary
    kstar_2 : int, list
        BSE stellar type for the secondary
        or minimum and maximum stellar types for the secondary

    Returns
    -------
    min_mass[0] : float
        minimum primary mass for initial sample
    max_mass[0] : float
        maximum primary mass for initial sample
    min_mass[1] : float
        minimum secondary mass for initial sample
    max_mass[1] : float
        maximum secondary mass for initial sample
    """

    primary_max = kwargs["m_max"] if "m_max" in kwargs.keys() else 150.0
    secondary_max = kwargs["m_max"] if "m_max" in kwargs.keys() else 150.0

    primary_min = kwargs["m1_min"] if "m1_min" in kwargs.keys() else 0.08
    secondary_min = kwargs["m2_min"] if "m2_min" in kwargs.keys() else 0.08

    if ((primary_min < 0.08) | (secondary_min < 0.08)):
        warnings.warn("Tread carefully, BSE is not equipped to handle stellar masses less than 0.08 Msun!")
    if primary_max > 150:
        warnings.warn("Tread carefully, BSE is not equipped to handle stellar masses greater than 150 Msun! And to be honest, we are extrapolating beyond 50 Msun :-/")

    min_mass = [primary_min, secondary_min]
    max_mass = [primary_max, secondary_max]

    if len(kstar_1) == 1:
        # there is a range of final kstar_1s to save
        kstar_1_lo = kstar_1[0]
        kstar_1_hi = kstar_1[0]
    else:
        kstar_1_lo = min(kstar_1)
        kstar_1_hi = max(kstar_1)

    if len(kstar_2) == 1:
        # there is a range of final kstar_1s to save
        kstar_2_lo = kstar_2[0]
        kstar_2_hi = kstar_2[0]
    else:
        kstar_2_lo = min(kstar_2)
        kstar_2_hi = max(kstar_2)

    kstar_lo = [kstar_1_lo, kstar_2_lo]
    kstar_hi = [kstar_1_hi, kstar_2_hi]

    ii = 0
    for k in kstar_lo:
        if k == 14.0:
            min_mass[ii] = 8.0
        elif k == 13.0:
            min_mass[ii] = 3.0
        elif k == 12.0:
            min_mass[ii] = 1.0
        elif k == 11.0:
            min_mass[ii] = 0.8
        elif k == 10.0:
            min_mass[ii] = 0.5
        ii += 1

    ii = 0
    for k in kstar_hi:
        if k == 13.0:
            max_mass[ii] = 60.0
        elif k == 12.0:
            max_mass[ii] = 20.0
        elif k == 11.0:
            max_mass[ii] = 20.0
        elif k == 10.0:
            max_mass[ii] = 20.0
        ii += 1

    return min_mass[0], max_mass[0], min_mass[1], max_mass[1]


def idl_tabulate(x, f, p=5):
    """Function that replicates the IDL int_tabulated function
    which performs a p-point integration on a tabulated set of data

    Parameters
    ----------
    x : array
        tabulated x-value data
    f : array
        tabulated f-value data, same size as x
    p : int
        number of chunks to divide tabulated data into
        Default: 5

    Returns
    -------
    ret : float
        Integration result
    """

    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)

    ret = 0
    for idx in range(0, x.shape[0], p - 1):
        ret += newton_cotes(x[idx: idx + p], f[idx: idx + p])
    return ret


def rndm(a, b, g, size):
    r"""Power-law generator for pdf(x)\propto x^{g} for a<=x<=b

    Parameters
    ----------
    a : float
        Minimum of range for power law
    b : float
        Maximum of range for power law
    g : float
        Index for power law
    size : int
        Number of data points to draw

    Returns
    -------
    power : array
        Array of data sampled from power law distribution with params
        fixed by inputs
    """

    if g == -1:
        raise ValueError("Power law index cannot be exactly -1")
    r = np.random.random(size=size)
    ag, bg = a ** (g + 1), b ** (g + 1)
    return (ag + (bg - ag) * r) ** (1.0 / (g + 1))


def param_transform(dat):
    """Transforms a data set to limits between zero and one
    Leaves some wiggle room on the edges of the data set

    Parameters
    ----------
    dat : array
        array of data to transform between 0 and 1

    Returns
    -------
    datTransformed : array
        array of data with limits between 0 and 1
    """

    datMax = max(dat)
    datMin = min(dat)
    datZeroed = dat - datMin

    datTransformed = datZeroed / ((datMax - datMin))
    if np.max(datTransformed) == 1.0:
        datTransformed[datTransformed == 1.0] = 1 - 1e-6
    if np.min(datTransformed) == 0.0:
        datTransformed[datTransformed == 0.0] = 1e-6
    return datTransformed


def dat_transform(dat, dat_list):
    """Transform a data set to have limits between zero and one using
    param_transform, then transform to log space

    Parameters
    ----------
    dat " DataFrame
        Data to transform to eventually perform KDE
    dat_list : list
        List of DataFrame columns to include in transformation

    Returns
    -------
    dat_trans : array
        Transformed data for columns in dat_list
    """

    dat_trans = []
    for column in dat_list:
        dat_trans.append(ss.logit(param_transform(dat[column])))
    dat_trans = np.vstack([dat_trans])

    return dat_trans


def dat_un_transform(dat_sample, dat_set, dat_list):
    """Un-transform data that was transformed in dat_transform

    Parameters
    ----------
    dat_sample : array
        Data sampled from kde generated with transformed data
    dat_set : DataFrame
        Un-transformed data (same as dat in dat_transform)
    dat_list : list
        List of DataFrame columns to include in transformation

    Returns
    -------
    dat : array
        Array of data sampled from kde that is transformed back to
        bounds of the un-transformed data set the kde is generated from
    """
    dat = []

    dat_exp = ss.expit(dat_sample)
    for ii, column in zip(range(len(dat_list)), dat_list):
        dat_untrans = dat_exp[ii, :] * (
            max(dat_set[column]) - min(dat_set[column])
        ) + min(dat_set[column])
        dat.append(dat_untrans)
    dat = np.vstack(dat)
    return dat


def knuth_bw_selector(dat_list):
    """Selects the kde bandwidth using Knuth's rule implemented in Astropy
    If Knuth's rule raises error, Scott's rule is used

    Parameters
    ----------
    dat_list : list
        List of data arrays that will be used to generate a kde

    Returns
    -------
    bw_min : float
        Minimum of bandwidths for all of the data arrays in dat_list
    """

    bw_list = []
    for dat in dat_list:
        try:
            bw = astrostats.knuth_bin_width(dat)
        except Exception:
            print("Using Scott Rule!!")
            bw = astrostats.scott_bin_width(dat)
        bw_list.append(bw)
    return np.mean(bw_list)


def get_Z_from_FeH(FeH, Z_sun=0.02):
    """
    Converts from FeH to Z under the assumption that
    all stars have the same abundance as the sun
    
    Parameters
    ----------
    FeH : array
        Fe/H values to convert
    Z_sun : float
        solar metallicity
    
    Returns
    -------
    Z : array
        metallicities corresponding to Fe/H
    """
    Z = 10**(FeH + np.log10(Z_sun))
    return Z


def get_FeH_from_Z(Z, Z_sun=0.02):
    """
    Converts from Z to FeH under the assumption that
    all stars have the same abundance as the sun
    
    Parameters
    ----------
    Z : array
        metallicities to convert to Fe/H
    Z_sun : float
        solar metallicity
    
    Returns
    -------
    FeH : array
        Fe/H corresponding to metallicities
    """
    FeH = np.log10(Z) - np.log10(Z_sun)
    return FeH


def get_binfrac_of_Z(Z):
    '''
    Calculates the theoretical binary fraction as a 
    function of metallicity. Following Moe+2019
    
    Parameters
    ----------
    Z : array
        metallicity Z values
    
    Returns
    -------
    binfrac : array
        binary fraction values
    '''
    FeH = get_FeH_from_Z(Z)
    FeH_low = FeH[np.where(FeH<=-1.0)]
    FeH_high = FeH[np.where(FeH>-1.0)]
    binfrac_low = -0.0648 * FeH_low + 0.3356
    binfrac_high = -0.1977 * FeH_high + 0.2025
    binfrac = np.append(binfrac_low, binfrac_high)
    return binfrac


def get_porb_norm(Z, close_logP=4.0, wide_logP=6.0, binfrac_tot_solar=0.66, Z_sun=0.02):
    '''Returns normalization constants to produce log normals consistent with Fig 19 of Moe+19
    for the orbital period distribution
                
    Parameters
    ----------
    Z : array
        metallicity values
    close_logP : float
        divding line beween close and intermediate orbits
    wide_logP : float
        dividing line between intermediate and wide orbits
    binfrac_tot : float
        integrated total binary fraction at solar metallicity
                
    Returns
    -------
    norm_wide : float
        normalization factor for kde for wide binaries
    norm_close : float
        normalization factor for kde for wide binaries
    '''
    from scipy.stats import norm
    from scipy.integrate import trapezoid
    from scipy.interpolate import interp1d
    
    # fix to values used in Moe+19
    logP_lo_lim=0
    logP_hi_lim=9
    log_P = np.linspace(logP_lo_lim, logP_hi_lim, 10000)
    
    logP_pdf = norm.pdf(log_P, loc=4.9, scale=2.3)
    
    # set up the wide binary fraction inflection point
    norm_wide = binfrac_tot_solar/trapezoid(logP_pdf, log_P)
    
    # set up the close binary fraction inflection point
    FeHclose = np.linspace(-3.0, 0.5, 100)
    fclose = -0.0648 * FeHclose + 0.3356
    fclose[FeHclose > -1.0] = -0.1977 * FeHclose[FeHclose > -1.0] + 0.2025
    Zclose = get_Z_from_FeH(FeHclose, Z_sun=Z_sun)
    
    fclose_interp = interp1d(Zclose, fclose)
    
    fclose_Z = fclose_interp(Z)
    norm_close = fclose_Z/trapezoid(logP_pdf[log_P < close_logP], log_P[log_P < close_logP])
    
    return norm_wide, norm_close


def get_met_dep_binfrac(met):
    '''Returns a population-wide binary fraction consistent with
    Moe+19 based on the supplied metallicity

    Parameters
    ----------
    met : float
        metallicity of the population

    Returns
    -------
    binfrac : float
        binary fraction of the population based on metallicity

    '''
    logP_hi_lim = 9
    logP_lo_lim = 0
    wide_logP = 6
    close_logP = 4
    neval = 5000

    from scipy.interpolate import interp1d
    from scipy.integrate import trapezoid
    from scipy.stats import norm

    norm_wide, norm_close = get_porb_norm(met)
    prob_wide = norm.pdf(np.linspace(wide_logP, logP_hi_lim, neval), loc=4.9, scale=2.3)*norm_wide
    prob_close = norm.pdf(np.linspace(logP_lo_lim, close_logP, neval), loc=4.9, scale=2.3)*norm_close
    slope = -(prob_close[-1] - prob_wide[0]) / (wide_logP - close_logP)
    prob_intermediate = slope * (np.linspace(close_logP, wide_logP, neval) - close_logP) + prob_close[-1]
    prob_interp_int = interp1d(np.linspace(close_logP, wide_logP, neval), prob_intermediate)

    x_dat = np.hstack([np.linspace(logP_lo_lim, close_logP, neval),
                       np.linspace(close_logP, wide_logP, neval),
                       np.linspace(wide_logP, logP_hi_lim, neval),])
    y_dat = np.hstack([prob_close, prob_interp_int(np.linspace(close_logP, wide_logP, neval)), prob_wide])

    binfrac = trapezoid(y_dat, x_dat)/0.66 * 0.5

    return float(np.round(binfrac, 2))

def error_check(BSEDict, filters=None, convergence=None, sampling=None):
    """Checks that values in BSEDict, filters, and convergence are viable"""
    if not isinstance(BSEDict, dict):
        raise ValueError("BSE flags must be supplied via a dictionary")

    if filters is not None:
        if not isinstance(filters, dict):
            raise ValueError("Filters criteria must be supplied via a dictionary")
        for option in ["binary_state", "timestep_conditions"]:
            if option not in filters.keys():
                raise ValueError(
                    "Inifile section filters must have option {0} supplied".format(
                        option
                    )
                )

    if convergence is not None:
        if not isinstance(convergence, dict):
            raise ValueError("Convergence criteria must be supplied via a dictionary")
        for option in [
            "pop_select",
            "convergence_params",
            "convergence_limits",
            "match",
            "apply_convergence_limits",
        ]:
            if option not in convergence.keys():
                raise ValueError(
                    "Inifile section convergence must have option {0} supplied".format(
                        option
                    )
                )

    if sampling is not None:
        if not isinstance(sampling, dict):
            raise ValueError("Sampling criteria must be supplied via a dictionary")
        for option in ["sampling_method", "SF_start", "SF_duration", "metallicity", "keep_singles"]:
            if option not in sampling.keys():
                raise ValueError(
                    "Inifile section sampling must have option {0} supplied".format(
                        option
                    )
                )
        if ("qmin" not in sampling.keys()) & ("m2_min" not in sampling.keys()) & (sampling["sampling_method"] == 'independent'):
            raise ValueError("You have not specified qmin or m2_min. At least one of these must be specified.")
    # filters
    if filters is not None:
        flag = "binary_state"
        if any(x not in [0, 1, 2] for x in filters[flag]):
            raise ValueError(
                "{0} needs to be a subset of [0,1,2] (you set it to {1})".format(
                    flag, filters[flag]
                )
            )
        flag = "timestep_conditions"
        if (type(filters[flag]) != str) and (type(filters[flag]) != list):
            raise ValueError(
                "{0} needs to either be a string like 'dtp=None' or a list of conditions like [['binstate==0', 'dtp=1.0']] (you set it to {1})".format(
                    flag, filters[flag]
                )
            )

    # convergence
    if convergence is not None:
        flag = "convergence_limits"
        if convergence[flag]:
            for item, key in zip(convergence.items(), convergence.keys()):
                if len(item) != 2:
                    raise ValueError(
                        "The value for key '{0:s}' needs to be a list of length 2, it is length: {1:i}".format(
                            key, len(item)
                        )
                    )
        flag = "pop_select"
        if not convergence[flag] in [
            "formation",
            "1_SN",
            "2_SN",
            "disruption",
            "final_state",
            "XRB_form",
        ]:
            raise ValueError(
                "{0} needs to be in the list: ['formation', '1_SN', '2_SN', 'disruption', 'final_state', 'XRB_form'] "
                "(you set it to {1})".format(
                                             flag, convergence[flag]
                                            )
            )

        flag = "match"
        if not isinstance(convergence[flag], float):
            raise ValueError(
                "{0} must be a float (you set it to {1})".format(
                    flag, convergence[flag]
                )
            )

        flag = "convergence_params"
        acceptable_convergence_params = [
            "mass_1",
            "mass_2",
            "sep",
            "porb",
            "ecc",
            "massc_1",
            "massc_2",
            "rad_1",
            "rad_2",
        ]
        for param in convergence[flag]:
            if param not in acceptable_convergence_params:
                raise ValueError(
                    "Supplied convergence parameter {0} is not in list of "
                    "acceptable convergence parameters {1}".format(
                        param, acceptable_convergence_params
                    )
                )

        flag = "convergence_limits"
        if type(convergence[flag]) != dict:
            raise ValueError(
                "Supplied convergence limits must be passed as a dict "
                "(you passed type {0})".format(type(convergence[flag]))
                )

        for key in convergence[flag].keys():
            if key not in convergence["convergence_params"]:
                raise ValueError(
                    "Supplied convergence limits must correspond to already "
                    "supplied convergence_params. The supplied convergence_params "
                    "are {0}, while you supplied {1}".format(
                        convergence["convergence_params"], key)
                    )
        flag = "apply_convergence_limits"
        if type(convergence[flag]) != bool:
            raise ValueError(
                "apply_convergence_limits must be either True or False, "
                "you supplied {}".format(
                    convergence[flag])
                )  
 
    # sampling
    if sampling is not None:
        flag = "sampling_method"
        acceptable_sampling = ["multidim", "independent"]
        if sampling[flag] not in acceptable_sampling:
            raise ValueError(
                "sampling_method must be one of {0} you supplied {1}.".format(
                    acceptable_sampling, sampling[flag]
                )
            )

        flag = "metallicity"
        if not isinstance(sampling[flag], float):
            raise ValueError(
                "{0} must be a float (you set it to {1})".format(flag, sampling[flag])
            )
        if sampling[flag] <= 0:
            raise ValueError(
                "{0} needs to be greater than or equal to 0 (you set it to {1})".format(
                    flag, sampling[flag]
                )
            )

    # use the cosmic-settings.json file to define the valid ranges for BSE flags
    settings_path = io_files("cosmic.data").joinpath('cosmic-settings.json')
    settings = json.loads(settings_path.read_text(encoding='utf-8'))

    handle_separately = ['qcrit_array', 'natal_kick_array', 'fprimc_array']

    # go through the different categories in the settings file
    for cat in settings:
        # ignore anything that's not BSE
        if cat['category'] != "bse":
            continue

        # go through each flag in the settings
        for flag in cat['settings']:
            # if the user has provided it
            if flag['name'] in BSEDict and flag['name'] not in handle_separately:
                user_val = BSEDict[flag['name']]

                # track the valid options and whether the flag has matched any of them
                options = [o["name"] for o in flag["options"]]
                flag_is_valid = False

                # check each option
                for opt in options:
                    # for strings, we do something more complex
                    if isinstance(opt, str):
                        if opt == "positive values" and user_val > 0:
                            flag_is_valid = True
                            break
                        elif opt == "negative values" and user_val < 0:
                            flag_is_valid = True
                            break
                        # for things of the form "range [a,b)"
                        elif opt.startswith("range"):
                            # strip to just the brackets
                            r = opt[5:].strip()

                            # get the brackets and ensure the format is correct
                            start_brac, end_brac = r[0], r[-1]
                            if start_brac not in ["[", "("] or end_brac not in ["]", ")"] or ',' not in r:
                                raise ValueError(
                                    f"Range option for {flag['name']} is not formatted correctly."
                                )
                            # get the range value and check if the user value is in range
                            r_lo, r_hi = map(float, r[1:-1].split(","))
                            lower_ok = user_val > r_lo if start_brac == "(" else user_val >= r_lo
                            upper_ok = user_val < r_hi if end_brac == ")" else user_val <= r_hi
                            if lower_ok and upper_ok:
                                flag_is_valid = True
                                break
                    # otherwise, just do a direct comparison
                    elif user_val == opt:
                        flag_is_valid = True
                        break

                # if we didn't find a match, raise an error
                if not flag_is_valid:
                    raise ValueError(
                        f"{flag['name']} must be one of {options} (you set it to '{user_val}')"
                    )

    if "dtp" in BSEDict.keys():
        if BSEDict["dtp"] < 0:
            raise ValueError(
                f"dtp needs to be greater than or equal to 0 (you set it to '{BSEDict['dtp']:0.2f}')"
            )

    if "kickflag" in BSEDict.keys():
        if BSEDict["kickflag"] in [-1, -2] and ((BSEDict['ecsn'] != 2.25) or (BSEDict['ecsn_mlow'] != 1.6)):
            warnings.warn("You have chosen a kick flag that assumes compact object formation "
                            "according to Giacobbo & Mapelli 2020, but supplied electron "
                            "capture SN (ECSN) flags that are inconsistent with this study. "
                            "To maintain consistency, COSMIC will update your "
                            "ECSN flags to be ecsn=2.25 and ecsn_mlow=1.6")
            BSEDict['ecsn'] = 2.25
            BSEDict['ecsn_mlow'] = 1.6

    if "ecsn_mlow" in BSEDict.keys() and "ecsn" in BSEDict.keys():
        if BSEDict["ecsn_mlow"] > BSEDict["ecsn"]:
            raise ValueError(
                f"`ecsn_mlow` needs to be less than `ecsn`, (you set `ecsn_mlow` to {BSEDict['ecsn_mlow']} "
                f"and `ecsn` to {BSEDict['ecsn']})"
            )
    
    # ensure the natal kick array is the correct shape and each value is in the valid range
    if "natal_kick_array" in BSEDict.keys():
        shape = np.array(BSEDict["natal_kick_array"]).shape
        if shape != (2, 5):
            raise ValueError(
                f"'natal_kick_array' must have shape (2,5) (you supplied list, or array with shape '{shape}')"
            )
        
        valid_ranges = [
            (0, np.inf),        # velocity magnitude
            (-90, 90),          # polar angle
            (0, 360),           # azimuthal angle
            (0, 360),           # mean anomaly
            (-np.inf, np.inf)   # random seed
        ]

        for i in range(2):
            for j in range(5):
                val = BSEDict["natal_kick_array"][i][j]
                low, high = valid_ranges[j]
                if not (low <= val <= high) and val != -100.0:
                    raise ValueError(
                        f"Value at position ({i},{j}) in 'natal_kick_array' must be in range [{low}, {high}] "
                        f"(you set it to '{val}')"
                    )


    if "fprimc_array" in BSEDict.keys():
        if np.any(np.array(BSEDict["fprimc_array"]) < 0.0) or len(BSEDict["fprimc_array"]) != 16:
            raise ValueError(
                f"fprimc_array values must be >= 0 and there must be 16 values "
                f'(you set them to {BSEDict["fprimc_array"]}], length={len(BSEDict["fprimc_array"])})'
            )
        
    if "qcrit_array" in BSEDict.keys():
        if np.any(np.array(BSEDict["qcrit_array"]) < 0.0) or len(BSEDict["qcrit_array"]) != 16:
            raise ValueError(
                f"qcrit_array values must be >= 0 and there must be 16 values "
                f'(you set them to {BSEDict["qcrit_array"]}], length={len(BSEDict["qcrit_array"])})'
            )

    return


def check_initial_conditions(full_initial_binary_table):
    """Checks initial conditions and reports warnings

    Only warning provided right now is if star begins in Roche lobe
    overflow
    """

    def rzamsf(m):
        """A function to evaluate Rzams
        ( from Tout et al., 1996, MNRAS, 281, 257 ).
        """
        mx = np.sqrt(m)
        rzams = (
            (a[7] * m ** 2 + a[8] * m ** 6) * mx
            + a[9] * m ** 11
            + (a[10] + a[11] * mx) * m ** 19
        ) / (a[12] + a[13] * m ** 2 + (a[14] * m ** 8 + m ** 18 + a[15] * m ** 19) * mx)

        return rzams

    no_singles = ((full_initial_binary_table["mass_1"] > 0.0)
                  & (full_initial_binary_table["mass_2"] > 0.0)
                  & (full_initial_binary_table["porb"] > 0.0))
    initial_binary_table = full_initial_binary_table[no_singles]

    z = np.asarray(initial_binary_table["metallicity"])
    zpars, a = zcnsts(z)

    mass1 = np.asarray(initial_binary_table["mass_1"])
    mass2 = np.asarray(initial_binary_table["mass_2"])

    if np.all(mass2 == 0.0):
        return
    else:
        rzams1 = rzamsf(mass1)
        rzams2 = rzamsf(mass2)

        # assume some time step in order to calculate sep
        yeardy = 365.24
        aursun = 214.95
        tb = np.asarray(initial_binary_table["porb"]) / yeardy
        sep = aursun * (tb * tb * (mass1 + mass2)) ** (1.0 / 3.0)

        rol1 = calc_Roche_radius(mass1, mass2, sep)
        rol2 = calc_Roche_radius(mass2, mass1, sep)

        # check for a ZAMS that starts in RFOL
        mask = ((np.array(initial_binary_table["kstar_1"]) == 1) & (rzams1 >= rol1)) | (
            (initial_binary_table["kstar_2"] == 1) & (rzams2 >= rol2)
        )
        if mask.any():
            warnings.warn(
                "At least one of your initial binaries is starting in Roche Lobe Overflow:\n{0}".format(
                    initial_binary_table[mask]
                )
            )

        return


def convert_kstar_evol_type(bpp):
    """Provides way to convert integer values to their string counterpart

    The underlying fortran code relies on integers to indicate
    things like the evoltuionary stage of the star as well as
    key moments in its evolutionary track. If you pass the
    data frame returned from running

        ```Evolve.evolve```

    you can convert the columns with these integer proxies
    to their true astrophysical meaning.
    """
    kstar_int_to_string_dict = {
        0: "Main Sequence (MS), < 0.7 M⊙",
        1: "MS, > 0.7 M⊙",
        2: "Hertzsprung Gap",
        3: "First Giant Branch",
        4: "Core Helium Burning",
        5: "Early Asymptotic Giant Branch (AGB)",
        6: "Thermally Pulsing AGB",
        7: "Naked Helium Star MS",
        8: "Naked Helium Star Hertzsprung Gap",
        9: "Naked Helium Star Giant Branch",
        10: "Helium White Dwarf",
        11: "Carbon/Oxygen White Dwarf",
        12: "Oxygen/Neon White Dwarf",
        13: "Neutron Star",
        14: "Black Hole",
        15: "Massless Remnant",
    }

    kstar_string_to_int_dict = {v: k for k, v in kstar_int_to_string_dict.items()}

    evolve_type_int_to_string_dict = {
        1: "initial state",
        2: "kstar change",
        3: "begin Roche lobe overflow",
        4: "end Roche lobe overlow",
        5: "contact",
        6: "coalescence",
        7: "begin common envelope",
        8: "end common envelope",
        9: "no remnant leftover",
        10: "max evolution time",
        11: "binary disruption",
        12: "begin symbiotic phase",
        13: "end symbiotic phase",
        14: "blue straggler",
        15: "supernova of primary",
        16: "supernova of secondary",
       100: "RLOF interpolation timeout error"
    }

    evolve_type_string_to_int_dict = {
        v: k for k, v in evolve_type_int_to_string_dict.items()
    }

    if bpp.kstar_1.dtype in [int, float]:
        # convert from integer to string
        bpp["kstar_1"] = bpp["kstar_1"].astype(int)
        bpp["kstar_1"] = bpp["kstar_1"].apply(lambda x: kstar_int_to_string_dict[x])
    else:
        # convert from string to integer
        bpp["kstar_1"] = bpp["kstar_1"].apply(lambda x: kstar_string_to_int_dict[x])

    if bpp.kstar_2.dtype in [int, float]:
        # convert from integer to string
        bpp["kstar_2"] = bpp["kstar_2"].astype(int)
        bpp["kstar_2"] = bpp["kstar_2"].apply(lambda x: kstar_int_to_string_dict[x])
    else:
        # convert from string to integer
        bpp["kstar_2"] = bpp["kstar_2"].apply(lambda x: kstar_string_to_int_dict[x])

    if bpp.evol_type.dtype in [int, float]:
        # convert from integer to string
        bpp["evol_type"] = bpp["evol_type"].astype(int)
        bpp["evol_type"] = bpp["evol_type"].apply(
            lambda x: evolve_type_int_to_string_dict[x]
        )
    else:
        # convert from string to integer
        bpp["evol_type"] = bpp["evol_type"].apply(
            lambda x: evolve_type_string_to_int_dict[x]
        )

    return bpp


def parse_inifile(inifile):
    """Provides a method for parsing the inifile and returning dicts of each section"""
    if inifile is None:
        raise ValueError("Please supply an inifile")
    elif not os.path.isfile(inifile):
        raise ValueError("inifile supplied does not exist")

    binOps = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Mod: operator.mod,
    }

    def arithmetic_eval(s):
        """Allows us to control how the strings from the inifile get parses"""
        node = ast.parse(s, mode="eval")

        def _eval(node):
            """Different strings receive different evaluation"""
            if isinstance(node, ast.Expression):
                return _eval(node.body)
            elif isinstance(node, ast.Str):
                return node.s
            elif isinstance(node, ast.Num):
                return node.n
            elif isinstance(node, ast.BinOp):
                return binOps[type(node.op)](_eval(node.left), _eval(node.right))
            elif isinstance(node, ast.List):
                return [_eval(x) for x in node.elts]
            elif isinstance(node, ast.Name):
                result = VariableKey(item=node)
                constants_lookup = {
                    "True": True,
                    "False": False,
                    "None": None,
                }
                value = constants_lookup.get(
                    result.name,
                    result,
                )
                if type(value) == VariableKey:
                    # return regular string
                    return value.name
                else:
                    # return special string like True or False
                    return value
            elif isinstance(node, ast.NameConstant):
                # None, True, False are nameconstants in python3, but names in 2
                return node.value
            else:
                raise Exception("Unsupported type {}".format(node))

        return _eval(node.body)

    # ---- Create configuration-file-parser object and read parameters file.
    cp = ConfigParser()
    cp.optionxform = str
    cp.read(inifile)

    # ---- Read needed variables from the inifile
    dictionary = {}
    for section in cp.sections():
        # for cosmic we skip any CMC stuff
        # if "cmc" is a section in the ini file, then we can optionally skip the 
        # COSMIC population sections (or not, if they exist)
        if section == "cmc":
            if "rand_seed" not in dictionary.keys():
                dictionary["rand_seed"] = {}
                dictionary["rand_seed"]["seed"] = 0
            if "filters" not in dictionary.keys():
                dictionary["filters"] = 0
            if "convergence" not in dictionary.keys():
                dictionary["convergence"] = 0
            if "sampling" not in dictionary.keys():
                dictionary["sampling"] = 0
            continue
        dictionary[section] = {}
        for option in cp.options(section):
            opt = cp.get(section, option)
            if "\n" in opt:
                raise ValueError("We have detected an error in your inifile. A parameter was read in with the following "
                                 "value: {0}. Likely, you have an unexpected syntax, such as a space before an parameter/option (i.e. "
                                 "the parameter must be flush to the far left of the file".format(opt))
            try:
                dictionary[section][option] = arithmetic_eval(opt)
            except Exception:
                dictionary[section][option] = json.loads(opt)
            finally:
                if option not in dictionary[section].keys():
                    raise ValueError("We have detected an error in your inifile. The folloiwng parameter failed to be read correctly: {0}".format(option))

    BSEDict = dictionary["bse"]
    seed_int = int(dictionary["rand_seed"]["seed"])
    filters = dictionary["filters"]
    convergence = dictionary["convergence"]
    sampling = dictionary["sampling"]

    return BSEDict, seed_int, filters, convergence, sampling


def explain_setting(setting):
    """Provides explanation for a BSE setting from the cosmic-settings.json file

    Parameters
    ----------
    setting : str
        Name of BSE setting to explain
    """
    # use the cosmic-settings.json file to define the valid ranges for BSE flags
    settings_path = io_files("cosmic.data").joinpath('cosmic-settings.json')
    settings = json.loads(settings_path.read_text(encoding='utf-8'))

    strip_tags = lambda s: re.sub(r'<[^>]+>', '', s).replace("&amp;", "&")

    BOLD = '\033[1m'
    GREEN = '\033[92m'
    END = '\033[0m'

    for cat in settings:
        # ignore anything that's not BSE
        if cat['category'] != "bse":
            continue

        # go through each flag in the settings
        for flag in cat['settings']:
            if flag['name'] == setting:
                print(f"\n{BOLD}{flag['name']}{END}")
                print("-" * len(flag['name']))
                print(f"{strip_tags(flag['description'])}")
                print("\nValid options (default marked in green and with *):")
                for opt in flag['options']:
                    print(f"  {f'{GREEN}*' if 'default' in opt else '-'} {opt['name']}: {strip_tags(opt['description'])}{END}")
                return
            
    raise ValueError(f"Unknown setting '{setting}'")


class VariableKey(object):
    """
    A dictionary key which is a variable.
    @ivar item: The variable AST object.
    """

    def __init__(self, item):
        self.name = item.id

    def __eq__(self, compare):
        return compare.__class__ == self.__class__ and compare.name == self.name

    def __hash__(self):
        return hash(self.name)
