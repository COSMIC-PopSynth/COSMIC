# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017 - 2019)
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
import scipy.integrate
import numpy as np
import scipy.special as ss
import astropy.stats as astrostats
import warnings
import ast
import operator
import json

from configparser import ConfigParser
from .bse_utils.zcnsts import zcnsts

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = ['Scott Coughlin <scott.coughlin@ligo.org>',
               'Michael Zevin <zevin@northwestern.edu>']
__all__ = ['filter_bpp_bcm', 'conv_select', 'mass_min_max_select',
           'idl_tabulate', 'rndm', 'param_transform', 'dat_transform',
           'dat_un_transform', 'knuth_bw_selector', 'error_check',
           'check_initial_conditions', 'convert_kstar_evol_type']

def filter_bpp_bcm(bcm, bpp, method, kstar1_range, kstar2_range):
    """Filter the output of bpp and bcm

    Parameters
    ----------
    bcm : `pandas.DataFrame`
        bcm dataframe

    bpp : `pandas.DataFrame`
        bpp dataframe

    method : `dict`,
        one or more methods by which to filter the
        bpp or bcm table, e.g. ``{'select_final_state' : False}``;
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
    _known_methods = ['select_final_state',
                      'binary_state']

    if not set(method.keys()).issubset(set(_known_methods)):
        raise ValueError("You have supplied an "
                         "unknown method to filter out "
                         "the bpp or bcm array. Known methods are "
                         "{0}".format(_known_methods))

    for meth, use in method.items():
        if (meth == 'select_final_state') and use:
            bcm = bcm.iloc[bcm.reset_index().groupby('bin_num').tphys.idxmax()]
        elif (meth == 'binary_state'):
            bin_num_save = []

            # in order to find the properities of disrupted or systems
            # that are alive today we can simply check the last entry in the bcm
            # array for the system and see what its properities are today
            bcm_0_2 = bcm.loc[(bcm.bin_state != 1)]
            bin_num_save.extend(bcm_0_2.loc[(bcm_0_2.kstar_1.isin(kstar1_range)) &
                                          (bcm_0_2.kstar_2.isin(kstar2_range))].bin_num.tolist())
            # in order to find the properities of merged systems
            # we actually need to search in the BPP array for the properities
            # of the objects right at merge because the bcm will report
            # the post merge object only
            bcm_1 = bcm.loc[bcm.bin_state == 1]
            bpp_1 = bpp.loc[bpp.bin_num.isin(bcm_1.bin_num)]
            bin_num_save.extend(bpp_1.loc[(bpp_1.kstar_1.isin(kstar1_range)) &
                                          (bpp_1.kstar_2.isin(kstar2_range)) &
                                          (bpp_1.evol_type == 3)].bin_num.tolist())

            bcm = bcm.loc[bcm.bin_num.isin(bin_num_save)]

            # this will tell use the binary state fraction of the systems with a certain final kstar type
            # before we throw out certain binary states if a user requested that.
            bin_state_fraction = bcm.groupby('bin_state').tphys.count()

            bcm = bcm.loc[bcm.bin_state.isin(use)]

    return bcm, bin_state_fraction

def conv_select(bcm_save_tot, bcm_save_last, bpp_save_tot, final_kstar_1, final_kstar_2, method):
    """Select bcm data for special convergence cases

    Parameters
    ----------
    bcm_save_tot : `pandas.DataFrame`
        bcm dataframe containing all saved bcm data

    bcm_save_last : `pandas.DataFrame`
        bcm dataframe containing bcm data from last
        iteration

    bpp_save_tot : `pandas.DataFrame`
        bpp dataframe containing all saved bpp data

    final_kstar_1 : `list`
        contains list of final primary kstars specified by user

    final_kstar_2 : `list`
        contains list of final primary kstars specified by user

    method : `dict`,
        one or more methods by which to filter the
        bcm table, e.g. ``{'formation' : True}``;
        This means you want to only compute the convergence
        of the population at formation, e.g. BH-BH formation

    Returns
    -------
    bcm_conv_tot : `pandas.DataFrame`
        filtered bcm dataframe containing all saved bcm
        data

    bcm_conv_last : `pandas.DataFrame`
        filtered bcm dataframe containing saved bcm
        data from last iteration

    """
    _known_methods = ['formation', '1_SN', '2_SN', 'disruption', 'final_state', 'XRB_form']

    if not method in _known_methods:
        raise ValueError("You have supplied an "
                         "unknown method to filter the "
                         "bcm array for convergence. Known methods are "
                         "{0}".format(_known_methods))

    if len(bcm_save_tot) == len(bcm_save_last):
        bcm_save_last = bcm_save_tot[:int(len(bcm_save_tot)/2)]
    bpp_save_last = bpp_save_tot.loc[bpp_save_tot.bin_num.isin(bcm_save_last.bin_num)]

    if method == 'formation':
        # filter the bpp array to find the systems that match the user-specified 
        # final kstars
        conv_tot = bpp_save_tot.loc[(bpp_save_tot.kstar_1.isin(final_kstar_1)) &\
                                    (bpp_save_tot.kstar_2.isin(final_kstar_2)) &\
                                    (bpp_save_tot.sep > 0)]
        conv_last = bpp_save_last.loc[(bpp_save_last.kstar_1.isin(final_kstar_1)) &\
                                      (bpp_save_last.kstar_2.isin(final_kstar_2)) &\
                                      (bpp_save_last.sep > 0)]

        # select the formation parameters
        conv_tot = conv_tot.groupby('bin_num').first().reset_index()
        conv_last = conv_last.groupby('bin_num').first().reset_index()

    elif method == '1_SN':
        # select out the systems which will undergo a supernova
        conv_tot_sn_ind = bpp_save_tot.loc[bpp_save_tot.evol_type.isin([15.0, 16.0])].bin_num
        conv_last_sn_ind = bpp_save_last.loc[bpp_save_last.evol_type.isin([15.0, 16.0])].bin_num

        # select out the systems which will produce the user specified final kstars
        # and undergo a supernova
        conv_tot_sn_ind = bpp_save_tot.loc[(bpp_save_tot.bin_num.isin(conv_tot_sn_ind)) &\
                                           (bpp_save_tot.kstar_1.isin(final_kstar_1)) &\
                                           (bpp_save_tot.kstar_2.isin(final_kstar_2)) &\
                                           (bpp_save_tot.sep > 0)].bin_num
        conv_last_sn_ind = bpp_save_last.loc[(bpp_save_last.bin_num.isin(conv_last_sn_ind)) &\
                                             (bpp_save_last.kstar_1.isin(final_kstar_1)) &\
                                             (bpp_save_last.kstar_2.isin(final_kstar_2)) &\
                                             (bpp_save_last.sep > 0)].bin_num            
        # select out the values just before the supernova(e)
        conv_tot_sn = bpp_save_tot.loc[(bpp_save_tot.bin_num.isin(conv_tot_sn_ind)) &\
                                       (bpp_save_tot.evol_type.isin([15.0, 16.0]))]
        conv_last_sn = bpp_save_last.loc[(bpp_save_last.bin_num.isin(conv_last_sn_ind)) &\
                                         (bpp_save_last.evol_type.isin([15.0, 16.0]))]

        # make sure to select out only the first supernova
        conv_tot = conv_tot_sn.groupby('bin_num').first().reset_index()
        conv_last = conv_last_sn.groupby('bin_num').first().reset_index()

    elif method == '2_SN':
        # select out the systems which will undergo a supernova
        conv_tot_sn_ind = bpp_save_tot.loc[bpp_save_tot.evol_type.isin([15.0, 16.0])].bin_num
        conv_last_sn_ind = bpp_save_last.loc[bpp_save_last.evol_type.isin([15.0, 16.0])].bin_num

        # select out the systems which will produce the user specified final kstars
        # and undergo a supernova
        conv_tot_sn_ind = bpp_save_tot.loc[(bpp_save_tot.bin_num.isin(conv_tot_sn_ind)) &\
                                           (bpp_save_tot.kstar_1.isin(final_kstar_1)) &\
                                           (bpp_save_tot.kstar_2.isin(final_kstar_2)) &\
                                           (bpp_save_tot.sep > 0)].bin_num
        conv_last_sn_ind = bpp_save_last.loc[(bpp_save_last.bin_num.isin(conv_last_sn_ind)) &\
                                             (bpp_save_last.kstar_1.isin(final_kstar_1)) &\
                                             (bpp_save_last.kstar_2.isin(final_kstar_2)) &\
                                             (bpp_save_last.sep > 0)].bin_num 
        # select out the values just before the supernova(e)
        conv_tot_sn = bpp_save_tot.loc[(bpp_save_tot.bin_num.isin(conv_tot_sn_ind)) &\
                                       (bpp_save_tot.evol_type.isin([15.0, 16.0]))]
        conv_last_sn = bpp_save_last.loc[(bpp_save_last.bin_num.isin(conv_last_sn_ind)) &\
                                         (bpp_save_last.evol_type.isin([15.0, 16.0]))]         

        # select out only the systems that go through 2 supernovae
        conv_tot_sn_2 = conv_tot_sn.loc[conv_tot_sn.groupby('bin_num').size() == 2]
        conv_last_sn_2 = conv_last_sn.loc[conv_last_sn.groupby('bin_num').size() == 2]

        # make sure to select out only the second supernova
        conv_tot = conv_tot_sn_2.groupby('bin_num').nth(1).reset_index()
        conv_last = conv_last_sn_2.groupby('bin_num').nth(1).reset_index()

    elif method == 'disruption':
        # filter the bpp array to find the systems that match the user-specified 
        # final kstars
        conv_tot_ind = bpp_save_tot.loc[(bpp_save_tot.kstar_1.isin(final_kstar_1)) &\
                                        (bpp_save_tot.kstar_2.isin(final_kstar_2))].bin_num.unique()
        conv_last_ind = bpp_save_last.loc[(bpp_save_last.kstar_1.isin(final_kstar_1)) &\
                                          (bpp_save_last.kstar_2.isin(final_kstar_2))].bin_num.unique()

        conv_tot = bpp_save_tot.loc[(bpp_save_tot.bin_num.isin(conv_tot_ind))]
        conv_last = bpp_save_last.loc[(bpp_save_last.bin_num.isin(conv_last_ind))]

        # select out the parameters just before disruption
        # first reset the index:
        conv_tot_reset = conv_tot.reset_index()
        conv_last_reset = conv_last.reset_index()

        # next select out the index for the disrupted systems using evol_type == 11
        conv_tot_reset_ind = conv_tot_reset.loc[conv_tot_reset.evol_type == 11.0].index
        conv_last_reset_ind = conv_last_reset.loc[conv_last_reset.evol_type == 11.0].index

        conv_tot = conv_tot_reset.iloc[conv_tot_reset_ind]
        conv_last = conv_last_reset.iloc[conv_last_reset_ind]

    elif method == 'final_state':
        # the bcm array is all that we need!
        conv_tot = bcm_save_tot
        conv_last = bcm_save_last

    elif method == 'XRB_form':
        # select out the systems which undergo a SN
        conv_tot_ind = bpp_save_tot.loc[bpp_save_tot.evol_type.isin([15.0, 16.0])].bin_num.unique()
        conv_last_ind = bpp_save_last.loc[bpp_save_last.evol_type.isin([15.0, 16.0])].bin_num.unique()
 
        conv_tot_sn = bpp_save_tot.loc[bpp_save_tot.bin_num.isin(conv_tot_ind)]
        conv_last_sn = bpp_save_last.loc[bpp_save_last.bin_num.isin(conv_last_ind)]

        # select out systems when they first enter RLO after the 1st SN
        conv_tot_xrb = conv_tot_sn.loc[(conv_tot_sn.kstar_1.isin(final_kstar_1)) &\
                                       (conv_tot_sn.kstar_2.isin(final_kstar_2)) &\
                                       (conv_tot_sn.RROL_2 >= 1.0) &\
                                       (conv_tot_sn.sep > 0)]
        conv_last_xrb = conv_last_sn.loc[(conv_last_sn.kstar_1.isin(final_kstar_1)) &\
                                         (conv_last_sn.kstar_2.isin(final_kstar_2)) &\
                                         (conv_last_sn.RROL_2 >= 1.0)& \
                                         (conv_last_sn.sep > 0)]
        conv_tot = conv_tot_xrb.groupby('bin_num').first().reset_index()
        conv_last = conv_last_xrb.groupby('bin_num').first().reset_index()

    return conv_tot, conv_last

def mass_min_max_select(kstar_1, kstar_2):
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

    primary_max = 150.0
    secondary_max = 150.0

    primary_min = 0.08
    secondary_min = 0.08

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
            min_mass[ii] = 10.0
        elif k == 13.0:
            min_mass[ii] = 6.0
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
        elif k== 11.0:
            max_mass[ii] = 20.0
        elif k <= 10.0:
            max_mass[ii] = 20.0
        ii += 1

    return min_mass[0], max_mass[0], min_mass[1], max_mass[1]


def idl_tabulate(x, f, p=5) :
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

    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in range(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


def rndm(a, b, g, size):
    """Power-law generator for pdf(x)\propto x^{g} for a<=x<=b

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

    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

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
    datZeroed = dat-datMin

    datTransformed = datZeroed/((datMax-datMin))
    if np.max(datTransformed) == 1.0:
        datTransformed[datTransformed == 1.0] = 1-1e-6
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
    for ii,column in zip(range(len(dat_list)),dat_list):
        dat_untrans = dat_exp[ii, :]*\
                   (max(dat_set[column]) - min(dat_set[column])) +\
                    min(dat_set[column])
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
        except:
            print('Using Scott Rule!!')
            bw = astrostats.scott_bin_width(dat)
        bw_list.append(bw)
    return np.mean(bw_list)


def error_check(BSEDict, filters=None, convergence=None):
    """Checks that values in BSEDict, filters, and convergence are viable
    """
    if not isinstance(BSEDict, dict):
        raise ValueError('BSE flags must be supplied via a dictionary')

    if filters is not None:
        if not isinstance(filters, dict):
            raise ValueError('Filters criteria must be supplied via a dictionary')

    if convergence is not None:
        if not isinstance(convergence, dict):
            raise ValueError('Convergence criteria must be supplied via a dictionary')

    # filters
    if filters is not None:
        flag='select_final_state'
        if flag in filters.keys():
            if filters[flag] not in [True,False]:
                raise ValueError("'{0:s}' needs to be either True or False (you set it to '{1:s}')".format(flag, filters[flag]))

        flag='binary_state'
        if flag in filters.keys():
            if any(x not in [0,1,2] for x in filters[flag]):
                raise ValueError("'{0:s}' needs to be a subset of [0,1,2] (you set it to '[{1:d}]')".format(flag, *filters[flag]))

    # convergence
    if convergence is not None:
        flag='conv_filter'
        if flag in convergence.keys():
            if not convergence[flag] in ['formation', '1_SN', '2_SN', 'disruption', 'final_state', 'XRB_form']:
                raise ValueError("'{0:s}' needs to be in the list: ['formation', '1_SN', '2_SN', 'disruption', 'final_state', 'XRB_form'] (you set it to '{1:s}'".format(flag, convergence[flag]))
        flag='match'
        if flag in convergence.keys():
            if not isinstance(convergence[flag], float):
                raise ValueError("'{0:s}' must be a float (you set it to '{1:0.2f}')".format(flag, convergence[flag]))


    # BSEDict
    flag='dtp'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater than or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='pts1'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='pts2'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='pts3'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))

    flag='windflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2,3]:
            raise ValueError("'{0:s}' needs to be set to either 0, 1, 2, or 3 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='eddlimflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise ValueError("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='neta'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='bwind'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='hewind'
    if flag in BSEDict.keys():
        if (BSEDict[flag] < 0) or (BSEDict[flag] > 1):
            raise ValueError("'{0:s}' needs to be between 0 and 1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='beta'
    # --- all numbers are valid
    flag='xi'
    if flag in BSEDict.keys():
        if (BSEDict[flag] < 0) or (BSEDict[flag] > 1):
            raise ValueError("'{0:s}' needs to be between 0 and 1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='acc2'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))

    flag='alpha1'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='lambdaf'
    if flag in BSEDict.keys():
        if (BSEDict[flag]>0) and (BSEDict[flag]!=1):
            raise ValueError("'{0:s}' needs to either be set to 1 for variable lambda or a negative number for fixed lambda (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='ceflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise ValueError("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='cekickflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2]:
            raise ValueError("'{0:s}' needs to be set to either 0, 1, or 2 (you set it to '{1:d}')".format(flag,BSEDict[flag]))
    flag='cemergeflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise ValueError("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='cehestarflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2]:
            raise ValueError("'{0:s}' needs to be set to either 0, 1, or 2 (you set it to '{1:d}')".format(flag,BSEDict[flag]))
    flag='qcflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2,3,4]:
            raise ValueError("'{0:s}' needs to be set to 0, 1, 2, or 3 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))

    flag='qcrit_array'
    if flag in BSEDict.keys():
        if any(x < 0.0 for x in BSEDict[flag]):
            raise ValueError("'{0:s}' values must be greater than or equal to zero (you set them to '[{1:d}]')".format(flag, *BSEDict[flag]))
        if len(BSEDict[flag]) != 16:
            raise ValueError("'{0:s}' must be supplied 16 values (you supplied '{1:d}')".format(flag, len(BSEDict[flag])))
        if (any( x != 0.0 for x in BSEDict[flag])) and (BSEDict['qcflag'] != 4):
            raise ValueError("If '{0:s}' is used, qcflag must be set to 4".format(flag))

    flag='sigma'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='bhflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2,3]:
            raise ValueError("'{0:s}' needs to be set to either 0, 1, 2, or 3 (you set it to '{1:d}')".format(flag,BSEDict[flag]))
    flag='ecsn'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='ecsn_mlow'
    if flag in BSEDict.keys():
        if (BSEDict[flag]>BSEDict['ecsn']) or (BSEDict[flag]<0.0):
            raise ValueError("'{0:s}' needs to be less than 'ecsn', and must be greater than or equal to 0 (you set it to '{0:0.2f}')".format(flag, BSEDict[flag]))
    flag='sigmadiv'
    if flag in BSEDict.keys():
        if BSEDict[flag] == 0:
            raise ValueError("'{0:s}' must be positive or negative (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='aic'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise valueerror("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='ussn'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise valueerror("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='pisn'
    if flag in BSEDict.keys():
        if not ((BSEDict[flag] > 0) or (BSEDict[flag] == -1) or (BSEDict[flag] == -2) or (BSEDict[flag] == -3)):
            raise ValueError("'{0:s}' needs to be set to either greater than 0 or equal to -1, -2, or -3 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='bhsigmafrac'
    if flag in BSEDict.keys():
        if (BSEDict[flag] <= 0) or (BSEDict[flag] > 1):
            raise ValueError("'{0:s}' needs to be greater than 0 and less than or equal to 1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='polar_kick_angle'
    if flag in BSEDict.keys():
        if (BSEDict[flag] < 0) or (BSEDict[flag] > 90):
            raise ValueError("'{0:s}' needs to be within the allowed range of [0,90] (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='natal_kick_array'
    if flag in BSEDict.keys():
        if len(BSEDict[flag]) != 6:
            raise ValueError("'{0:s}' must be supplied 6 values (you supplied '{1:d}')".format(flag, len(BSEDict[flag])))

    flag='nsflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2,3,4]:
            raise ValueError("'{0:s}' needs to be set to either 0, 1, 2, 3, or 4 (you set it to '{1:d}')".format(flag,BSEDict[flag]))
    flag='mxns'
    if flag in BSEDict.keys():
        if BSEDict[flag] <= 0:
            raise ValueError("'{0:s}' needs to be greater than 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))

    flag='eddfac'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='gamma'
    if flag in BSEDict.keys():
        if (BSEDict[flag]<=0) and (BSEDict[flag]!=-1) and (BSEDict[flag]!=-2):
            raise ValueError("'{0:s}' needs to either be set to -2, -1, or a positive number (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))

    flag='tflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1]:
            raise ValueError("'{0:s}' needs to be set to either 0 or 1 (you set it to '{1:d}')".format(flag, BSEDict[flag]))
    flag='ifflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='wdflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] < 0:
            raise ValueError("'{0:s}' needs to be greater or equal to 0 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='qcflag'
    if flag in BSEDict.keys():
        if BSEDict[flag] not in [0,1,2,3]:
            raise ValueError("'{0:s}' needs to be set to 0, 1, 2, or 3 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='epsnov'
    if flag in BSEDict.keys():
        if (BSEDict[flag] < 0) or (BSEDict[flag] > 1):
            raise ValueError("'{0:s}' needs to be between 0 and 1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='bconst'
    # --- all numbers are valid
    flag='ck'
    # --- all numbers are valid

    flag='fprimc_array'
    if flag in BSEDict.keys():
        if any(x < 0.0 for x in BSEDict[flag]):
            raise ValueError("'{0:s}' values must be greater than or equal to zero (you set them to '[{1:d}]')".format(flag, *BSEDict[flag]))
        if len(BSEDict[flag]) != 16:
            raise ValueError("'{0:s}' must be supplied 16 values (you supplied '{1:d}')".format(flag, len(BSEDict[flag])))

    return

def check_initial_conditions(initial_binary_table):
    """Checks initial conditions and reports warnings

        Only warning provided right now is if star begins in Roche lobe
        overflow
    """
    def rzamsf(m):
        """A function to evaluate Rzams
        ( from Tout et al., 1996, MNRAS, 281, 257 ).
        """
        mx = np.sqrt(m)
        rzams = (((a[7]*m**2 + a[8]*m**6)*mx + a[9]*m**11 +
                  (a[10] + a[11]*mx)*m**19)/
                  (a[12] + a[13]*m**2 + (a[14]*m**8 + m**18 + a[15]*m**19)*mx))

        return rzams

    def rl(Q):
        """A function to evaluate R_L/a(q), Eggleton 1983."""
        P = Q**(1.0/3.0)
        RL = 0.49*P*P/(0.6*P*P + np.log(1.0+P))
        return RL

    z = np.asarray(initial_binary_table['metallicity'])
    zpars, a = zcnsts(z)

    mass1 = np.asarray(initial_binary_table['mass1_binary'])
    mass2 = np.asarray(initial_binary_table['mass2_binary'])

    rzams1 = rzamsf(mass1)
    rzams2 = rzamsf(mass2)

    # assume some time step in order to calculate sep
    yeardy = 365.24
    aursun = 214.95
    tb = np.asarray(initial_binary_table['porb'])/yeardy
    sep = aursun*(tb*tb*(mass1 + mass2))**(1.0/3.0)

    q1 = mass1/mass2
    q2 = mass2/mass1
    rol1 = rl(q1)*sep
    rol2 = rl(q2)*sep

    # check for a ZAMS that starts in RFOL
    mask = ((np.array(initial_binary_table['kstar_1'])==1) & (rzams1 >= rol1)) | ((initial_binary_table['kstar_2']==1) & (rzams2 >= rol2))
    if mask.any():
        warnings.warn("At least one of your initial binaries is starting in Roche Lobe Overflow:\n{0}".format(initial_binary_table[mask]))

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
        0 : 'Main Sequence (MS), < 0.7 M⊙',
        1 : 'MS, > 0.7 M⊙',
        2 : 'Hertzsprung Gap',
        3 : 'First Giant Branch',
        4 : 'Core Helium Burning',
        5 : 'Early Asymptotic Giant Branch (AGB)',
        6 : 'Thermally Pulsing AGB',
        7 : 'Naked Helium Star MS',
        8 : 'Naked Helium Star Hertzsprung Gap',
        9 : 'Naked Helium Star Giant Branch',
        10 : 'Helium White Dwarf',
        11 : 'Carbon/Oxygen White Dwarf',
        12 : 'Oxygen/Neon White Dwarf',
        13 : 'Neutron Star',
        14 : 'Black Hole',
        15 : 'Massless Remnant',
    }

    kstar_string_to_int_dict = {v:k for k,v in kstar_int_to_string_dict.items()}

    evolve_type_int_to_string_dict = {
        1 : 'initial state',
        2 : 'kstar change',
        3 : 'begin Roche lobe overflow',
        4 : 'end Roche lobe overlow',
        5 : 'contact',
        6 : 'coalescence',
        7 : 'begin common envelope',
        8 : 'end common envelope',
        9 : 'no remnant leftover',
        10 : 'max evolution time',
        11 : 'binary disruption',
        12 : 'begin symbiotic phase',
        13 : 'end symbiotic phase',
        14 : 'blue straggler',
        15 : 'supernova of primary',
        16 : 'supernova of secondary',
    }

    evolve_type_string_to_int_dict = {v:k for k,v in evolve_type_int_to_string_dict.items()}

    if bpp.kstar_1.dtype in [int,float]:
        # convert from integer to string
        bpp['kstar_1'] = bpp['kstar_1'].astype(int)
        bpp['kstar_1'] = bpp['kstar_1'].apply(lambda x: kstar_int_to_string_dict[x])
    else:
        # convert from string to integer
        bpp['kstar_1'] = bpp['kstar_1'].apply(lambda x: kstar_string_to_int_dict[x])

    if bpp.kstar_2.dtype in [int,float]:
        # convert from integer to string
        bpp['kstar_2'] = bpp['kstar_2'].astype(int)
        bpp['kstar_2'] = bpp['kstar_2'].apply(lambda x: kstar_int_to_string_dict[x])
    else:
        # convert from string to integer
        bpp['kstar_2'] = bpp['kstar_2'].apply(lambda x: kstar_string_to_int_dict[x])

    if bpp.evol_type.dtype in [int,float]:
        # convert from integer to string
        bpp['evol_type'] = bpp['evol_type'].astype(int)
        bpp['evol_type'] = bpp['evol_type'].apply(lambda x: evolve_type_int_to_string_dict[x])
    else:
        # convert from string to integer
        bpp['evol_type'] = bpp['evol_type'].apply(lambda x: evolve_type_string_to_int_dict[x])

    return bpp

def parse_inifile(inifile):
    """Provides a method for parsing the inifile and returning dicts of each section
    """
    binOps = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Mod: operator.mod
    }

    def arithmetic_eval(s):
        """Allows us to control how the strings from the inifile get parses"""
        node = ast.parse(s, mode='eval')

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
                    'True': True,
                    'False': False,
                    'None': None,
                }
                return constants_lookup.get(
                    result.name,
                    result,
                )
            elif isinstance(node, ast.NameConstant):
                # None, True, False are nameconstants in python3, but names in 2
                return node.value
            else:
                raise Exception('Unsupported type {}'.format(node))

        return _eval(node.body)

    # ---- Create configuration-file-parser object and read parameters file.
    cp = ConfigParser()
    cp.optionxform = str
    cp.read(inifile)

    # ---- Read needed variables from the inifile
    dictionary = {}
    for section in cp.sections():
        dictionary[section] = {}
        for option in cp.options(section):
            opt = cp.get(section, option)
            try:
                dictionary[section][option] = arithmetic_eval(opt)
            except:
                dictionary[section][option] = json.loads(opt)
            if type(dictionary[section][option]) == VariableKey:
                dictionary[section][option] = dictionary[section][option].name

    BSEDict = dictionary['bse']
    seed_int = int(dictionary['rand_seed']['seed'])
    filters = dictionary['filters']
    convergence = dictionary['convergence']

    return BSEDict, seed_int, filters, convergence

class VariableKey(object):
    """
    A dictionary key which is a variable.
    @ivar item: The variable AST object.
    """
    def __init__(self, item):
        self.name = item.id

    def __eq__(self, compare):
        return (
            compare.__class__ == self.__class__
            and compare.name == self.name
        )

    def __hash__(self):
        return hash(self.name)
