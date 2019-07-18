# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
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


def calc_Roche_radius(M1, M2, A):
    """ Get Roche lobe radius (Eggleton 1983)

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
        in units of input, A
    """
    q = M1 / M2
    return A * 0.49*q**(2.0/3.0) / (0.6*q**(2.0/3.0) + np.log(1.0 + q**(1.0/3.0)))


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
        bpp or bcm table, e.g. ``{'disrupted_binaries' : False}``;
        This means you do *not* want disrupted_binaries in your table

    kstar1_range : `list`
        list containing all kstar1 values to retain

    kstar2_range : `list`
        list containing all kstar2 values to retain

    Returns
    -------
    bcm : `pandas.DataFrame`
        filtered bcm dataframe
    """
    _known_methods = ['mass_transfer_white_dwarf_to_co',
                      'select_final_state',
                      'binary_state',
                      'lisa_sources']

    if not set(method.keys()).issubset(set(_known_methods)):
        raise ValueError("You have supplied an "
                         "unknown method to filter out "
                         "the bpp or bcm array. Known methods are "
                         "{0}".format(_known_methods))

    for meth, use in method.items():
        if meth == 'mass_transfer_white_dwarf_to_co' and not use:
            idx_save_1 = bpp.loc[~(bpp.kstar_1.isin([10,11,12,13,14])) &
                                  (bpp.kstar_2.isin([10,11,12])) &
                                  (bpp.RROL_1 > 1)].bin_num.tolist()
            idx_save_2 = bpp.loc[~(bpp.kstar_1.isin([10,11,12,13,14])) &
                                  (bpp.kstar_2.isin([10,11,12])) &
                                  (bpp.RROL_1 > 2)].bin_num.tolist()
            idx_save = np.intersect1d(idx_save_1, idx_save_2)
            bcm = bcm.loc[~bcm.bin_num.isin(idx_save)]
        elif (meth == 'select_final_state') and use:
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

        elif (meth == 'lisa_sources') and use:
            if 0 in method['binary_state']:
                bcm_0 = bcm.loc[bcm.bin_state==0]
                bcm_0_LISAflag = bcm_0.loc[bcm_0.porb > 5].bin_num
                bcm = bcm.loc[~bcm.bin_num.isin(bcm_0_LISAflag)]
            else:
                raise ValueError("You must have bin state = 0 for lisa"
                                 "sources filter")
    return bcm, bin_state_fraction

def bcm_conv_select(bcm_save_tot, bcm_save_last, method):
    """Select bcm data for special convergence cases

    Parameters
    ----------
    bcm_save_tot : `pandas.DataFrame`
        bcm dataframe containing all saved bcm data

    bcm_save_last : `pandas.DataFrame`
        bcm dataframe containing bcm data from last
        iteration

    method : `dict`,
        one or more methods by which to filter the
        bcm table, e.g. ``{'lisa_convergence' : True}``;
        This means you want to only compute the convergence
        over the region specified for the lisa_convergence
        method below

    Returns
    -------
    bcm_conv_tot : `pandas.DataFrame`
        filtered bcm dataframe containing all saved bcm
        data

    bcm_conv_last : `pandas.DataFrame`
        filtered bcm dataframe containing saved bcm
        data from last iteration

    """
    _known_methods = ['lisa_convergence']

    if not set(method.keys()).issubset(set(_known_methods)):
        raise ValueError("You have supplied an "
                         "unknown method to filter the "
                         "bcm array for convergence. Known methods are "
                         "{0}".format(_known_methods))
    bcm_conv_tot = bcm_save_tot
    if len(bcm_save_tot) == len(bcm_save_last):
        bcm_conv_last = bcm_save_last
    else:
        bcm_conv_last = bcm_save_tot.iloc[:len(bcm_save_tot)-len(bcm_save_last)]
    for meth, use in method.items():
        if meth == 'lisa_convergence' and use:
            bcm_conv_tot = bcm_conv_tot.loc[bcm_conv_tot.porb < np.log10(5000)]
            bcm_conv_last = bcm_conv_last.loc[bcm_conv_last.porb < np.log10(5000)]

    # If it is the first iteration, divide the bcm array in two
    # for convergence computation
    if len(bcm_conv_tot) == len(bcm_conv_last):
        bcm_conv_last = bcm_conv_tot.iloc[:int(len(bcm_conv_tot)/2)]

    return bcm_conv_tot, bcm_conv_last

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
            min_mass[ii] = 5.0
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
        flag='mass_transfer_white_dwarf_to_co'
        if flag in filters.keys():
            if filters[flag] not in [True,False]:
                raise ValueError("'{0:s}' needs to be either True or False (you set it to '{1:s}')".format(flag, filters[flag]))

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
        flag='lisa_convergence'
        if flag in convergence.keys():
            if convergence[flag] not in [True,False]:
                raise ValueError("'{0:s}' needs to be either True or False (you set it to '{1:s}')".format(flag, convergence[flag]))

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
    flag='qcrit_array'
    if flag in BSEDict.keys():
        if any(x < 0.0 for x in BSEDict[flag]):
            raise ValueError("'{0:s}' values must be greater than or equal to zero (you set them to '[{1:d}]')".format(flag, *BSEDict[flag]))
        if len(BSEDict[flag]) != 16:
            raise ValueError("'{0:s}' must be supplied 16 values (you supplied '{1:d}')".format(flag, len(BSEDict[flag])))

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
        if not ((BSEDict[flag] > 0) or (BSEDict[flag] == -1)):
            raise ValueError("'{0:s}' needs to be set to either greater than 0 or -1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
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
    flag='epsnov'
    if flag in BSEDict.keys():
        if (BSEDict[flag] < 0) or (BSEDict[flag] > 1):
            raise ValueError("'{0:s}' needs to be between 0 and 1 (you set it to '{1:0.2f}')".format(flag, BSEDict[flag]))
    flag='bconst'
    # --- all numbers are valid
    flag='ck'
    # --- all numbers are valid

    return
