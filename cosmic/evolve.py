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

"""`evolve`
"""

from cosmic import _evolvebin
from . import utils
from .sample import initialbinarytable

from configparser import ConfigParser
from .mp import mp as mp_utils

import numpy as np
import pandas as pd
import json
import warnings
import os
import sys

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__credits__ = ['Katelyn Breivik <katie.breivik@gmail.com>',
               'Michael Zevin <zevin@northwestern.edu>']
__all__ = ['Evolve']


BPP_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2' ,
               'sep', 'porb', 'ecc', 'RROL_1', 'RROL_2', 'evol_type',
               'Vsys_1', 'Vsys_2', 'SNkick', 'SNtheta',
               'aj_1', 'aj_2', 'tms_1', 'tms_2',
               'massc_1', 'massc_2', 'rad_1', 'rad_2',
               'bin_num']

BCM_COLUMNS = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lumin_1', 'rad_1',
               'teff_1', 'massc_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'ospin_1', 'deltam_1', 'RROL_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lumin_2', 'rad_2', 'teff_2', 'massc_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'ospin_2', 'deltam_2', 'RROL_2',
               'porb', 'sep', 'ecc', 'B_0_1', 'B_0_2',
               'SNkick_1', 'SNkick_2', 'Vsys_final', 'SNtheta_final',
               'SN_1', 'SN_2', 'bin_state', 'merger_type', 'bin_num']

# We use the list of column in the initialbinarytable function to initialize
# the list of columns that we will send to the fortran evolv2 function.
# we also send this in a specific order so this he3lp ensures that the list that
# is created at the end has a consistent order
if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_BSE_COLUMNS = ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                             'ceflag', 'tflag', 'ifflag', 'wdflag', 'pisn', 'bhflag', 'nsflag',
                             'cekickflag', 'cemergeflag', 'cehestarflag',
                             'mxns', 'pts1', 'pts2', 'pts3',
                             'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma', 'sigmadiv', 'bhsigmafrac', 'polar_kick_angle',
                             'natal_kick_array', 'qcrit_array',
                             'beta', 'xi', 'acc2', 'epsnov',
                             'eddfac', 'gamma', 'bdecayfac', 'bconst', 'ck', 
                             'windflag', 'qcflag', 'eddlimflag',
                             'fprimc_array', 'dtp', 'randomseed',
                             'bhspinflag','bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb',
                             'ST_cr', 'ST_tide']

INITIAL_CONDITIONS_MISC_COLUMN = ['bin_num']

# Add the BSE COLUMSN and MISC COLUMN to the PASS_COLUMNS list
INITIAL_CONDITIONS_PASS_COLUMNS.extend(INITIAL_CONDITIONS_BSE_COLUMNS)
INITIAL_CONDITIONS_PASS_COLUMNS.extend(INITIAL_CONDITIONS_MISC_COLUMN)

if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_BINARY_TABLE_SAVE_COLUMNS = INITIAL_CONDITIONS_PASS_COLUMNS[:]
else:
    INITIAL_BINARY_TABLE_SAVE_COLUMNS = INITIAL_CONDITIONS_PASS_COLUMNS.copy()

for col in ['natal_kick_array', 'qcrit_array', 'fprimc_array']:
    INITIAL_BINARY_TABLE_SAVE_COLUMNS.remove(col)

NATAL_KICK_COLUMNS = ['SNkick_1', 'SNkick_2', 'phi_1', 'phi_2', 'theta_1', 'theta_2']
QCRIT_COLUMNS = ['qcrit_{0}'.format(kstar) for kstar in range(0,16)]
FPRIMC_COLUMNS = ['fprimc_{0}'.format(kstar) for kstar in range(0,16)]

INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(NATAL_KICK_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(QCRIT_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(FPRIMC_COLUMNS)

# BSE doesn't need the binary fraction, so just add to columns for saving
INITIAL_BINARY_TABLE_SAVE_COLUMNS.insert(7, 'binfrac')

class Evolve(object):
    def __init__():
        '''
        initialize Evolve
        '''

    @classmethod
    def evolve(cls, initialbinarytable, **kwargs):
        """After setting a number of initial conditions we evolve the system.

        Parameters
        ----------
        initialbinarytable : DataFrame
            Initial conditions of the binary

        **kwargs:
            There are three ways to tell evolve and thus the fortran
            what you want all the flags and other BSE specific
            parameters to be. If you pass both a dictionary of flags and/or a inifile
            and a table with the BSE parameters in the columns,
            the column values will be overwritten by
            what is in the dictionary or ini file.

            NUMBER 1: PASS A DICTIONARY OF FLAGS

                 BSEDict

            NUMBER 2: PASS A PANDAS DATA FRAME WITH PARAMS DEFINED AS COLUMNS

                 All you need is the initialbinarytable if the all
                 the BSE parameters are defined as columns

            NUMBER 3: PASS PATH TO A INI FILE WITH THE FLAGS DEFINED

                params

        randomseed : `int`, optional, default let numpy choose for you
            If you would like the random seed that the underlying fortran code
            uses to be the same for all of the initial conditions you passed
            then you can send this keyword argument in. It is recommended
            to just let numpy choose a random number as the Fortran random seed
            and then this number will be returned as a column in the
            initial binary table so that you can reproduce the results.

        nproc : `int`, optional, default: 1
            number of CPUs to use to evolve systems
            in parallel

        idx : `int`, optional, default: 0
            initial index of the bcm/bpp arrays

        dtp : `float`, optional: default: tphysf
            timestep size in Myr for bcm output where tphysf
            is total evolution time in Myr

        n_per_block : `int`, optional, default: -1
            number of systems to evolve in a block with 
            _evolve_multi_system, to allow larger multiprocessing
            queues and reduced overhead. If less than 1 use _evolve_single_system

        Returns
        -------
        output_bpp : DataFrame
            Evolutionary history of each binary

        output_bcm : DataFrame
            Final state of each binary

        initialbinarytable : DataFrame
            Initial conditions for each binary
        """
        idx = kwargs.pop('idx', 0)
        nproc = min(kwargs.pop('nproc', 1), len(initialbinarytable))
        n_per_block = kwargs.pop('n_per_block',-1)

        # There are three ways to tell evolve and thus the fortran
        # what you want all the flags and other BSE specific
        # parameters to be

        # NUMBER 1: PASS A DICTIONARY OF FLAGS
        BSEDict = kwargs.pop('BSEDict', {})

        # NUMBER 2: PASS A PANDAS DATA FRAME WITH PARAMS DEFINED AS COLUMNS

            # All you need is the initialbinarytable with columns,
            # If you pass both a dictionary of flags and/or a inifile
            # and a table with the columns, the column values will be
            # overwritten by what is in the dictionary or ini file

        # NUMBER 3: PASS PATH TO A INI FILE WITH THE FLAGS DEFINED
        params = kwargs.pop('params', None)

        if BSEDict and params is not None:
            raise ValueError('Please pass either a dictionary '
                             'of BSE flags or a path to an inifle not both.')

        if params is not None:
            if not os.path.isfile(params):
                raise ValueError("File does not exist, probably supplied incorrect "
                                 "path to the inifile.")
            BSEDict, _, _, _, _ = utils.parse_inifile(params)

        # error check the parameters you are trying to pass to BSE
        # if we sent in a table with the parameter names
        # then we will temporarily create a dictionary
        # in order to verify that the values in the table
        # are valid
        utils.error_check(BSEDict)

        # check the initial conditions of the system and warn user if
        # anything is weird about them, such as the star starts
        # in Roche Lobe overflow
        utils.check_initial_conditions(initialbinarytable)


        # assign some columns based on keyword arguments but that
        # can be overwritten by the params or BSEDict
        if 'dtp' not in initialbinarytable.keys():
            initialbinarytable = initialbinarytable.assign(dtp=kwargs.pop('dtp', initialbinarytable['tphysf']))
        if 'randomseed' not in initialbinarytable.keys():
            initialbinarytable = initialbinarytable.assign(randomseed=kwargs.pop('randomseed',
                                                                                 np.random.randint(np.iinfo(np.int32).min,
                                                                                 np.iinfo(np.int32).max,
                                                                                 size=len(initialbinarytable))
                                                                                 )
                                                           )
        if 'bin_num' not in initialbinarytable.keys():
            initialbinarytable = initialbinarytable.assign(bin_num=np.arange(idx, idx + len(initialbinarytable)))

        for k,v in BSEDict.items():
            if k in initialbinarytable.keys():
                warnings.warn("The value for {0} in initial binary table is being "
                              "overwritten by the value of {0} from either the params "
                              "file or the BSEDict.".format(k))
            # special columns that need to be handled differently
            if k == 'natal_kick_array':
                initialbinarytable = initialbinarytable.assign(natal_kick_array=[BSEDict['natal_kick_array']] * len(initialbinarytable))
                for idx, column_name in enumerate(NATAL_KICK_COLUMNS):
                    kwargs1 = {column_name : pd.Series([BSEDict['natal_kick_array'][idx]] * len(initialbinarytable), index=initialbinarytable.index, name=column_name)}
                    initialbinarytable = initialbinarytable.assign(**kwargs1)
            elif k == 'qcrit_array':
                initialbinarytable = initialbinarytable.assign(qcrit_array=[BSEDict['qcrit_array']] * len(initialbinarytable))
                for kstar in range(0,16):
                    initialbinarytable.loc[:, 'qcrit_{0}'.format(kstar)] = pd.Series([BSEDict['qcrit_array'][kstar]]* len(initialbinarytable), index=initialbinarytable.index, name='qcrit_{0}'.format(kstar))
            elif k == 'fprimc_array':
                initialbinarytable = initialbinarytable.assign(fprimc_array=[BSEDict['fprimc_array']] * len(initialbinarytable))
                for kstar in range(0,16):
                    initialbinarytable.loc[:, 'fprimc_{0}'.format(kstar)] = pd.Series([BSEDict['fprimc_array'][kstar]]* len(initialbinarytable), index=initialbinarytable.index, name='fprimc_{0}'.format(kstar))
            else:
                # assigning values this way work for most of the parameters.
                kwargs1 = {k:v}
                initialbinarytable = initialbinarytable.assign(**kwargs1)

        # Here we perform two checks
        # First, if the BSE parameters are not in the initial binary table
        # and either a dictionary or an inifile was not provided
        # then we need to raise an ValueError and tell the user to provide
        # either a dictionary or an inifile or add more columns
        if not BSEDict:
            if ((not set(INITIAL_BINARY_TABLE_SAVE_COLUMNS).issubset(initialbinarytable.columns)) and
                (not set(INITIAL_CONDITIONS_PASS_COLUMNS).issubset(initialbinarytable.columns))):
                raise ValueError("You are passing BSE parameters as columns in the "
                                 "initial binary table but not all BSE parameters are defined. "
                                 "Please pass a BSEDict or a params file or make sure "
                                 "you have all BSE parameters as columns {0} or {1}.".format(
                                  INITIAL_BINARY_TABLE_SAVE_COLUMNS, INITIAL_CONDITIONS_PASS_COLUMNS))

        # If you did not supply the natal kick or qcrit_array or fprimc_array in the BSEdict then we construct
        # it from the initial conditions table
        if (pd.Series(NATAL_KICK_COLUMNS).isin(initialbinarytable.keys()).all()) and ('natal_kick_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(natal_kick_array=initialbinarytable[NATAL_KICK_COLUMNS].values.tolist())

        if (pd.Series(QCRIT_COLUMNS).isin(initialbinarytable.keys()).all()) and ('qcrit_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(qcrit_array=initialbinarytable[QCRIT_COLUMNS].values.tolist())

        if (pd.Series(FPRIMC_COLUMNS).isin(initialbinarytable.keys()).all()) and ('fprimc_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(fprimc_array=initialbinarytable[FPRIMC_COLUMNS].values.tolist())

        # need to ensure that the order of parameters that we pass to BSE
        # is correct
        initial_conditions = initialbinarytable[INITIAL_CONDITIONS_PASS_COLUMNS].values

        # we use different columns to save the BSE parameters because some
        # of the parameters are list/arrays which we instead save as
        # individual values because it makes saving to HDF5 easier/more efficient.
        initialbinarytable = initialbinarytable[INITIAL_BINARY_TABLE_SAVE_COLUMNS]

        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                _evolvebin.windvars.neta = f[37]
                _evolvebin.windvars.bwind = f[38]
                _evolvebin.windvars.hewind = f[39]
                _evolvebin.cevars.alpha1 = f[40]
                _evolvebin.cevars.lambdaf = f[41]
                _evolvebin.ceflags.ceflag = f[42]
                _evolvebin.flags.tflag = f[43]
                _evolvebin.flags.ifflag = f[44]
                _evolvebin.flags.wdflag = f[45]
                _evolvebin.snvars.pisn = f[46]
                _evolvebin.flags.bhflag = f[47]
                _evolvebin.flags.nsflag = f[48]
                _evolvebin.ceflags.cekickflag = f[49]
                _evolvebin.ceflags.cemergeflag = f[50]
                _evolvebin.ceflags.cehestarflag = f[51]
                _evolvebin.windvars.mxns = f[52]
                _evolvebin.points.pts1 = f[53]
                _evolvebin.points.pts2 = f[54]
                _evolvebin.points.pts3 = f[55]
                _evolvebin.snvars.ecsn = f[56]
                _evolvebin.snvars.ecsn_mlow = f[57]
                _evolvebin.flags.aic = f[58]
                _evolvebin.ceflags.ussn = f[59]
                _evolvebin.snvars.sigma = f[60]
                _evolvebin.snvars.sigmadiv = f[61]
                _evolvebin.snvars.bhsigmafrac = f[62]
                _evolvebin.snvars.polar_kick_angle = f[63]
                _evolvebin.snvars.natal_kick_array = f[64]
                _evolvebin.cevars.qcrit_array = f[65]
                _evolvebin.windvars.beta = f[66]
                _evolvebin.windvars.xi = f[67]
                _evolvebin.windvars.acc2 = f[68]
                _evolvebin.windvars.epsnov = f[69]
                _evolvebin.windvars.eddfac = f[70]
                _evolvebin.windvars.gamma = f[71]
                _evolvebin.flags.bdecayfac = f[72]
                _evolvebin.magvars.bconst = f[73]
                _evolvebin.magvars.ck = f[74]
                _evolvebin.flags.windflag = f[75]
                _evolvebin.flags.qcflag = f[76]
                _evolvebin.windvars.eddlimflag = f[77]
                _evolvebin.tidalvars.fprimc_array = f[78]
                _evolvebin.rand1.idum1 = f[80]
                _evolvebin.flags.bhspinflag = f[81]
                _evolvebin.snvars.bhspinmag = f[82]
                _evolvebin.mixvars.rejuv_fac = f[83]
                _evolvebin.flags.rejuvflag = f[84]
                _evolvebin.flags.htpmb = f[85]
                _evolvebin.flags.st_cr = f[86]
                _evolvebin.flags.st_tide = f[87]
                _evolvebin.cmcpass.using_cmc = 0

                [bpp, bcm] = _evolvebin.evolv2([f[0],f[1]], [f[2],f[3]], f[4], f[5], f[6], f[7], f[79],
                                                [f[8],f[9]], [f[10],f[11]], [f[12],f[13]],
                                                [f[14],f[15]], [f[16],f[17]], [f[18],f[19]],
                                                [f[20],f[21]], [f[22],f[23]], [f[24],f[25]],
                                                [f[26],f[27]], [f[28],f[29]], [f[30],f[31]],
                                                [f[32],f[33]], [f[34],f[35]], f[36],
                                                np.zeros(20),np.zeros(20))

                try:
                    bpp = bpp[:np.argwhere(bpp[:,0] == -1)[0][0]]
                    bcm = bcm[:np.argwhere(bcm[:,0] == -1)[0][0]]
                except IndexError:
                    bpp = bpp[:np.argwhere(bpp[:,0] > 0)[0][0]]
                    bcm = bcm[:np.argwhere(bcm[:,0] > 0)[0][0]]
                    raise Warning('bpp overload: mass1 = {0}, mass2 = {1}, porb = {2}, ecc = {3}, tphysf = {4}, metallicity = {5}'\
                                   .format(f[2], f[3], f[4], f[5], f[7], f[6]))

                bpp_bin_numbers = np.atleast_2d(np.array([f[88]] * len(bpp))).T
                bcm_bin_numbers = np.atleast_2d(np.array([f[88]] * len(bcm))).T

                bpp = np.hstack((bpp, bpp_bin_numbers))
                bcm = np.hstack((bcm, bcm_bin_numbers))

                return f, bpp, bcm

            except Exception as e:
                raise

        #define multiprocessing method to process the systems in batches 
        def _evolve_multi_system(f):
            try:
                res_bcm = np.zeros(f.shape[0],dtype=object)
                res_bpp = np.zeros(f.shape[0],dtype=object)
                for i in range(0,f.shape[0]):
                    # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                    _evolvebin.windvars.neta = f[i,37]
                    _evolvebin.windvars.bwind = f[i,38]
                    _evolvebin.windvars.hewind = f[i,39]
                    _evolvebin.cevars.alpha1 = f[i,40]
                    _evolvebin.cevars.lambdaf = f[i,41]
                    _evolvebin.ceflags.ceflag = f[i,42]
                    _evolvebin.flags.tflag = f[i,43]
                    _evolvebin.flags.ifflag = f[i,44]
                    _evolvebin.flags.wdflag = f[i,45]
                    _evolvebin.snvars.pisn = f[i,46]
                    _evolvebin.flags.bhflag = f[i,47]
                    _evolvebin.flags.nsflag = f[i,48]
                    _evolvebin.ceflags.cekickflag = f[i,49]
                    _evolvebin.ceflags.cemergeflag = f[i,50]
                    _evolvebin.ceflags.cehestarflag = f[i,51]
                    _evolvebin.windvars.mxns = f[i,52]
                    _evolvebin.points.pts1 = f[i,53]
                    _evolvebin.points.pts2 = f[i,54]
                    _evolvebin.points.pts3 = f[i,55]
                    _evolvebin.snvars.ecsn = f[i,56]
                    _evolvebin.snvars.ecsn_mlow = f[i,57]
                    _evolvebin.flags.aic = f[i,58]
                    _evolvebin.ceflags.ussn = f[i,59]
                    _evolvebin.snvars.sigma = f[i,60]
                    _evolvebin.snvars.sigmadiv = f[i,61]
                    _evolvebin.snvars.bhsigmafrac = f[i,62]
                    _evolvebin.snvars.polar_kick_angle = f[i,63]
                    _evolvebin.snvars.natal_kick_array = f[i,64]
                    _evolvebin.cevars.qcrit_array = f[i,65]
                    _evolvebin.windvars.beta = f[i,66]
                    _evolvebin.windvars.xi = f[i,67]
                    _evolvebin.windvars.acc2 = f[i,68]
                    _evolvebin.windvars.epsnov = f[i,69]
                    _evolvebin.windvars.eddfac = f[i,70]
                    _evolvebin.windvars.gamma = f[i,71]
                    _evolvebin.flags.bdecayfac = f[i,72]
                    _evolvebin.magvars.bconst = f[i,73]
                    _evolvebin.magvars.ck = f[i,74]
                    _evolvebin.flags.windflag = f[i,75]
                    _evolvebin.flags.qcflag = f[i,76]
                    _evolvebin.windvars.eddlimflag = f[i,77]
                    _evolvebin.tidalvars.fprimc_array = f[i,78]
                    _evolvebin.rand1.idum1 = f[i,80]
                    _evolvebin.flags.bhspinflag = f[i,81]
                    _evolvebin.snvars.bhspinmag = f[i,82]
                    _evolvebin.mixvars.rejuv_fac = f[i,83]
                    _evolvebin.flags.rejuvflag = f[i,84]
                    _evolvebin.flags.htpmb = f[i,85]
                    _evolvebin.flags.st_cr = f[i,86]
                    _evolvebin.flags.st_tide = f[i,87]
                    _evolvebin.cmcpass.using_cmc = 0 
                    [bpp, bcm] = _evolvebin.evolv2([f[i,0],f[i,1]], [f[i,2],f[i,3]], f[i,4], f[i,5], f[i,6], f[i,7], f[i,79],
                                                    [f[i,8],f[i,9]], [f[i,10],f[i,11]], [f[i,12],f[i,13]],
                                                    [f[i,14],f[i,15]], [f[i,16],f[i,17]], [f[i,18],f[i,19]],
                                                    [f[i,20],f[i,21]], [f[i,22],f[i,23]], [f[i,24],f[i,25]],
                                                    [f[i,26],f[i,27]], [f[i,28],f[i,29]], [f[i,30],f[i,31]],
                                                    [f[i,32],f[i,33]], [f[i,34],f[i,35]], f[i,36],
                                                    np.zeros(20),np.zeros(20))
                    try:
                        idx1 = np.argmax(bpp[:,0] == -1)
                        idx2 = np.argmax(bcm[:,0] == -1)
                        bpp = bpp[:idx1]
                        bcm = bcm[:idx2]
                    except IndexError:
                        bpp = bpp[:np.argwhere(bpp[:,0] > 0)[0][0]]
                        bcm = bcm[:np.argwhere(bcm[:,0] > 0)[0][0]]
                        raise Warning('bpp overload: mass1 = {0}, mass2 = {1}, porb = {2}, ecc = {3}, tphysf = {4}, metallicity = {5}'\
                                       .format(f[i,2], f[i,3], f[i,4], f[i,5], f[i,7], f[i,6]))
                    bpp_bin_numbers = np.atleast_2d(np.array([f[i,88]] * len(bpp))).T
                    bcm_bin_numbers = np.atleast_2d(np.array([f[i,88]] * len(bcm))).T

                    res_bpp[i] = np.hstack((bpp, bpp_bin_numbers))
                    res_bcm[i] = np.hstack((bcm, bcm_bin_numbers))

                return f, np.vstack(res_bpp), np.vstack(res_bcm)

            except Exception as e:
                raise

        # evolve systems
        if n_per_block > 0:
            n_tot = initial_conditions.shape[0]
            initial_conditions_blocked = []
            itr_block = 0
            while itr_block < n_tot:
                itr_next = np.min([n_tot,itr_block+n_per_block])
                initial_conditions_blocked.append(initial_conditions[itr_block:itr_next,:])
                itr_block = itr_next
            output = mp_utils.multiprocess_with_queues(
                nproc, _evolve_multi_system, initial_conditions_blocked, raise_exceptions=False)
        else:
            output = mp_utils.multiprocess_with_queues(
                nproc, _evolve_single_system, initial_conditions, raise_exceptions=False)

        output = np.array(output)
        bpp_arrays = np.vstack(output[:, 1])
        bcm_arrays = np.vstack(output[:, 2])

        bpp = pd.DataFrame(bpp_arrays,
                           columns=BPP_COLUMNS,
                           index=bpp_arrays[:, -1].astype(int))

        bcm = pd.DataFrame(bcm_arrays,
                           columns=BCM_COLUMNS,
                           index=bcm_arrays[:, -1].astype(int))

        bcm.merger_type = bcm.merger_type.astype(int).astype(str).apply(lambda x: x.zfill(4))
        bcm.bin_state = bcm.bin_state.astype(int)
        bpp.bin_num = bpp.bin_num.astype(int)
        bcm.bin_num = bcm.bin_num.astype(int)

        return bpp, bcm, initialbinarytable
