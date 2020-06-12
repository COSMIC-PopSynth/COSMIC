# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017 - 2020)
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
from .checkstate import set_checkstates

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
               'Michael Zevin <zevin@northwestern.edu>',
               'digman.12@osu.edu']
__all__ = ['Evolve']


BPP_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2' ,
               'sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type',
               'aj_1', 'aj_2', 'tms_1', 'tms_2',
               'massc_1', 'massc_2', 'rad_1', 'rad_2',
               'mass0_1', 'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2',
               'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
               'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2',
               'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
               'bhspin_1','bhspin_2', 'bin_num']

BCM_COLUMNS = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lum_1', 'rad_1',
               'teff_1', 'massc_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'omega_spin_1', 'deltam_1', 'RRLO_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lum_2', 'rad_2', 'teff_2', 'massc_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'omega_spin_2', 'deltam_2', 'RRLO_2',
               'porb', 'sep', 'ecc', 'B_1', 'B_2',
               'SN_1', 'SN_2', 'bin_state', 'merger_type', 'bin_num']

KICK_COLUMNS = ['star', 'disrupted', 'natal_kick', 'phi', 'theta', 'eccentric_anomaly',
                'delta_vsysx_1', 'delta_vsysy_1', 'delta_vsysz_1', 'vsys_1_total',
                'delta_vsysx_2', 'delta_vsysy_2', 'delta_vsysz_2', 'vsys_2_total',
                'delta_theta_total', 'omega', 'randomseed', 'bin_num']

# We use the list of column in the initialbinarytable function to initialize
# the list of columns that we will send to the fortran evolv2 function.
# we also send this in a specific order so this he3lp ensures that the list that
# is created at the end has a consistent order
if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_BSE_COLUMNS = ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                             'ceflag', 'tflag', 'ifflag', 'wdflag', 'pisn', 'bhflag', 'remnantflag',
                             'cekickflag', 'cemergeflag', 'cehestarflag',
                             'mxns', 'pts1', 'pts2', 'pts3',
                             'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma', 'sigmadiv', 'bhsigmafrac', 'polar_kick_angle',
                             'natal_kick_array', 'qcrit_array',
                             'beta', 'xi', 'acc2', 'epsnov',
                             'eddfac', 'gamma', 'bdecayfac', 'bconst', 'ck', 
                             'windflag', 'qcflag', 'eddlimflag',
                             'fprimc_array', 'dtp', 'randomseed',
                             'bhspinflag','bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb',
                             'ST_cr', 'ST_tide', 'rembar_massloss', 'zsun','kickflag']

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

NATAL_KICK_COLUMNS = ['natal_kick',
                      'phi',
                      'theta',
                      'eccentric_anomaly',
                      'randomseed']

FLATTENED_NATAL_KICK_COLUMNS = []
for sn_idx in range(2):
    for idx, column_name in enumerate(NATAL_KICK_COLUMNS):
        FLATTENED_NATAL_KICK_COLUMNS.append(column_name + '_{0}'.format(sn_idx + 1))

QCRIT_COLUMNS = ['qcrit_{0}'.format(kstar) for kstar in range(0,16)]
FPRIMC_COLUMNS = ['fprimc_{0}'.format(kstar) for kstar in range(0,16)]

INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(FLATTENED_NATAL_KICK_COLUMNS)
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
                    for sn_idx in range(2):
                        column_name_sn = column_name + '_{0}'.format(sn_idx + 1)
                        kwargs1 = {column_name_sn : pd.Series([BSEDict['natal_kick_array'][sn_idx][idx]] * len(initialbinarytable), index=initialbinarytable.index, name=column_name_sn)}
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
        if (pd.Series(FLATTENED_NATAL_KICK_COLUMNS).isin(initialbinarytable.keys()).all()) and ('natal_kick_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(natal_kick_array=initialbinarytable[FLATTENED_NATAL_KICK_COLUMNS].values.reshape(-1,2,len(NATAL_KICK_COLUMNS)).tolist())

        if (pd.Series(QCRIT_COLUMNS).isin(initialbinarytable.keys()).all()) and ('qcrit_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(qcrit_array=initialbinarytable[QCRIT_COLUMNS].values.tolist())

        if (pd.Series(FPRIMC_COLUMNS).isin(initialbinarytable.keys()).all()) and ('fprimc_array' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(fprimc_array=initialbinarytable[FPRIMC_COLUMNS].values.tolist())

        # need to ensure that the order of parameters that we pass to BSE
        # is correct
        initial_conditions = initialbinarytable[INITIAL_CONDITIONS_PASS_COLUMNS].to_dict('records')

        # we use different columns to save the BSE parameters because some
        # of the parameters are list/arrays which we instead save as
        # individual values because it makes saving to HDF5 easier/more efficient.
        initialbinarytable = initialbinarytable[INITIAL_BINARY_TABLE_SAVE_COLUMNS]

        # Allow a user to specify a custom time step sampling for certain parts of the evolution
        timestep_conditions = kwargs.pop('timestep_conditions', [])
        set_checkstates(timestep_conditions=timestep_conditions)

        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                f['kick_info'] = np.zeros((2,len(KICK_COLUMNS)-1))
                # determine if we already have a compact object, if yes than one SN has already occured
                if (f['kstar_1'] in range(10,15)) or (f['kstar_2'] in range(10,15)):
                    f['kick_info'][0,0] = 1
                # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                _evolvebin.windvars.neta = f['neta']
                _evolvebin.windvars.bwind = f['bwind']
                _evolvebin.windvars.hewind = f['hewind']
                _evolvebin.cevars.alpha1 = f['alpha1']
                _evolvebin.cevars.lambdaf = f['lambdaf']
                _evolvebin.ceflags.ceflag = f['ceflag']
                _evolvebin.flags.tflag = f['tflag']
                _evolvebin.flags.ifflag = f['ifflag']
                _evolvebin.flags.wdflag = f['wdflag']
                _evolvebin.snvars.pisn = f['pisn']
                _evolvebin.flags.bhflag = f['bhflag']
                _evolvebin.flags.remnantflag = f['remnantflag']
                _evolvebin.ceflags.cekickflag = f['cekickflag']
                _evolvebin.ceflags.cemergeflag = f['cemergeflag']
                _evolvebin.ceflags.cehestarflag = f['cehestarflag']
                _evolvebin.snvars.mxns = f['mxns']
                _evolvebin.points.pts1 = f['pts1']
                _evolvebin.points.pts2 = f['pts2']
                _evolvebin.points.pts3 = f['pts3']
                _evolvebin.snvars.ecsn = f['ecsn']
                _evolvebin.snvars.ecsn_mlow = f['ecsn_mlow']
                _evolvebin.flags.aic = f['aic']
                _evolvebin.ceflags.ussn = f['ussn']
                _evolvebin.snvars.sigma = f['sigma']
                _evolvebin.snvars.sigmadiv = f['sigmadiv']
                _evolvebin.snvars.bhsigmafrac = f['bhsigmafrac']
                _evolvebin.snvars.polar_kick_angle = f['polar_kick_angle']
                _evolvebin.snvars.natal_kick_array = f['natal_kick_array']
                _evolvebin.cevars.qcrit_array = f['qcrit_array']
                _evolvebin.windvars.beta = f['beta']
                _evolvebin.windvars.xi = f['xi']
                _evolvebin.windvars.acc2 = f['acc2']
                _evolvebin.windvars.epsnov = f['epsnov']
                _evolvebin.windvars.eddfac = f['eddfac']
                _evolvebin.windvars.gamma = f['gamma']
                _evolvebin.flags.bdecayfac = f['bdecayfac']
                _evolvebin.magvars.bconst = f['bconst']
                _evolvebin.magvars.ck = f['ck']
                _evolvebin.flags.windflag = f['windflag']
                _evolvebin.flags.qcflag = f['qcflag']
                _evolvebin.flags.eddlimflag = f['eddlimflag']
                _evolvebin.tidalvars.fprimc_array = f['fprimc_array']
                _evolvebin.rand1.idum1 = f['randomseed']
                _evolvebin.flags.bhspinflag = f['bhspinflag']
                _evolvebin.snvars.bhspinmag = f['bhspinmag']
                _evolvebin.mixvars.rejuv_fac = f['rejuv_fac']
                _evolvebin.flags.rejuvflag = f['rejuvflag']
                _evolvebin.flags.htpmb = f['htpmb']
                _evolvebin.flags.st_cr = f['ST_cr']
                _evolvebin.flags.st_tide = f['ST_tide']
                _evolvebin.snvars.rembar_massloss = f['rembar_massloss']
                _evolvebin.metvars.zsun = f['zsun']
                _evolvebin.snvars.kickflag = f['kickflag']
                _evolvebin.cmcpass.using_cmc = 0

                [bpp, bcm, bpp_index, bcm_index, kick_info] = _evolvebin.evolv2([f['kstar_1'], f['kstar_2']],
                                                                                [f['mass_1'], f['mass_2']],
                                                                                f['porb'], f['ecc'], f['metallicity'], f['tphysf'], f['dtp'],
                                                                                [f['mass0_1'], f['mass0_2']],
                                                                                [f['rad_1'], f['rad_2']],
                                                                                [f['lum_1'], f['lum_2']],
                                                                                [f['massc_1'], f['massc_2']],
                                                                                [f['radc_1'], f['radc_2']],
                                                                                [f['menv_1'], f['menv_2']],
                                                                                [f['renv_1'], f['renv_2']],
                                                                                [f['omega_spin_1'], f['omega_spin_2']],
                                                                                [f['B_1'], f['B_2']],
                                                                                [f['bacc_1'], f['bacc_2']],
                                                                                [f['tacc_1'], f['tacc_2']],
                                                                                [f['epoch_1'], f['epoch_2']],
                                                                                [f['tms_1'], f['tms_2']],
                                                                                [f['bhspin_1'], f['bhspin_2']],
                                                                                f['tphys'],
                                                                                np.zeros(20),
                                                                                f['kick_info'])

                bcm = bcm[:bcm_index]
                bpp = bpp[:bpp_index]

                bpp = np.hstack((bpp, np.ones((bpp.shape[0], 1))*f['bin_num']))
                bcm = np.hstack((bcm, np.ones((bcm.shape[0], 1))*f['bin_num']))
                kick_info = np.hstack((kick_info, np.ones((kick_info.shape[0], 1))*f['bin_num']))

                return f, bpp, bcm, kick_info, _evolvebin.snvars.natal_kick_array

            except Exception as e:
                raise

        #define multiprocessing method to process the systems in batches 
        def _evolve_multi_system(f):
            try:
                res_bcm = np.zeros(f.shape[0],dtype=object)
                res_bpp = np.zeros(f.shape[0],dtype=object)
                res_kick_info = np.zeros(f.shape[0],dtype=object)
                res_natal_kick_array = np.zeros(f.shape[0],dtype=object)
                for i in range(0,f.shape[0]):
                    f[i]['kick_info'] = np.zeros((2,len(KICK_COLUMNS)-1))
                    # determine if we already have a compact object, if yes than one SN has already occured
                    if (f[i]['kstar_1'] in range(10,15)) or (f[i]['kstar_2'] in range(10,15)):
                        f[i]['kick_info'][0,0] = 1
                    _evolvebin.windvars.neta = f[i]['neta']
                    _evolvebin.windvars.bwind = f[i]['bwind']
                    _evolvebin.windvars.hewind = f[i]['hewind']
                    _evolvebin.cevars.alpha1 = f[i]['alpha1']
                    _evolvebin.cevars.lambdaf = f[i]['lambdaf']
                    _evolvebin.ceflags.ceflag = f[i]['ceflag']
                    _evolvebin.flags.tflag = f[i]['tflag']
                    _evolvebin.flags.ifflag = f[i]['ifflag']
                    _evolvebin.flags.wdflag = f[i]['wdflag']
                    _evolvebin.snvars.pisn = f[i]['pisn']
                    _evolvebin.flags.bhflag = f[i]['bhflag']
                    _evolvebin.flags.remnantflag = f[i]['remnantflag']
                    _evolvebin.ceflags.cekickflag = f[i]['cekickflag']
                    _evolvebin.ceflags.cemergeflag = f[i]['cemergeflag']
                    _evolvebin.ceflags.cehestarflag = f[i]['cehestarflag']
                    _evolvebin.snvars.mxns = f[i]['mxns']
                    _evolvebin.points.pts1 = f[i]['pts1']
                    _evolvebin.points.pts2 = f[i]['pts2']
                    _evolvebin.points.pts3 = f[i]['pts3']
                    _evolvebin.snvars.ecsn = f[i]['ecsn']
                    _evolvebin.snvars.ecsn_mlow = f[i]['ecsn_mlow']
                    _evolvebin.flags.aic = f[i]['aic']
                    _evolvebin.ceflags.ussn = f[i]['ussn']
                    _evolvebin.snvars.sigma = f[i]['sigma']
                    _evolvebin.snvars.sigmadiv = f[i]['sigmadiv']
                    _evolvebin.snvars.bhsigmafrac = f[i]['bhsigmafrac']
                    _evolvebin.snvars.polar_kick_angle = f[i]['polar_kick_angle']
                    _evolvebin.snvars.natal_kick_array = f[i]['natal_kick_array']
                    _evolvebin.cevars.qcrit_array = f[i]['qcrit_array']
                    _evolvebin.windvars.beta = f[i]['beta']
                    _evolvebin.windvars.xi = f[i]['xi']
                    _evolvebin.windvars.acc2 = f[i]['acc2']
                    _evolvebin.windvars.epsnov = f[i]['epsnov']
                    _evolvebin.windvars.eddfac = f[i]['eddfac']
                    _evolvebin.windvars.gamma = f[i]['gamma']
                    _evolvebin.flags.bdecayfac = f[i]['bdecayfac']
                    _evolvebin.magvars.bconst = f[i]['bconst']
                    _evolvebin.magvars.ck = f[i]['ck']
                    _evolvebin.flags.windflag = f[i]['windflag']
                    _evolvebin.flags.qcflag = f[i]['qcflag']
                    _evolvebin.flags.eddlimflag = f[i]['eddlimflag']
                    _evolvebin.tidalvars.fprimc_array = f[i]['fprimc_array']
                    _evolvebin.rand1.idum1 = f[i]['randomseed']
                    _evolvebin.flags.bhspinflag = f[i]['bhspinflag']
                    _evolvebin.snvars.bhspinmag = f[i]['bhspinmag']
                    _evolvebin.mixvars.rejuv_fac = f[i]['rejuv_fac']
                    _evolvebin.flags.rejuvflag = f[i]['rejuvflag']
                    _evolvebin.flags.htpmb = f[i]['htpmb']
                    _evolvebin.flags.st_cr = f[i]['ST_cr']
                    _evolvebin.flags.st_tide = f[i]['ST_tide']
                    _evolvebin.snvars.rembar_massloss = f[i]['rembar_massloss']
                    _evolvebin.metvars.zsun = f[i]['zsun']
                    _evolvebin.snvars.kickflag = f[i]['kickflag']
                    _evolvebin.cmcpass.using_cmc = 0

                    [bpp, bcm, bpp_index, bcm_index, kick_info_out] = _evolvebin.evolv2([f[i]['kstar_1'], f[i]['kstar_2']],
                                                                                    [f[i]['mass_1'], f[i]['mass_2']],
                                                                                    f[i]['porb'], f[i]['ecc'], f[i]['metallicity'], f[i]['tphysf'], f[i]['dtp'],
                                                                                    [f[i]['mass0_1'], f[i]['mass0_2']],
                                                                                    [f[i]['rad_1'], f[i]['rad_2']],
                                                                                    [f[i]['lum_1'], f[i]['lum_2']],
                                                                                    [f[i]['massc_1'], f[i]['massc_2']],
                                                                                    [f[i]['radc_1'], f[i]['radc_2']],
                                                                                    [f[i]['menv_1'], f[i]['menv_2']],
                                                                                    [f[i]['renv_1'], f[i]['renv_2']],
                                                                                    [f[i]['omega_spin_1'], f[i]['omega_spin_2']],
                                                                                    [f[i]['B_1'], f[i]['B_2']],
                                                                                    [f[i]['bacc_1'], f[i]['bacc_2']],
                                                                                    [f[i]['tacc_1'], f[i]['tacc_2']],
                                                                                    [f[i]['epoch_1'], f[i]['epoch_2']],
                                                                                    [f[i]['tms_1'], f[i]['tms_2']],
                                                                                    [f[i]['bhspin_1'], f[i]['bhspin_2']],
                                                                                    f[i]['tphys'],
                                                                                    np.zeros(20),
                                                                                    f[i]['kick_info'])

                    bpp = bpp[:bpp_index]
                    bcm = bcm[:bcm_index]

                    bpp_bin_numbers = np.atleast_2d(np.array([f[i]['bin_num']] * len(bpp))).T
                    bcm_bin_numbers = np.atleast_2d(np.array([f[i]['bin_num']] * len(bcm))).T
                    kick_info_out_bin_numbers = np.atleast_2d(np.array([f[i]['bin_num']] * len(kick_info_out))).T

                    res_bpp[i] = np.hstack((bpp, bpp_bin_numbers))
                    res_bcm[i] = np.hstack((bcm, bcm_bin_numbers))
                    res_kick_info[i] = np.hstack((kick_info_out, kick_info_out_bin_numbers))
                    res_natal_kick_array[i] = _evolvebin.snvars.natal_kick_array

                return f, np.vstack(res_bpp), np.vstack(res_bcm), np.vstack(res_kick_info), np.vstack(res_natal_kick_array)

            except Exception as e:
                raise

        # evolve systems
        if n_per_block > 0:
            initial_conditions = np.asarray(initial_conditions)
            n_tot = initial_conditions.shape[0]
            initial_conditions_blocked = []
            itr_block = 0
            while itr_block < n_tot:
                itr_next = np.min([n_tot,itr_block+n_per_block])
                initial_conditions_blocked.append(initial_conditions[itr_block:itr_next])
                itr_block = itr_next
            output = mp_utils.multiprocess_with_queues(
                nproc, _evolve_multi_system, initial_conditions_blocked, raise_exceptions=False)
        else:
            output = mp_utils.multiprocess_with_queues(
                nproc, _evolve_single_system, initial_conditions, raise_exceptions=False)

        output = np.array(output)
        bpp_arrays = np.vstack(output[:, 1])
        bcm_arrays = np.vstack(output[:, 2])
        kick_info_arrays = np.vstack(output[:, 3])

        natal_kick_arrays = np.vstack(output[:, 4])
        natal_kick_arrays = natal_kick_arrays.reshape(-1,1,len(FLATTENED_NATAL_KICK_COLUMNS))
        for idx, column in enumerate(FLATTENED_NATAL_KICK_COLUMNS):
                # assigning values this way work for most of the parameters.
                kwargs1 = {column : natal_kick_arrays[:,:,idx]}
                initialbinarytable = initialbinarytable.assign(**kwargs1)

        kick_info = pd.DataFrame(kick_info_arrays,
                           columns=KICK_COLUMNS,
                           index=kick_info_arrays[:, -1].astype(int))

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

        return bpp, bcm, initialbinarytable, kick_info
