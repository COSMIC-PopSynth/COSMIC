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

"""`evolve`
"""

from cosmic import _evolvebin
from . import utils
from .sample import initialbinarytable
from .checkstate import set_checkstates

from schwimmbad import MultiPool

import numpy as np
import pandas as pd
import warnings
import os
import sys
try:
    import multiprocessing
    multiprocessing.set_start_method("fork")
except RuntimeError:
    pass


__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__credits__ = ['Katelyn Breivik <katie.breivik@gmail.com>',
               'Michael Zevin <zevin@northwestern.edu>',
               'digman.12@osu.edu',
               'Tom Wagg <tomjwagg@gmail.com>']
__all__ = ['Evolve']


# Make this match the ordering of all_cols in bpp_array.f
ALL_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'porb',
               'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 'aj_1', 'aj_2', 'tms_1',
               'tms_2', 'massc_1', 'massc_2', 'rad_1', 'rad_2', 'mass0_1',
               'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2', 'radc_1',
               'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2', 'omega_spin_1',
               'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2', 'tacc_1',
               'tacc_2', 'epoch_1', 'epoch_2', 'bhspin_1', 'bhspin_2',
               'deltam_1', 'deltam_2', 'SN_1', 'SN_2', 'bin_state', 'merger_type']

INTEGER_COLUMNS = ["bin_state", "bin_num", "kstar_1", "kstar_2", "SN_1", "SN_2", "evol_type"]


BPP_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2',
               'sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type',
               'aj_1', 'aj_2', 'tms_1', 'tms_2',
               'massc_1', 'massc_2', 'rad_1', 'rad_2',
               'mass0_1', 'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2',
               'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
               'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2',
               'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
               'bhspin_1', 'bhspin_2']

BCM_COLUMNS = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lum_1', 'rad_1',
               'teff_1', 'massc_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'omega_spin_1', 'deltam_1', 'RRLO_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lum_2', 'rad_2', 'teff_2', 'massc_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'omega_spin_2', 'deltam_2', 'RRLO_2',
               'porb', 'sep', 'ecc', 'B_1', 'B_2',
               'SN_1', 'SN_2', 'bin_state', 'merger_type']

KICK_COLUMNS = ['star', 'disrupted', 'natal_kick', 'phi', 'theta', 'mean_anomaly',
                'delta_vsysx_1', 'delta_vsysy_1', 'delta_vsysz_1', 'vsys_1_total',
                'delta_vsysx_2', 'delta_vsysy_2', 'delta_vsysz_2', 'vsys_2_total',
                'theta_euler', 'phi_euler', 'psi_euler', 'randomseed', 'bin_num']

# We use the list of column in the initialbinarytable function to initialize
# the list of columns that we will send to the fortran evolv2 function.
# we also send this in a specific order so this help ensures that the list that
# is created at the end has a consistent order
if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_BSE_COLUMNS = ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                                  'ceflag', 'tflag', 'ifflag', 'wdflag', 'pisn', 'rtmsflag',
                                  'bhflag', 'remnantflag', 'grflag', 'bhms_coll_flag', 'wd_mass_lim',
                                  'cekickflag', 'cemergeflag', 'cehestarflag',
                                  'mxns', 'pts1', 'pts2', 'pts3',
                                  'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma', 'sigmadiv',
                                  'bhsigmafrac', 'polar_kick_angle',
                                  'natal_kick_array', 'qcrit_array',
                                  'beta', 'xi', 'acc2', 'epsnov',
                                  'eddfac', 'gamma', 'don_lim', 'acc_lim', 
                                  'bdecayfac', 'bconst', 'ck',
                                  'windflag', 'qcflag', 'eddlimflag',
                                  'fprimc_array', 'dtp', 'randomseed',
                                  'bhspinflag', 'bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb',
                                  'ST_cr', 'ST_tide', 'rembar_massloss', 'zsun', 'kickflag']

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
                      'mean_anomaly',
                      'randomseed']

FLATTENED_NATAL_KICK_COLUMNS = []
for sn_idx in range(2):
    for idx, column_name in enumerate(NATAL_KICK_COLUMNS):
        FLATTENED_NATAL_KICK_COLUMNS.append(column_name + '_{0}'.format(sn_idx + 1))

QCRIT_COLUMNS = ['qcrit_{0}'.format(kstar) for kstar in range(0, 16)]
FPRIMC_COLUMNS = ['fprimc_{0}'.format(kstar) for kstar in range(0, 16)]

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
    def evolve(self, initialbinarytable, pool=None, bpp_columns=None, bcm_columns=None, **kwargs):
        """After setting a number of initial conditions we evolve the system.

        Parameters
        ----------
        initialbinarytable : DataFrame
            Initial conditions of the binary

        pool : Multiprocessing pool
            Pool of workers to use to evolve systems in parallel

        bpp_columns : list, optional, default: None
            Columns to save in the bpp table (key evolutionary stage table)

        bcm_columns : list, optional, default: None
            Columns to save in the bcm table (detailed evolution table)

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
        output_bpp : :class:`pandas.DataFrame`
            Table of key evolutionary stages for each binary

        output_bcm : :class:`pandas.DataFrame`
            Table of detailed evolution for each binary

        initialbinarytable : DataFrame
            Initial conditions for each binary
        """
        idx = kwargs.pop('idx', 0)
        nproc = min(kwargs.pop('nproc', 1), len(initialbinarytable))
        n_per_block = kwargs.pop('n_per_block', -1)

        if bpp_columns is None:
            bpp_columns = BPP_COLUMNS
        if bcm_columns is None:
            bcm_columns = BCM_COLUMNS

        # There are three ways to tell evolve and thus the fortran
        # what you want all the flags and other BSE specific
        # parameters to be

        # NUMBER 1: PASS A DICTIONARY OF FLAGS
        BSEDict = kwargs.pop('BSEDict', {})

        # NUMBER 2: PASS A PANDAS DATA FRAME WITH PARAMS DEFINED AS COLUMNS

        #     All you need is the initialbinarytable with columns,
        #     If you pass both a dictionary of flags and/or a inifile
        #     and a table with the columns, the column values will be
        #     overwritten by what is in the dictionary or ini file

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
            seed = np.random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max, size=len(initialbinarytable))
            initialbinarytable = initialbinarytable.assign(randomseed=kwargs.pop('randomseed', seed))
        if 'bin_num' not in initialbinarytable.keys():
            initialbinarytable = initialbinarytable.assign(bin_num=np.arange(idx, idx + len(initialbinarytable)))

        for k, v in BSEDict.items():
            if k in initialbinarytable.keys():
                warnings.warn("The value for {0} in initial binary table is being "
                              "overwritten by the value of {0} from either the params "
                              "file or the BSEDict.".format(k))
            # special columns that need to be handled differently
            if k == 'natal_kick_array':
                assign_natal_kick_array = [BSEDict['natal_kick_array']] * len(initialbinarytable)
                initialbinarytable = initialbinarytable.assign(natal_kick_array=assign_natal_kick_array)
                for idx, column_name in enumerate(NATAL_KICK_COLUMNS):
                    for sn_idx in range(2):
                        column_name_sn = column_name + '_{0}'.format(sn_idx + 1)
                        column_values = pd.Series([BSEDict['natal_kick_array'][sn_idx][idx]] * len(initialbinarytable),
                                                  index=initialbinarytable.index,
                                                  name=column_name_sn)
                        kwargs1 = {column_name_sn: column_values}
                        initialbinarytable = initialbinarytable.assign(**kwargs1)
            elif k == 'qcrit_array':
                initialbinarytable = initialbinarytable.assign(qcrit_array=[BSEDict['qcrit_array']] * len(initialbinarytable))
                for kstar in range(0, 16):
                    columns_values = pd.Series([BSEDict['qcrit_array'][kstar]] * len(initialbinarytable),
                                               index=initialbinarytable.index,
                                               name='qcrit_{0}'.format(kstar))
                    initialbinarytable.loc[:, 'qcrit_{0}'.format(kstar)] = columns_values
            elif k == 'fprimc_array':
                columns_values = [BSEDict['fprimc_array']] * len(initialbinarytable)
                initialbinarytable = initialbinarytable.assign(fprimc_array=columns_values)
                for kstar in range(0, 16):
                    columns_values = pd.Series([BSEDict['fprimc_array'][kstar]] * len(initialbinarytable),
                                               index=initialbinarytable.index,
                                               name='fprimc_{0}'.format(kstar))
                    initialbinarytable.loc[:, 'fprimc_{0}'.format(kstar)] = columns_values
            else:
                # assigning values this way work for most of the parameters.
                kwargs1 = {k: v}
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
        if ((pd.Series(FLATTENED_NATAL_KICK_COLUMNS).isin(initialbinarytable.keys()).all()) and
           ('natal_kick_array' not in BSEDict)):
            column_values = initialbinarytable[FLATTENED_NATAL_KICK_COLUMNS].values.reshape(-1,
                                                                                            2,
                                                                                            len(NATAL_KICK_COLUMNS)).tolist()
            initialbinarytable = initialbinarytable.assign(natal_kick_array=column_values)

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

        # set the indices of the columns to include in bpp table (+1 because fortran is 1-indexed)
        col_inds_bpp = np.zeros(len(ALL_COLUMNS), dtype=int)
        col_inds_bpp[:len(bpp_columns)] = [ALL_COLUMNS.index(col) + 1 for col in bpp_columns]

        # save bpp column information in the initial conditions
        for i in range(len(initial_conditions)):
            initial_conditions[i]["n_col_bpp"] = len(bpp_columns)
            initial_conditions[i]["col_inds_bpp"] = col_inds_bpp

        # same for bcm
        col_inds_bcm = np.zeros(len(ALL_COLUMNS), dtype=int)
        col_inds_bcm[:len(bcm_columns)] = [ALL_COLUMNS.index(col) + 1 for col in bcm_columns]
        for i in range(len(initial_conditions)):
            initial_conditions[i]["n_col_bcm"] = len(bcm_columns)
            initial_conditions[i]["col_inds_bcm"] = col_inds_bcm

        # check if a pool was passed
        if pool is None:
            with MultiPool(processes=nproc) as pool:
                # evolve systems
                if n_per_block > 0:
                    initial_conditions = np.asarray(initial_conditions)
                    n_tot = initial_conditions.shape[0]
                    initial_conditions_blocked = []
                    itr_block = 0
                    while itr_block < n_tot:
                        itr_next = np.min([n_tot, itr_block+n_per_block])
                        initial_conditions_blocked.append(initial_conditions[itr_block:itr_next])
                        itr_block = itr_next
                    output = list(pool.map(_evolve_multi_system, initial_conditions_blocked))
                else:
                    output = list(pool.map(_evolve_single_system, initial_conditions))
        else:
            # evolve systems
            if n_per_block > 0:
                initial_conditions = np.asarray(initial_conditions)
                n_tot = initial_conditions.shape[0]
                initial_conditions_blocked = []
                itr_block = 0
                while itr_block < n_tot:
                    itr_next = np.min([n_tot, itr_block+n_per_block])
                    initial_conditions_blocked.append(initial_conditions[itr_block:itr_next])
                    itr_block = itr_next
                output = list(pool.map(_evolve_multi_system, initial_conditions_blocked))
            else:
                output = list(pool.map(_evolve_single_system, initial_conditions))

        output = np.array(output, dtype=object)
        bpp_arrays = np.vstack(output[:, 1])
        bcm_arrays = np.vstack(output[:, 2])
        kick_info_arrays = np.vstack(output[:, 3])

        natal_kick_arrays = np.vstack(output[:, 4])
        natal_kick_arrays = natal_kick_arrays.reshape(-1, 1, len(FLATTENED_NATAL_KICK_COLUMNS))
        for idx, column in enumerate(FLATTENED_NATAL_KICK_COLUMNS):
            # assigning values this way work for most of the parameters.
            kwargs1 = {column: natal_kick_arrays[:, :, idx]}
            initialbinarytable = initialbinarytable.assign(**kwargs1)

        kick_info = pd.DataFrame(kick_info_arrays,
                                 columns=KICK_COLUMNS,
                                 index=kick_info_arrays[:, -1].astype(int))

        bpp = pd.DataFrame(bpp_arrays,
                           columns=bpp_columns + ["bin_num"],
                           index=bpp_arrays[:, -1].astype(int))

        bcm = pd.DataFrame(bcm_arrays,
                           columns=bcm_columns + ["bin_num"],
                           index=bcm_arrays[:, -1].astype(int))

        # convert a subset of columns to integers
        for col in INTEGER_COLUMNS:
            if col in bpp.columns:
                bpp[col] = bpp[col].astype(int)
            if col in bcm.columns:
                bcm[col] = bcm[col].astype(int)

        # convert merger type to a padded string
        if 'merger_type' in bpp.columns:
            bpp.merger_type = bpp.merger_type.astype(int).astype(str).apply(lambda x: x.zfill(4))
        if 'merger_type' in bcm.columns:
            bcm.merger_type = bcm.merger_type.astype(int).astype(str).apply(lambda x: x.zfill(4))

        return bpp, bcm, initialbinarytable, kick_info


def _evolve_single_system(f):
    try:
        f["kick_info"] = np.zeros((2, len(KICK_COLUMNS)-1))
        # determine if we already have a compact object, if yes than one SN has already occured
        if (f["kstar_1"] in range(10, 15)) or (f["kstar_2"] in range(10, 15)):
            f["kick_info"][0, 0] = 1
        # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
        _evolvebin.windvars.neta = f["neta"]
        _evolvebin.windvars.bwind = f["bwind"]
        _evolvebin.windvars.hewind = f["hewind"]
        _evolvebin.cevars.alpha1 = f["alpha1"]
        _evolvebin.cevars.lambdaf = f["lambdaf"]
        _evolvebin.ceflags.ceflag = f["ceflag"]
        _evolvebin.flags.tflag = f["tflag"]
        _evolvebin.flags.ifflag = f["ifflag"]
        _evolvebin.flags.wdflag = f["wdflag"]
        _evolvebin.flags.rtmsflag = f["rtmsflag"]
        _evolvebin.snvars.pisn = f["pisn"]
        _evolvebin.flags.bhflag = f["bhflag"]
        _evolvebin.flags.remnantflag = f["remnantflag"]
        _evolvebin.ceflags.cekickflag = f["cekickflag"]
        _evolvebin.ceflags.cemergeflag = f["cemergeflag"]
        _evolvebin.ceflags.cehestarflag = f["cehestarflag"]
        _evolvebin.flags.grflag = f["grflag"]
        _evolvebin.flags.bhms_coll_flag = f["bhms_coll_flag"]
        _evolvebin.flags.wd_mass_lim = f["wd_mass_lim"]
        _evolvebin.snvars.mxns = f["mxns"]
        _evolvebin.points.pts1 = f["pts1"]
        _evolvebin.points.pts2 = f["pts2"]
        _evolvebin.points.pts3 = f["pts3"]
        _evolvebin.snvars.ecsn = f["ecsn"]
        _evolvebin.snvars.ecsn_mlow = f["ecsn_mlow"]
        _evolvebin.flags.aic = f["aic"]
        _evolvebin.ceflags.ussn = f["ussn"]
        _evolvebin.snvars.sigma = f["sigma"]
        _evolvebin.snvars.sigmadiv = f["sigmadiv"]
        _evolvebin.snvars.bhsigmafrac = f["bhsigmafrac"]
        _evolvebin.snvars.polar_kick_angle = f["polar_kick_angle"]
        _evolvebin.snvars.natal_kick_array = f["natal_kick_array"]
        _evolvebin.cevars.qcrit_array = f["qcrit_array"]
        _evolvebin.mtvars.don_lim = f["don_lim"]
        _evolvebin.mtvars.acc_lim = f["acc_lim"]
        _evolvebin.windvars.beta = f["beta"]
        _evolvebin.windvars.xi = f["xi"]
        _evolvebin.windvars.acc2 = f["acc2"]
        _evolvebin.windvars.epsnov = f["epsnov"]
        _evolvebin.windvars.eddfac = f["eddfac"]
        _evolvebin.windvars.gamma = f["gamma"]
        _evolvebin.flags.bdecayfac = f["bdecayfac"]
        _evolvebin.magvars.bconst = f["bconst"]
        _evolvebin.magvars.ck = f["ck"]
        _evolvebin.flags.windflag = f["windflag"]
        _evolvebin.flags.qcflag = f["qcflag"]
        _evolvebin.flags.eddlimflag = f["eddlimflag"]
        _evolvebin.tidalvars.fprimc_array = f["fprimc_array"]
        _evolvebin.rand1.idum1 = f["randomseed"]
        _evolvebin.flags.bhspinflag = f["bhspinflag"]
        _evolvebin.snvars.bhspinmag = f["bhspinmag"]
        _evolvebin.mixvars.rejuv_fac = f["rejuv_fac"]
        _evolvebin.flags.rejuvflag = f["rejuvflag"]
        _evolvebin.flags.htpmb = f["htpmb"]
        _evolvebin.flags.st_cr = f["ST_cr"]
        _evolvebin.flags.st_tide = f["ST_tide"]
        _evolvebin.snvars.rembar_massloss = f["rembar_massloss"]
        _evolvebin.metvars.zsun = f["zsun"]
        _evolvebin.snvars.kickflag = f["kickflag"]
        _evolvebin.cmcpass.using_cmc = 0

        _evolvebin.col.n_col_bpp = f["n_col_bpp"]
        _evolvebin.col.col_inds_bpp = f["col_inds_bpp"]
        _evolvebin.col.n_col_bcm = f["n_col_bcm"]
        _evolvebin.col.col_inds_bcm = f["col_inds_bcm"]

        [bpp_index, bcm_index, kick_info] = _evolvebin.evolv2([f["kstar_1"], f["kstar_2"]],
                                                              [f["mass_1"], f["mass_2"]],
                                                              f["porb"], f["ecc"], f["metallicity"], 
                                                              f["tphysf"], f["dtp"],
                                                              [f["mass0_1"], f["mass0_2"]],
                                                              [f["rad_1"], f["rad_2"]],
                                                              [f["lum_1"], f["lum_2"]],
                                                              [f["massc_1"], f["massc_2"]],
                                                              [f["radc_1"], f["radc_2"]],
                                                              [f["menv_1"], f["menv_2"]],
                                                              [f["renv_1"], f["renv_2"]],
                                                              [f["omega_spin_1"], f["omega_spin_2"]],
                                                              [f["B_1"], f["B_2"]],
                                                              [f["bacc_1"], f["bacc_2"]],
                                                              [f["tacc_1"], f["tacc_2"]],
                                                              [f["epoch_1"], f["epoch_2"]],
                                                              [f["tms_1"], f["tms_2"]],
                                                              [f["bhspin_1"], f["bhspin_2"]],
                                                              f["tphys"],
                                                              np.zeros(20),
                                                              np.zeros(20),
                                                              f["kick_info"])
        bpp = _evolvebin.binary.bpp[:bpp_index, :f["n_col_bpp"]].copy()
        _evolvebin.binary.bpp[:bpp_index, :f["n_col_bpp"]] = np.zeros(bpp.shape)
        bcm = _evolvebin.binary.bcm[:bcm_index, :f["n_col_bcm"]].copy()
        _evolvebin.binary.bcm[:bcm_index, :f["n_col_bcm"]] = np.zeros(bcm.shape)

        bpp = np.hstack((bpp, np.ones((bpp.shape[0], 1))*f["bin_num"]))
        bcm = np.hstack((bcm, np.ones((bcm.shape[0], 1))*f["bin_num"]))
        kick_info = np.hstack((kick_info, np.ones((kick_info.shape[0], 1))*f["bin_num"]))

        return f, bpp, bcm, kick_info, _evolvebin.snvars.natal_kick_array.copy()

    except Exception as e:
        print(e)
        raise


def _evolve_multi_system(f):
    try:
        res_bcm = np.zeros(f.shape[0], dtype=object)
        res_bpp = np.zeros(f.shape[0], dtype=object)
        res_kick_info = np.zeros(f.shape[0], dtype=object)
        res_natal_kick_array = np.zeros(f.shape[0], dtype=object)
        for i in range(0, f.shape[0]):

            # call evolve single system
            _, bpp, bcm, kick_info, _ = _evolve_single_system(f[i])

            # add results to pre-allocated list
            res_bpp[i] = bpp
            res_bcm[i] = bcm
            res_kick_info[i] = kick_info
            res_natal_kick_array[i] = _evolvebin.snvars.natal_kick_array

        return f, np.vstack(res_bpp), np.vstack(res_bcm), np.vstack(res_kick_info), np.vstack(res_natal_kick_array)

    except Exception as e:
        print(e)
        raise
