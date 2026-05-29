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
import tqdm
from functools import partial
from pathlib import Path
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
               'tms_2', 'massc_he_layer_1', 'massc_he_layer_2', 'massc_co_layer_1', 'massc_co_layer_2',
               'rad_1', 'rad_2', 'mass0_1',
               'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2', 'radc_1',
               'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2', 'omega_spin_1',
               'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2', 'tacc_1',
               'tacc_2', 'epoch_1', 'epoch_2', 'bhspin_1', 'bhspin_2',
               'deltam_1', 'deltam_2', 'SN_1', 'SN_2', 'bin_state', 'merger_type', 'metallicity']

INTEGER_COLUMNS = ["bin_state", "bin_num", "kstar_1", "kstar_2", "SN_1", "SN_2", "evol_type"]


BPP_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2',
               'sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type',
               'aj_1', 'aj_2', 'tms_1', 'tms_2',
               'massc_he_layer_1', 'massc_he_layer_2', 'massc_co_layer_1', 'massc_co_layer_2', 'rad_1', 'rad_2',
               'mass0_1', 'mass0_2', 'lum_1', 'lum_2', 'teff_1', 'teff_2',
               'radc_1', 'radc_2', 'menv_1', 'menv_2', 'renv_1', 'renv_2',
               'omega_spin_1', 'omega_spin_2', 'B_1', 'B_2', 'bacc_1', 'bacc_2',
               'tacc_1', 'tacc_2', 'epoch_1', 'epoch_2',
               'bhspin_1', 'bhspin_2']

BCM_COLUMNS = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lum_1', 'rad_1',
               'teff_1', 'massc_he_layer_1', 'massc_co_layer_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'omega_spin_1', 'deltam_1', 'RRLO_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lum_2', 'rad_2', 'teff_2', 'massc_he_layer_2', 'massc_co_layer_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'omega_spin_2', 'deltam_2', 'RRLO_2',
               'porb', 'sep', 'ecc', 'B_1', 'B_2',
               'SN_1', 'SN_2', 'bin_state', 'merger_type']

KICK_COLUMNS = ['star', 'disrupted', 'natal_kick', 'phi', 'theta', 'mean_anomaly',
                'delta_vsysx_1', 'delta_vsysy_1', 'delta_vsysz_1', 'vsys_1_total',
                'delta_vsysx_2', 'delta_vsysy_2', 'delta_vsysz_2', 'vsys_2_total',
                'theta_euler', 'phi_euler', 'psi_euler', 'randomseed', 'tphys', 'bin_num']

# We use the list of column in the initialbinarytable function to initialize
# the list of columns that we will send to the fortran evolv2 function.
# we also send this in a specific order so this help ensures that the list that
# is created at the end has a consistent order
if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_PASS_COLUMNS = initialbinarytable.INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_BSE_COLUMNS = ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                                  'ce2stageflag', 'ceflag', 'tflag', 'ifflag', 'wdflag',
                                  'pisn', 'ppi_co_shift', 'ppi_extra_ml',
                                  'rtmsflag',
                                  'bhflag', 'remnantflag', 'fryer_mass_limit',
                                  'maltsev_mode', 'maltsev_fallback', 'maltsev_pf_prob',
                                  'grflag', 'bhms_coll_flag', 'wd_mass_lim',
                                  'cekickflag', 'cemergeflag', 'cehestarflag',
                                  'mxns', 'pts1', 'pts2', 'pts3',
                                  "fryer_fmix", "fryer_mcrit_nsbh",
                                  'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma', 'sigmadiv',
                                  'bhsigmafrac', 'polar_kick_angle', 'mm_mu_ns', 'mm_mu_bh',
                                  'natal_kick_array', 'qcrit_array',
                                  'beta', 'xi', 'acc2', 'epsnov',
                                  'eddfac', 'gamma', 'don_lim', 'acc_lim', 'smt_periastron_check',
                                  'bdecayfac', 'bconst', 'ck',
                                  'windflag', 'qcflag', 'eddlimflag', 'LBV_flag',
                                  'fprimc_array', 'dtp', 'randomseed',
                                  'bhspinflag', 'bhspinmag', 'rejuv_fac', 'rejuvflag', 'htpmb',
                                  'ST_cr', 'ST_tide', 'rembar_massloss', 'zsun', 'kickflag']

INITIAL_CONDITIONS_MISC_COLUMN = ['bin_num']

INITIAL_CONDITIONS_SSE_COLUMN = ['stellar_engine','path_to_tracks','path_to_he_tracks','z_accuracy_limit']

# Add the BSE COLUMSN and MISC COLUMN to the PASS_COLUMNS list
INITIAL_CONDITIONS_PASS_COLUMNS.extend(INITIAL_CONDITIONS_BSE_COLUMNS)
INITIAL_CONDITIONS_PASS_COLUMNS.extend(INITIAL_CONDITIONS_MISC_COLUMN)
INITIAL_CONDITIONS_PASS_COLUMNS.extend(INITIAL_CONDITIONS_SSE_COLUMN)

if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_BINARY_TABLE_SAVE_COLUMNS = INITIAL_CONDITIONS_PASS_COLUMNS[:]
else:
    INITIAL_BINARY_TABLE_SAVE_COLUMNS = INITIAL_CONDITIONS_PASS_COLUMNS.copy()

for col in ['natal_kick_array', 'qcrit_array', 'fprimc_array', 'alpha1', 'acc_lim']:
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
ALPHA_COLUMNS = ['alpha1_{0}'.format(star) for star in range(0, 2)]
ACCLIM_COLUMNS = ['acc_lim_{0}'.format(star) for star in range(0, 2)]

INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(FLATTENED_NATAL_KICK_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(QCRIT_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(FPRIMC_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(ALPHA_COLUMNS)
INITIAL_BINARY_TABLE_SAVE_COLUMNS.extend(ACCLIM_COLUMNS)

# BSE doesn't need the binary fraction, so just add to columns for saving
INITIAL_BINARY_TABLE_SAVE_COLUMNS.insert(7, 'binfrac')


class Evolve(object):
    def __init__():
        '''
        initialize Evolve
        '''

    @classmethod
    def evolve(self, initialbinarytable, pool=None, bpp_columns=None, bcm_columns=None,
               dt_mass_modifiers=[(40, 70, 0.3), (70, np.inf, 0.1)], **kwargs):
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

        dt_modifiers : list of tuples, optional, default: [(40, 70, 0.3), (70, np.inf, 0.1)]
            List of tuples specifying the mass ranges and corresponding modifiers for the timestep size.
            Our recommended default improves the numerical stability at higher masses.
            Each tuple should be of the form (m_low, m_high, mod) and will modify the default timestep
            by a factor of mod for systems with a *primary* mass in the range m_low <= mass_1 < m_high.
            For example, (40, 70, 0.3) would multiply the default timestep size by 0.3 for systems with
            primary mass between [40, 70) solar masses. We apply the modifier to the pts1, pts2, and pts3
            parameters which control the timestep size in different evolutionary phases. These changes
            are logged in the initial conditions table so you can keep track of which systems had their
            timesteps modified. Avoid overlapping mass ranges for different modifiers as this will result
            in multiple modifiers being applied in the overlap region.
            NOTE: these modifiers are only applied to columns which aren't present in the initialbinarytable
            that is passed in (i.e. they only modify values provided by a BSEDict or params.ini file)

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

            You can also add a progress bar by setting: progress=True

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
        progress = kwargs.pop('progress', False)

        if bpp_columns is None:
            bpp_columns = BPP_COLUMNS
        if bcm_columns is None:
            bcm_columns = BCM_COLUMNS

        columns_in_passed_initC = set(initialbinarytable.columns)

        # There are three ways to tell evolve and thus the fortran
        # what you want all the flags and other BSE specific
        # parameters to be

        # NUMBER 1: PASS A DICTIONARY OF FLAGS
        BSEDict = kwargs.pop('BSEDict', {})
        SSEDict = kwargs.pop('SSEDict', {})


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
            BSEDict, SSEDict, _, _, _, _ = utils.parse_inifile(params)

        # default to SSE when no SSEDict is provided
        if BSEDict and not SSEDict:
            SSEDict = {'stellar_engine': 'sse'}

        # error check the parameters you are trying to pass to BSE
        # if we sent in a table with the parameter names
        # then we will temporarily create a dictionary
        # in order to verify that the values in the table
        # are valid
        utils.error_check(BSEDict, SSEDict)
        
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

        # ensure SSEDict keys are unique in the initial binary table and warn otherwise
        for col in INITIAL_CONDITIONS_SSE_COLUMN:
            if col in initialbinarytable.columns and initialbinarytable[col].nunique() > 1:
                raise ValueError(f"The value for {col} in the initial binary table is not unique. "
                                 f"Please make sure that the value for {col} is the same for all rows in the initial binary table.")

        # if user passed an SSEDict, we may need to update the initial binary table
        if SSEDict:
            for k, v in SSEDict.items():
                if k in initialbinarytable.keys():
                    warnings.warn(f"The value for {k} in initial binary table is being "
                                  f"overwritten by the value of {k} from either the params "
                                  f"file or the SSEDict.")
                kwargs1 = {k: v}
                initialbinarytable = initialbinarytable.assign(**kwargs1)

        # if the user wants to use METISSE then we need to load the necessary tracks
        if initialbinarytable['stellar_engine'].iloc[0] == 'metisse':
            _evolvebin.se_flags.using_metisse = 1
            _evolvebin.se_flags.using_sse = 0

            # make sure all of the SSE columns are there
            if not set(['path_to_tracks', 'path_to_he_tracks', 'z_accuracy_limit']).issubset(initialbinarytable.columns):
                raise ValueError("If you want to use the METISSE stellar engine, you must provide the following in the SSEDict, initial binary table, or params file: path_to_tracks, path_to_he_tracks, z_accuracy_limit.")

            #check if the metallicity for the initialbinarytable changes
            # raise an error if all the metallicities are not the same
            if initialbinarytable['metallicity'].nunique() > 1:
                raise ValueError("All the metallicities in the initial binary table "
                                 "must be the same if you are using the METISSE stellar engine. ")
            
            # load in the METISSE files
            m_min, m_max = read_tracks_for_METISSE(
                path_to_tracks=initialbinarytable['path_to_tracks'].iloc[0], 
                IBT_Z=initialbinarytable['metallicity'].iloc[0],
                z_accuracy_limit=initialbinarytable['z_accuracy_limit'].iloc[0],
                is_he=False
            )

            if (initialbinarytable['path_to_he_tracks'].iloc[0] != ''):
                read_tracks_for_METISSE(
                    path_to_tracks=initialbinarytable['path_to_he_tracks'].iloc[0],
                    IBT_Z=initialbinarytable['metallicity'].iloc[0], 
                    z_accuracy_limit=initialbinarytable['z_accuracy_limit'].iloc[0],
                    is_he=True
                )

        else:
            # default to SSE if stellar engine is SSE or no stellar engine is specified
            _evolvebin.se_flags.using_sse = 1
            _evolvebin.se_flags.using_metisse = 0
           
        # go through each item in the BSEDict and update the initialbinarytable
        new_cols = {}
        n = len(initialbinarytable)
        idx = initialbinarytable.index
        for k, v in list(BSEDict.items()):
            # warn the user if they are overwriting a value
            if k in initialbinarytable.columns:
                warnings.warn(
                    "The value for {0} in initial binary table is being overwritten by the value of {0} "
                    "from either the params file or the BSEDict.".format(k)
                )

            # handle special cases where we need to expand arrays into multiple columns
            if k == 'natal_kick_array':
                initialbinarytable["natal_kick_array"] = [BSEDict['natal_kick_array']] * n
                for j, column_name in enumerate(NATAL_KICK_COLUMNS):
                    for sn in range(2):
                        col = f"{column_name}_{sn+1}"
                        if col in initialbinarytable.columns:
                            initialbinarytable[col] = BSEDict['natal_kick_array'][sn][j]
                        else:
                            new_cols[col] = BSEDict['natal_kick_array'][sn][j]

            elif k == 'qcrit_array':
                initialbinarytable["qcrit_array"] = [BSEDict['qcrit_array']] * n
                for kstar in range(16):
                    col = f"qcrit_{kstar}"
                    if col in initialbinarytable.columns:
                        initialbinarytable[col] = BSEDict['qcrit_array'][kstar]
                    else:
                        new_cols[col] = BSEDict['qcrit_array'][kstar]

            elif k == 'fprimc_array':
                initialbinarytable["fprimc_array"] = [BSEDict['fprimc_array']] * n
                for kstar in range(16):
                    col = f"fprimc_{kstar}"
                    if col in initialbinarytable.columns:
                        initialbinarytable[col] = BSEDict['fprimc_array'][kstar]
                    else:
                        new_cols[col] = BSEDict['fprimc_array'][kstar]
            elif k == 'alpha1':
                columns_values = [BSEDict['alpha1']] * len(initialbinarytable)
                initialbinarytable = initialbinarytable.assign(alpha1=columns_values)
                for kstar in range(0,2):
                    columns_values = pd.Series([BSEDict['alpha1'][kstar]] * len(initialbinarytable),
                                               index=initialbinarytable.index,
                                               name='alpha1_{0}'.format(kstar))
                    initialbinarytable.loc[:, 'alpha1_{0}'.format(kstar)] = columns_values
            elif k == 'acc_lim':
                columns_values = [BSEDict['acc_lim']] * len(initialbinarytable)
                initialbinarytable = initialbinarytable.assign(acc_lim=columns_values)
                for kstar in range(0,2):
                    columns_values = pd.Series([BSEDict['acc_lim'][kstar]] * len(initialbinarytable),
                                               index=initialbinarytable.index,
                                               name='acc_lim_{0}'.format(kstar))
                    initialbinarytable.loc[:, 'acc_lim_{0}'.format(kstar)] = columns_values
            else:
                # base case: if it's present, overwrite, if not, add to a list of new columns (see below)
                if k in initialbinarytable.columns:
                    initialbinarytable[k] = v
                else:
                    new_cols[k] = v

        # for columns that are new to the initial binary table, concat once
        if new_cols:
            new_df = pd.DataFrame(new_cols, index=idx)
            initialbinarytable = pd.concat([initialbinarytable, new_df], axis=1)



        # Here we perform two checks
        # First, if the BSE parameters are not in the initial binary table
        # and either a dictionary or an inifile was not provided
        # then we need to raise an ValueError and tell the user to provide
        # either a dictionary or an inifile or add more columns
        if BSEDict and SSEDict is None:
            if ((not set(INITIAL_BINARY_TABLE_SAVE_COLUMNS).issubset(initialbinarytable.columns)) and
               (not set(INITIAL_CONDITIONS_PASS_COLUMNS).issubset(initialbinarytable.columns))):
                raise ValueError("You are passing BSE parameters as columns in the "
                                 "initial binary table but not all BSE parameters are defined. "
                                 "Please pass a BSEDict or a params file or make sure "
                                 "you have all BSE parameters as columns {0} or {1}.".format(
                                  INITIAL_BINARY_TABLE_SAVE_COLUMNS, INITIAL_CONDITIONS_PASS_COLUMNS))
            
        if (BSEDict and not SSEDict) or (SSEDict and not BSEDict):
            raise ValueError("If you are passing BSE parameters as columns in the "
                             "initial binary table you must also pass SSE parameters "
                             "in the initial binary table.")

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

        if (pd.Series(ALPHA_COLUMNS).isin(initialbinarytable.keys()).all()) and ('alpha1' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(alpha1=initialbinarytable[ALPHA_COLUMNS].values.tolist())
        
        if (pd.Series(ACCLIM_COLUMNS).isin(initialbinarytable.keys()).all()) and ('acc_lim' not in BSEDict):
            initialbinarytable = initialbinarytable.assign(acc_lim=initialbinarytable[ACCLIM_COLUMNS].values.tolist())

        # update timesteps based on mass modifier
        mass_modifier_cols = set(['pts1', 'pts2', 'pts3']).difference(columns_in_passed_initC)
        if dt_mass_modifiers and len(mass_modifier_cols) != 0:
            # warn the user if their mass ranges overlap
            for i in range(len(dt_mass_modifiers)):
                for j in range(i + 1, len(dt_mass_modifiers)):
                    m_low_i, m_high_i, _ = dt_mass_modifiers[i]
                    m_low_j, m_high_j, _ = dt_mass_modifiers[j]
                    if (m_low_i < m_high_j) and (m_low_j < m_high_i):
                        overlap_range = (max(m_low_i, m_low_j), min(m_high_i, m_high_j))
                        warnings.warn(
                            f"Mass ranges for timestep modifiers overlap. You passed {dt_mass_modifiers[i]} "
                            f"and {dt_mass_modifiers[j]} which have overlapping mass ranges in {overlap_range}. "
                            f"This will result in *both* timestep modifiers being applied in the overlap region."
                            "If intentional, separate the overlap region into its own mass range with its "
                            "own modifier to avoid this warning."
                        )

            # apply the modifiers to the appropriate systems based on the primary mass, left->right
            for m_low, m_high, mod in dt_mass_modifiers:
                if mod <= 0:
                    raise ValueError(f"Timestep modifiers must be positive. You passed {mod} for the "
                                     f"mass range {m_low} to {m_high}.")
                mask = (initialbinarytable['mass_1'] >= m_low) & (initialbinarytable['mass_1'] < m_high)
                for col in mass_modifier_cols:
                    initialbinarytable.loc[mask, col] *= mod

        # if stellar engine is METISSE then check all of the SSE columns are present and if not raise an error
        if initialbinarytable['stellar_engine'].iloc[0] == 'metisse':
            if not set(INITIAL_CONDITIONS_SSE_COLUMN).issubset(initialbinarytable.columns):
                raise ValueError("If you want to use the METISSE stellar engine, you must provide the following in the SSEDict, initial binary table, or params file: path_to_tracks, path_to_he_tracks, z_accuracy_limit.")
        else:
            # if not using METISSE, set default values for the SSE columns if they are not present in the initial binary table
            if 'stellar_engine' not in initialbinarytable.columns:
                initialbinarytable = initialbinarytable.assign(stellar_engine='sse')
            if 'path_to_tracks' not in initialbinarytable.columns:
                initialbinarytable = initialbinarytable.assign(path_to_tracks='')
            if 'path_to_he_tracks' not in initialbinarytable.columns:
                initialbinarytable = initialbinarytable.assign(path_to_he_tracks='')
            if 'z_accuracy_limit' not in initialbinarytable.columns:
                initialbinarytable = initialbinarytable.assign(z_accuracy_limit=1e-2)

        # need to ensure that the order of parameters that we pass to BSE is correct
        initial_conditions = initialbinarytable[INITIAL_CONDITIONS_PASS_COLUMNS].to_dict('records')

        # ensure that metallicity is in the valid range (Z in [1e-4, 0.03])
        low_met_mask = (initialbinarytable["metallicity"] < 1e-4)
        high_met_mask = (initialbinarytable["metallicity"] > 0.03)
        if any(low_met_mask | high_met_mask) and not initialbinarytable["stellar_engine"].values[0]:
            raise ValueError(
                f"COSMIC-SSE only supports metallicities in the range [1e-4, 0.03]. You have {sum(low_met_mask)} "
                f"systems with metallicity below 1e-4 and {sum(high_met_mask)} systems with metallicity "
                "above 0.03. Some examples of problematic binaries have the following bin_nums: "
                f"{initialbinarytable['bin_num'][low_met_mask | high_met_mask].values[:5]}."
            )
        
        # ensure that the initial masses are in the valid range for the loaded tracks
        if initialbinarytable["stellar_engine"].values[0] == "metisse":
            low_mass_mask = (initialbinarytable["mass_1"] < m_min) | (initialbinarytable["mass_2"] < m_min)
            high_mass_mask = (initialbinarytable["mass_1"] > m_max) | (initialbinarytable["mass_2"] > m_max)
            if any(low_mass_mask | high_mass_mask):
                raise ValueError(
                    f"COSMIC-METISSE only supports initial masses in the range specified by the loaded tracks [{m_min}, {m_max}]. You have {sum(low_mass_mask)} "
                    f"systems with mass below {m_min} and {sum(high_mass_mask)} systems with mass above {m_max}. "
                "Some examples of problematic binaries have the following bin_nums: "
                f"{initialbinarytable['bin_num'][low_mass_mask | high_mass_mask].values[:5]}."
            )

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

        # evolve one system to get zpars
        _, _, _, _, _, zpars = _evolve_single_system(initial_conditions[0], None)

        # helper to collect results with an optional tqdm progress bar
        def _collect(pool, func, items, total=None):
            if progress:
                return list(tqdm.tqdm(pool.imap(func, items), total=total, desc='Evolving', unit='sys'))
            return list(pool.map(func, items))

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
                    output = _collect(pool, _evolve_multi_system, initial_conditions_blocked,
                                      total=len(initial_conditions_blocked))
                else:
                    evolve_args = partial(_evolve_single_system, zpars=zpars)
                    output = _collect(pool, evolve_args, initial_conditions,
                                      total=len(initial_conditions))
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
                output = _collect(pool, _evolve_multi_system, initial_conditions_blocked,
                                  total=len(initial_conditions_blocked))
            else:
                evolve_args = partial(_evolve_single_system, zpars=zpars)
                output = _collect(pool, evolve_args, initial_conditions,
                                  total=len(initial_conditions))

        output = np.array(output, dtype=object)
        bpp_arrays = np.vstack(output[:, 1])
        bcm_arrays = np.vstack(output[:, 2])
        kick_info_arrays = np.vstack(output[:, 3])

        natal_kick_arrays = np.vstack(output[:, 4])
        natal_kick_arrays = natal_kick_arrays.reshape(-1, 1, len(FLATTENED_NATAL_KICK_COLUMNS))

        # update initial table with sampled kicks
        to_add = {}
        for idx, column in enumerate(FLATTENED_NATAL_KICK_COLUMNS):
            if column not in initialbinarytable.columns:
                to_add[column] = natal_kick_arrays[:, 0, idx]
            else:
                initialbinarytable[column] = natal_kick_arrays[:, 0, idx]

        # if kicks weren't already present, add them
        if to_add:
            natal_kick_df = pd.DataFrame(to_add, index=initialbinarytable.index)
            initialbinarytable = pd.concat([initialbinarytable, natal_kick_df], axis=1)

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


def _evolve_single_system(f, zpars=None):
    if zpars is None:
        zpars = np.zeros(20, dtype=float)
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
        _evolvebin.ceflags.ce2stageflag = f["ce2stageflag"]
        _evolvebin.ceflags.ceflag = f["ceflag"]
        _evolvebin.flags.tflag = f["tflag"]
        _evolvebin.flags.ifflag = f["ifflag"]
        _evolvebin.flags.wdflag = f["wdflag"]
        _evolvebin.flags.rtmsflag = f["rtmsflag"]
        _evolvebin.snvars.pisn = f["pisn"]
        _evolvebin.snvars.ppi_co_shift = f["ppi_co_shift"]
        _evolvebin.snvars.ppi_extra_ml = f["ppi_extra_ml"]
        _evolvebin.flags.bhflag = f["bhflag"]
        _evolvebin.flags.remnantflag = f["remnantflag"]
        _evolvebin.flags.maltsev_mode = f["maltsev_mode"]
        _evolvebin.snvars.maltsev_fallback = f["maltsev_fallback"]
        _evolvebin.snvars.maltsev_pf_prob = f["maltsev_pf_prob"]
        _evolvebin.snvars.fryer_mass_limit = f["fryer_mass_limit"]
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
        _evolvebin.snvars.fryer_fmix = f["fryer_fmix"]
        _evolvebin.snvars.fryer_mcrit_nsbh = f["fryer_mcrit_nsbh"]
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
        _evolvebin.mtvars.smt_periastron_check = f["smt_periastron_check"]
        _evolvebin.windvars.beta = f["beta"]
        _evolvebin.windvars.xi = f["xi"]
        _evolvebin.windvars.acc2 = f["acc2"]
        _evolvebin.windvars.epsnov = f["epsnov"]
        _evolvebin.windvars.eddfac = f["eddfac"]
        _evolvebin.windvars.gamma = f["gamma"]
        _evolvebin.windvars.lbv_flag = f["LBV_flag"]
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
        _evolvebin.snvars.mm_mu_ns = f["mm_mu_ns"]
        _evolvebin.snvars.mm_mu_bh = f["mm_mu_bh"]
        _evolvebin.cmcpass.using_cmc = 0
        
        if f["stellar_engine"] == "sse":
            _evolvebin.se_flags.using_sse = 1
            _evolvebin.se_flags.using_metisse = 0
            _evolvebin.metissevars.path_to_tracks = ""
            _evolvebin.metissevars.path_to_he_tracks = ""
            _evolvebin.metissevars.z_match_limit = f["z_accuracy_limit"]
            _evolvebin.metissevars.METISSE_verbose = False
        elif f["stellar_engine"] == "metisse":
            _evolvebin.se_flags.using_metisse = 1
            _evolvebin.se_flags.using_sse = 0
            _evolvebin.metissevars.path_to_tracks = f["path_to_tracks"]
            _evolvebin.metissevars.path_to_he_tracks = f["path_to_he_tracks"]
            _evolvebin.metissevars.z_match_limit = f["z_accuracy_limit"]
            _evolvebin.metissevars.METISSE_verbose = False
        else:
            raise ValueError("Use either 'sse' or 'metisse' as stellar engine")

        _evolvebin.col.n_col_bpp = f["n_col_bpp"]
        _evolvebin.col.col_inds_bpp = f["col_inds_bpp"]
        _evolvebin.col.n_col_bcm = f["n_col_bcm"]
        _evolvebin.col.col_inds_bcm = f["col_inds_bcm"]

        [zpars, kick_info, bpp_index, bcm_index] = _evolvebin.evolv2([f["kstar_1"], f["kstar_2"]],
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
                                                              zpars,
                                                              f["kick_info"])
        if bpp_index<0:
            raise ValueError("Failed in METISSE_zcnsts")
        else:
            bpp = _evolvebin.binary.bpp[:bpp_index, :f["n_col_bpp"]].copy()
            _evolvebin.binary.bpp[:bpp_index, :f["n_col_bpp"]] = np.zeros(bpp.shape)
            bcm = _evolvebin.binary.bcm[:bcm_index, :f["n_col_bcm"]].copy()
            _evolvebin.binary.bcm[:bcm_index, :f["n_col_bcm"]] = np.zeros(bcm.shape)

            bpp = np.hstack((bpp, np.ones((bpp.shape[0], 1))*f["bin_num"]))
            bcm = np.hstack((bcm, np.ones((bcm.shape[0], 1))*f["bin_num"]))
            kick_info = np.hstack((kick_info, np.ones((kick_info.shape[0], 1))*f["bin_num"]))

        return f, bpp, bcm, kick_info, _evolvebin.snvars.natal_kick_array.copy(), zpars

    except Exception as e:
        print(e)
        raise


def _evolve_multi_system(f):
    try:
        zpars = np.zeros(20, dtype=float)
        res_bcm = np.zeros(f.shape[0], dtype=object)
        res_bpp = np.zeros(f.shape[0], dtype=object)
        res_kick_info = np.zeros(f.shape[0], dtype=object)
        res_natal_kick_array = np.zeros(f.shape[0], dtype=object)
        for i in range(0, f.shape[0]):

            # call evolve single system
            _, bpp, bcm, kick_info, _, zpars = _evolve_single_system(f[i], zpars=zpars)

            # add results to pre-allocated list
            res_bpp[i] = bpp
            res_bcm[i] = bcm
            res_kick_info[i] = kick_info
            res_natal_kick_array[i] = _evolvebin.snvars.natal_kick_array

        return f, np.vstack(res_bpp), np.vstack(res_bcm), np.vstack(res_kick_info), np.vstack(res_natal_kick_array)

    except Exception as e:
        print(e)
        raise

def read_tracks_for_METISSE(path_to_tracks,IBT_Z,z_accuracy_limit,is_he):

    """load in the metallicity, format, and eep files
    
    Parameters
    ----------
    path_to_tracks : str
        Direct path to where all single star data and metallicty/format files are stored
        for hydrogen-rich stars
    IBT_Z : float
        Metallicity from the initialbinarytable
    z_accuracy_limit : float
        Tolerance for match in metallicity value
    is_he : bool, default=False
        Indicates whether the tracks are helium-enriched.
    Returns
    -------
    m_min : float
        Minimum mass in the tracks that were loaded in
    m_max : float
        Maximum mass in the tracks that were loaded in
    
    """

    # load in the METISSE files
    met_files = utils.get_METISSE_metallicity_files(path_to_tracks)

    # loop over the hydrogen metallicity files to find the one that is the closest
    # to the metallicity in the initial binary table
    met_dict_keep = None
    fmt_dict_keep = None
    mets = []
    for i, m in enumerate(met_files):
        met_dict = utils.read_metallicity_file(m)
        mets.append(met_dict['z_files'])
        if abs(met_dict['z_files'] - IBT_Z)/min(met_dict['z_files'], IBT_Z) <= z_accuracy_limit:
            met_dict_keep = met_dict
            fmt_dict_keep = utils.read_format_file(met_dict['format_file'])
            Z_idx = i
    if met_dict_keep is None: 
        raise ValueError("No metallicity file found that matches the metallicity "
                         "in the initial binary table. Please check the metallicity "
                         "and supply one that is in this list: {0}".format(mets))

    assert fmt_dict_keep is not None 


    if is_he:
        # Pass the format dictionaries for helium tracks:
        _evolvebin.c_m_interface.set_format_controls_he(
            read_eep=fmt_dict_keep['read_eep_files'],
            bgb=fmt_dict_keep['bgb_eep'],
            cheburn=fmt_dict_keep['cheburn_eep'],
            ta_cheb=fmt_dict_keep['ta_cheb_eep'],
            tpagb=fmt_dict_keep['tpagb_eep'],
            ccburn=fmt_dict_keep['ccburn_eep'],
            postagb=fmt_dict_keep['post_agb_eep'],
            initeep=fmt_dict_keep['initial_eep'],
            finaleep=fmt_dict_keep['final_eep'],
            fixtrack=fmt_dict_keep['fix_track'],
            loweep=fmt_dict_keep['low_mass_final_eep'],
            higheep=fmt_dict_keep['high_mass_final_eep'],
            age_col=fmt_dict_keep['age_colname'],
            mass_col=fmt_dict_keep['mass_colname'],
            logl_col=fmt_dict_keep['log_l_colname'],
            logt_col=fmt_dict_keep['log_t_colname'],
            logr_col=fmt_dict_keep['log_r_colname'],
            he_mass_col=fmt_dict_keep['he_core_mass'],
            co_mass_col=fmt_dict_keep['co_core_mass'],
            he_radius_col=fmt_dict_keep['he_core_radius'],
            co_radius_col=fmt_dict_keep['co_core_radius'],
            mass_env_col=fmt_dict_keep['mass_conv_envelope'],
            radius_env_col=fmt_dict_keep['radius_conv_envelope'],
            logtc_col=fmt_dict_keep['log_tc'],
            he4_col=fmt_dict_keep['he4_mass_frac'],
            c12_col=fmt_dict_keep['c12_mass_frac'],
            o16_col=fmt_dict_keep['o16_mass_frac']
        )
    else:
        # Pass the format dictionaries:
        _evolvebin.c_m_interface.set_format_controls_h(
            read_eep=fmt_dict_keep['read_eep_files'],
            prems=fmt_dict_keep['prems_eep'],
            zams=fmt_dict_keep['zams_eep'],
            iams=fmt_dict_keep['iams_eep'],
            tams=fmt_dict_keep['tams_eep'],
            bgb=fmt_dict_keep['bgb_eep'],
            cheign=fmt_dict_keep['cheignition_eep'],
            cheburn=fmt_dict_keep['cheburn_eep'],
            ta_cheb=fmt_dict_keep['ta_cheb_eep'],
            tpagb=fmt_dict_keep['tpagb_eep'],
            ccburn=fmt_dict_keep['ccburn_eep'],
            postagb=fmt_dict_keep['post_agb_eep'],
            initeep=fmt_dict_keep['initial_eep'],
            finaleep=fmt_dict_keep['final_eep'],
            fixtrack=fmt_dict_keep['fix_track'],
            loweep=fmt_dict_keep['low_mass_final_eep'],
            higheep=fmt_dict_keep['high_mass_final_eep'],
            age_col=fmt_dict_keep['age_colname'],
            mass_col=fmt_dict_keep['mass_colname'],
            logl_col=fmt_dict_keep['log_l_colname'],
            logt_col=fmt_dict_keep['log_t_colname'],
            logr_col=fmt_dict_keep['log_r_colname'],
            he_mass_col=fmt_dict_keep['he_core_mass'],
            co_mass_col=fmt_dict_keep['co_core_mass'],
            he_radius_col=fmt_dict_keep['he_core_radius'],
            co_radius_col=fmt_dict_keep['co_core_radius'],
            mass_env_col=fmt_dict_keep['mass_conv_envelope'],
            radius_env_col=fmt_dict_keep['radius_conv_envelope'],
            logtc_col=fmt_dict_keep['log_tc'],
            he4_col=fmt_dict_keep['he4_mass_frac'],
            c12_col=fmt_dict_keep['c12_mass_frac'],
            o16_col=fmt_dict_keep['o16_mass_frac']
        )

    # Finally, load in the EEPs!
    track_list = utils.read_eep_directory(
                met_dict_keep['eep_tracks_dir'],
                fmt_dict_keep)
    m_min, m_max = populate_tracks(track_list, is_he)
    return m_min, m_max


def populate_tracks(track_list, is_he=False):
    """
    Populate Fortran track data structures from a list of Python track dictionaries
    and pass them to the COSMIC Fortran backend.

    Parameters
    ----------
    track_list : list of dict
        Each dictionary must contain the following keys:
            - 'filename' : str
            - 'initial_mass' : float
            - 'initial_Y' : float
            - 'initial_Z' : float
            - 'Fe_div_H' : float
            - 'alpha_div_Fe' : float
            - 'v_div_vcrit' : float
            - 'ntrack' : int
            - 'neep' : int
            - 'ncol' : int
            - 'eep' : array-like of shape (neep,)
            - 'tr' : array-like of shape (ncol, ntrack)
            - 'cols' : list of str of length ncol

    is_he : bool, default=False
        Indicates whether the tracks are helium-enriched.

    Returns
    -------
    tuple of float
        The minimum and maximum initial masses from the track list.

        The function calls the Fortran subroutine `_evolvebin.c_m_interface.set_tracks_from_python`
        and populates the Fortran-side track arrays. The Python-side arrays are used only
        as temporary buffers for the call.
    """
    ntracks = len(track_list)
    col_width = 32  # must match Fortran CHARACTER(len=32)

    # Allocate arrays
    filenames = np.array([t['filename'].encode('ascii') for t in track_list], dtype='S256')
    initial_mass = np.array([t['initial_mass'] for t in track_list], dtype=np.float64)
    initial_Y = np.array([t['initial_Y'] for t in track_list], dtype=np.float64)
    initial_Z = np.array([t['initial_Z'] for t in track_list], dtype=np.float64)
    Fe_div_H = np.array([t['Fe_div_H'] for t in track_list], dtype=np.float64)
    alpha_div_Fe = np.array([t['alpha_div_Fe'] for t in track_list], dtype=np.float64)
    v_div_vcrit = np.array([t['v_div_vcrit'] for t in track_list], dtype=np.float64)
    ntrack_arr = np.array([t['ntrack'] for t in track_list], dtype=np.int32)
    neep_arr = np.array([t['neep'] for t in track_list], dtype=np.int32)
    ncol_arr = np.array([t['ncol'] for t in track_list], dtype=np.int32)

    # Determine max sizes
    max_neep = max(neep_arr)
    max_ncol = max(ncol_arr)
    max_ntrack = sum(ntrack_arr)

    # Prepare 2D arrays

    if max_neep<0:
        max_neep = 0
    eep_data = np.zeros((max_neep, ntracks), dtype=np.int32, order='F')
    tr_data = np.zeros((max_ncol, max_ntrack), dtype=np.float64, order='F')
    col_names = np.full((max_ncol, ntracks), b' ' * col_width, dtype=f'S{col_width}', order='F')

    # Fill arrays
    offset = 0
    for i, t in enumerate(track_list):
        if max_neep>0:
            eep_data[:t['neep'], i] = t['eep']
        tr_data[:t['ncol'], offset:offset+t['ntrack']] = t['tr']

        for j, col in enumerate(t['cols']):
            s = col.encode('ascii')[:col_width]      # truncate if too long
            col_names[j, i] = s.ljust(col_width, b' ')  # pad with spaces
        offset += t['ntrack']

    _evolvebin.c_m_interface.set_tracks_from_python(
        filenames, initial_mass, initial_Y, initial_Z,
        Fe_div_H, alpha_div_Fe, v_div_vcrit,
        ntrack_arr, neep_arr, ncol_arr, 
        eep_data, tr_data, col_names, is_he
    )

    return np.min(initial_mass), np.max(initial_mass) 


    
    
