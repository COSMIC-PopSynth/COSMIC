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

"""`evolve`
"""

import numpy as np
from gwpy.utils import mp as mp_utils
from cosmic import _evolvebin
import pandas as pd
from astropy.table import Table

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['Evolve']


bpp_columns = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2' ,
               'sep', 'porb', 'ecc', 'RROL_1', 'RROL_2', 'evol_type',
               'Vsys_1', 'Vsys_2', 'SNkick', 'SNtheta', 'bin_num']

bcm_columns = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lumin_1', 'rad_1',
               'teff_1', 'massc_1', 'radc_1', 'menv_1', 'renv_1', 'epoch_1',
               'ospin_1', 'deltam_1', 'RROL_1', 'kstar_2', 'mass0_2', 'mass_2',
               'lumin_2', 'rad_2', 'teff_2', 'massc_2', 'radc_2', 'menv_2',
               'renv_2', 'epoch_2', 'ospin_2', 'deltam_2', 'RROL_2',
               'porb', 'sep', 'ecc', 'B_0_1', 'B_0_2',
               'SNkick_1', 'SNkick_2', 'Vsys_final', 'SNtheta_final',
               'SN_1', 'SN_2', 'bin_state', 'merger_type', 'bin_num']

class Evolve(Table):
    def __init__():
        '''
        initialize Evolve
        '''

    @classmethod
    def evolve(cls, initialbinarytable, BSEDict, **kwargs):
        """After setting a number of initial conditions
        we evolve the system.

        Parameters
        ----------
        nproc : `int`, optional, default: 1
            number of CPUs to use to evolve systems
            in parallel
        idx : `int`, optional, default: 0
            initial index of the bcm/bpp arrays
        dtp : `float`, optional: default: tphysf
            timestep size in Myr for bcm output where tphysf
            is total evolution time in Myr

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

        if 'neta' not in initialbinarytable.keys():
            initialbinarytable['neta'] = BSEDict['neta']
        if 'bwind' not in initialbinarytable.keys():
            initialbinarytable['bwind'] = BSEDict['bwind']
        if 'hewind' not in initialbinarytable.keys():
            initialbinarytable['hewind'] = BSEDict['hewind']
        if 'alpha1' not in initialbinarytable.keys():
            initialbinarytable['alpha1'] = BSEDict['alpha1']
        if 'lambdaf' not in initialbinarytable.keys():
            initialbinarytable['lambdaf'] = BSEDict['lambdaf']
        if 'cekickflag' not in initialbinarytable.keys():
            initialbinarytable['cekickflag'] = BSEDict['cekickflag']
        if 'cemergeflag' not in initialbinarytable.keys():
            initialbinarytable['cemergeflag'] = BSEDict['cemergeflag']
        if 'cehestarflag' not in initialbinarytable.keys():
            initialbinarytable['cehestarflag'] = BSEDict['cehestarflag']
        if 'ceflag' not in initialbinarytable.keys():
            initialbinarytable['ceflag'] = BSEDict['ceflag']
        if 'tflag' not in initialbinarytable.keys():
            initialbinarytable['tflag'] = BSEDict['tflag']
        if 'ifflag' not in initialbinarytable.keys():
            initialbinarytable['ifflag'] = BSEDict['ifflag']
        if 'wdflag' not in initialbinarytable.keys():
            initialbinarytable['wdflag'] = BSEDict['wdflag']
        if 'ppsn' not in initialbinarytable.keys():
            initialbinarytable['ppsn'] = BSEDict['ppsn']
        if 'bhflag' not in initialbinarytable.keys():
            initialbinarytable['bhflag'] = BSEDict['bhflag']
        if 'nsflag' not in initialbinarytable.keys():
            initialbinarytable['nsflag'] = BSEDict['nsflag']
        if 'mxns' not in initialbinarytable.keys():
            initialbinarytable['mxns'] = BSEDict['mxns']
        if 'pts1' not in initialbinarytable.keys():
            initialbinarytable['pts1'] = BSEDict['pts1']
        if 'pts2' not in initialbinarytable.keys():
            initialbinarytable['pts2'] = BSEDict['pts2']
        if 'pts3' not in initialbinarytable.keys():
            initialbinarytable['pts3'] = BSEDict['pts3']
        if 'ecsnp' not in initialbinarytable.keys():
            initialbinarytable['ecsnp'] = BSEDict['ecsnp']
        if 'ecsn_mlow' not in initialbinarytable.keys():
            initialbinarytable['ecsn_mlow'] = BSEDict['ecsn_mlow']
        if 'aic' not in initialbinarytable.keys():
            initialbinarytable['aic'] = BSEDict['aic']
        if 'sigma' not in initialbinarytable.keys():
            initialbinarytable['sigma'] = BSEDict['sigma']
        if 'sigmadiv' not in initialbinarytable.keys():
            initialbinarytable['sigmadiv'] = BSEDict['sigmadiv']
        if 'bhsigmafrac' not in initialbinarytable.keys():
            initialbinarytable['bhsigmafrac'] = BSEDict['bhsigmafrac']
        if 'polar_kick_angle' not in initialbinarytable.keys():
            initialbinarytable['polar_kick_angle'] = BSEDict['polar_kick_angle']
        if 'beta' not in initialbinarytable.keys():
            initialbinarytable['beta'] = BSEDict['beta']
        if 'xi' not in initialbinarytable.keys():
            initialbinarytable['xi'] = BSEDict['xi']
        if 'acc2' not in initialbinarytable.keys():
            initialbinarytable['acc2'] = BSEDict['acc2']
        if 'epsnov' not in initialbinarytable.keys():
            initialbinarytable['epsnov'] = BSEDict['epsnov']
        if 'eddfac' not in initialbinarytable.keys():
            initialbinarytable['eddfac'] = BSEDict['eddfac']
        if 'gamma' not in initialbinarytable.keys():
            initialbinarytable['gamma'] = BSEDict['gamma']
        if 'bconst' not in initialbinarytable.keys():
            initialbinarytable['bconst'] = BSEDict['bconst']
        if 'ck' not in initialbinarytable.keys():
            initialbinarytable['ck'] = BSEDict['ck']
        if 'merger' not in initialbinarytable.keys():
            initialbinarytable['merger'] = BSEDict['merger']
        if 'windflag' not in initialbinarytable.keys():
            initialbinarytable['windflag'] = BSEDict['windflag']
        if 'dtp' not in initialbinarytable.keys():
            initialbinarytable['dtp'] = kwargs.pop('dtp', initialbinarytable['tphysf'])
        if 'randomseed' not in initialbinarytable.keys():
            initialbinarytable['randomseed'] = np.random.randint(1, 1000000, size=len(initialbinarytable))
        if 'bin_num' not in initialbinarytable.keys():
            initialbinarytable['bin_num'] = np.arange(idx, idx + len(initialbinarytable))


        natal_kick_columns = ['SNkick_1', 'SNkick_2', 'phi_1', 'phi_2', 'theta_1', 'theta_2']
        if pd.Series(natal_kick_columns).isin(initialbinarytable.keys()).all() and 'natal_kick_array' not in initialbinarytable.keys():
            initialbinarytable['natal_kick_array'] = initialbinarytable[natal_kick_columns].values.tolist()
        if 'natal_kick_array' not in initialbinarytable.keys():
            initialbinarytable['natal_kick_array'] = [BSEDict['natal_kick_array']] * len(initialbinarytable)
            for idx, column_name in enumerate(natal_kick_columns):
                initialbinarytable.loc[:, column_name] = pd.Series([BSEDict['natal_kick_array'][idx]] * len(initialbinarytable), index=initialbinarytable.index, name=column_name)

        qcrit_columns = ['qcrit_{0}'.format(kstar) for kstar in range(0,16)]
        if pd.Series(qcrit_columns).isin(initialbinarytable.keys()).all() and 'qcrit_array' not in initialbinarytable.keys():
            initialbinarytable['qcrit_array'] = initialbinarytable[qcrit_columns].values.tolist()

        if 'qcrit_array' not in initialbinarytable.keys():
            initialbinarytable['qcrit_array'] = [BSEDict['qcrit_array']] * len(initialbinarytable)
            for kstar in range(0,16):
                initialbinarytable.loc[:, 'qcrit_{0}'.format(kstar)] = pd.Series([BSEDict['qcrit_array'][kstar]]* len(initialbinarytable), index=initialbinarytable.index, name='qcrit_{0}'.format(kstar))

        # need to ensure that the order of variables is correct
        initial_conditions = initialbinarytable[['kstar_1', 'kstar_2', 'mass1_binary', 'mass2_binary', 'porb', 'ecc',
                                                'metallicity', 'tphysf', 'neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                                                'ceflag', 'tflag', 'ifflag', 'wdflag', 'ppsn', 'bhflag', 'nsflag',
                                                'cekickflag', 'cemergeflag', 'cehestarflag',
                                                'mxns', 'pts1', 'pts2', 'pts3',
                                                'ecsnp', 'ecsn_mlow', 'aic', 'sigma', 'sigmadiv', 'bhsigmafrac', 'polar_kick_angle',
                                                'natal_kick_array', 'qcrit_array',
                                                'beta', 'xi', 'acc2', 'epsnov',
                                                'eddfac', 'gamma', 'bconst', 'ck', 'merger', 'windflag', 'dtp',
                                                'randomseed', 'bin_num']].values

        initial_binary_table_column_names = ['kstar_1', 'kstar_2', 'mass1_binary', 'mass2_binary', 'porb', 'ecc',
                                             'metallicity', 'tphysf', 'neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
                                             'ceflag', 'tflag', 'ifflag', 'wdflag', 'ppsn', 'bhflag', 'nsflag',
                                             'cekickflag', 'cemergeflag', 'cehestarflag',
                                             'mxns', 'pts1', 'pts2', 'pts3',
                                             'ecsnp', 'ecsn_mlow', 'aic', 'sigma', 'sigmadiv', 'bhsigmafrac', 'polar_kick_angle',
                                             'beta', 'xi', 'acc2', 'epsnov',
                                             'eddfac', 'gamma', 'bconst', 'ck', 'merger', 'windflag', 'dtp',
                                             'randomseed', 'bin_num']

        initial_binary_table_column_names.extend(natal_kick_columns)
        initial_binary_table_column_names.extend(qcrit_columns)

        initialbinarytable = initialbinarytable[initial_binary_table_column_names]

        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                [bpp, bcm] = _evolvebin.evolv2(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9],
                                               f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17], f[18], f[19],
                                               f[20], f[21], f[22], f[23], f[24], f[25], f[26], f[27], f[28], f[29],
                                               f[30], f[31], f[32], f[33], f[34], f[35], f[36], f[37], f[38], f[39],
                                               f[40], f[41], f[42], f[43], f[44], f[45], f[46], f[47])

                try:
                    bpp = bpp[:np.argwhere(bpp[:,0] == -1)[0][0]]
                    bcm = bcm[:np.argwhere(bcm[:,0] == -1)[0][0]]
                except IndexError:
                    bpp = bpp[:np.argwhere(bpp[:,0] > 0)[0][0]]
                    bcm = bcm[:np.argwhere(bcm[:,0] > 0)[0][0]]
                    raise Warning('bpp overload: mass1 = {0}, mass2 = {1}, porb = {2}, ecc = {3}, tphysf = {4}, metallicity = {5}'\
                                   .format(f[2], f[3], f[4], f[5], f[7], f[6]))

                bpp_bin_numbers = np.atleast_2d(np.array([f[48]] * len(bpp))).T
                bcm_bin_numbers = np.atleast_2d(np.array([f[48]] * len(bcm))).T

                bpp = np.hstack((bpp, bpp_bin_numbers))
                bcm = np.hstack((bcm, bcm_bin_numbers))

                return f, bpp, bcm

            except Exception as e:
                raise

        # evolve sysyems
        output = mp_utils.multiprocess_with_queues(
            nproc, _evolve_single_system, initial_conditions, raise_exceptions=False)

        output = np.array(output)
        bpp_arrays = np.vstack(output[:, 1])
        bcm_arrays = np.vstack(output[:, 2])

        bpp = pd.DataFrame(bpp_arrays,
                           columns=bpp_columns,
                           index=bpp_arrays[:, -1].astype(int))

        bcm = pd.DataFrame(bcm_arrays,
                           columns=bcm_columns,
                           index=bcm_arrays[:, -1].astype(int))

        bcm.merger_type = bcm.merger_type.astype(int).astype(str).apply(lambda x: x.zfill(4))
        bcm.bin_state = bcm.bin_state.astype(int)
        bpp.bin_num = bpp.bin_num.astype(int)
        bcm.bin_num = bcm.bin_num.astype(int)

        return bpp, bcm, initialbinarytable
