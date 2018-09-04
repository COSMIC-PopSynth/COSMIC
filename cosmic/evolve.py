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


bpp_columns = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2' , 'sep', 'ecc', 'RROL_1', 'RROL_2', 'evol_type']
bcm_columns = ['tphys', 'kstar_1', 'mass0_1', 'mass_1', 'lumin_1', 'rad_1', 'teff_1', 'massc_1',
'radc_1', 'menv_1', 'renv_1', 'epoch_1', 'ospin_1', 'deltam_1', 'RROL_1', 'kstar_2', 'mass0_2', 'mass_2',
'lumin_2', 'rad_2', 'teff_2', 'massc_2', 'radc_2', 'menv_2', 'renv_2', 'epoch_2', 'ospin_2', 'deltam_2', 'RROL_2',
'porb', 'sep', 'ecc', 'B_0_1', 'B_0_2', 'formation_1', 'formation_2']

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
            number of CPUs to use for parallel file reading

        kwargs: 

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

        initialbinarytable['neta'] = BSEDict['neta']
        initialbinarytable['bwind'] = BSEDict['bwind']
        initialbinarytable['hewind'] = BSEDict['hewind']
        initialbinarytable['alpha1'] = BSEDict['alpha1']
        initialbinarytable['lambdaf'] = BSEDict['lambdaf']
        initialbinarytable['ceflag'] = BSEDict['ceflag']
        initialbinarytable['tflag'] = BSEDict['tflag']
        initialbinarytable['ifflag'] = BSEDict['ifflag']
        initialbinarytable['wdflag'] = BSEDict['wdflag']
        initialbinarytable['bhflag'] = BSEDict['bhflag']
        initialbinarytable['nsflag'] = BSEDict['nsflag']
        initialbinarytable['mxns'] = BSEDict['mxns']
        initialbinarytable['pts1'] = BSEDict['pts1']
        initialbinarytable['pts2'] = BSEDict['pts2']
        initialbinarytable['pts3'] = BSEDict['pts3']
        initialbinarytable['sigma'] = BSEDict['sigma']
        initialbinarytable['beta'] = BSEDict['beta']
        initialbinarytable['xi'] = BSEDict['xi']
        initialbinarytable['acc2'] = BSEDict['acc2']
        initialbinarytable['epsnov'] = BSEDict['epsnov']
        initialbinarytable['eddfac'] = BSEDict['eddfac']
        initialbinarytable['gamma'] = BSEDict['gamma']
        initialbinarytable['bconst'] = BSEDict['bconst']
        initialbinarytable['CK'] = BSEDict['CK']
        initialbinarytable['merger'] = BSEDict['merger']
        initialbinarytable['windflag'] = BSEDict['windflag']
        initialbinarytable['dtp'] = initialbinarytable['tphysf']
        initialbinarytable['bin_num'] = np.arange(idx, idx + len(initialbinarytable))
        initialbinarytable['randomseed'] = np.random.randint(1, 1000000, size=len(initialbinarytable))

        initial_conditions = np.array(initialbinarytable) 

        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                [tmp1, tmp2] = _evolvebin.evolv2(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9],
                                        f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17], f[18], f[19],
                                        f[20], f[21], f[22], f[23], f[24], f[25], f[26], f[27], f[28], f[29],
                                        f[30], f[31], f[32], f[33], f[34], f[36])

                bpp_tmp = tmp1[np.argwhere(tmp1[:,0]>0),:].squeeze(1)
                bcm_tmp = tmp2[np.argwhere(tmp2[:,0]>1),:].squeeze(1)

                bpp_tmp = pd.DataFrame(bpp_tmp, columns=bpp_columns, index=[int(f[35])] * len(bpp_tmp))
                bpp_tmp['bin_num'] = int(f[35])
                bpp_tmp.set_index('bin_num')

                bcm_tmp = pd.DataFrame(bcm_tmp, columns=bcm_columns, index=[int(f[35])] * len(bcm_tmp))
                bcm_tmp['bin_num'] = int(f[35])
                bcm_tmp.set_index('bin_num')

                return f, bpp_tmp, bcm_tmp

            except Exception as e:
                if nproc == 1:
                    raise
                else:
                    return f, e
        # evolve sysyems
        output = mp_utils.multiprocess_with_queues(
            nproc, _evolve_single_system, initial_conditions, raise_exceptions=False)

        # raise exceptions (from multiprocessing, single process raises inline)
        for f, x, y in output:
            if isinstance(x, Exception):
                x.args = ('Failed to evolve %s: %s' % (f, str(x)),)
                raise x
            if isinstance(y, Exception):
                y.args = ('Failed to evolve %s: %s' % (f, str(y)),)
                raise y

        output_bpp = pd.DataFrame()
        output_bcm = pd.DataFrame()
        for f, x, y in output:
            output_bpp = output_bpp.append(x)
            output_bcm = output_bcm.append(y)        

        initialbinarytable.set_index('bin_num')
        return output_bpp, output_bcm, initialbinarytable 
