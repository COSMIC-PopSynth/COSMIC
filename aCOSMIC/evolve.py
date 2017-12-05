# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of aCOSMIC.
#
# aCOSMIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# aCOSMIC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with aCOSMIC.  If not, see <http://www.gnu.org/licenses/>.

"""`sample`
"""

import numpy as np
from gwpy.utils import mp as mp_utils
from aCOSMIC import _evolvebin

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'Evolve'

class Evolve:
    def __init__(self, sample):
        '''
        initialize Evolve
        '''
        self.initial_conditions = sample


    def evolve(self, **kwargs):
        """After setting a number of initial conditions
        we evolve the system.

        Parameters
        ----------
        nproc : `int`, optional, default: 1
            number of CPUs to use for parallel file reading

        kwargs: 

        Returns
        -------
        An evolved binary
        """
        # Populate BSEDict
        BSEDict = kwargs.pop('BSEDict')
        self.initial_conditions.neta = np.repeat(BSEDict['neta'], self.initial_conditions.kstar1.size)
        self.initial_conditions.bwind = np.repeat(BSEDict['bwind'], self.initial_conditions.kstar1.size)
        self.initial_conditions.hewind = np.repeat(BSEDict['hewind'], self.initial_conditions.kstar1.size)
        self.initial_conditions.alpha1 = np.repeat(BSEDict['alpha1'], self.initial_conditions.kstar1.size)
        self.initial_conditions.lambdaf = np.repeat(BSEDict['lambdaf'], self.initial_conditions.kstar1.size)
        self.initial_conditions.ceflag = np.repeat(BSEDict['ceflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.tflag = np.repeat(BSEDict['tflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.ifflag = np.repeat(BSEDict['ifflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.wdflag = np.repeat(BSEDict['wdflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.bhflag = np.repeat(BSEDict['bhflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.nsflag = np.repeat(BSEDict['nsflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.mxns = np.repeat(BSEDict['mxns'], self.initial_conditions.kstar1.size)
        self.initial_conditions.pts1 = np.repeat(BSEDict['pts1'], self.initial_conditions.kstar1.size)
        self.initial_conditions.pts2 = np.repeat(BSEDict['pts2'], self.initial_conditions.kstar1.size)
        self.initial_conditions.pts3 = np.repeat(BSEDict['pts3'], self.initial_conditions.kstar1.size)
        self.initial_conditions.sigma = np.repeat(BSEDict['sigma'], self.initial_conditions.kstar1.size)
        self.initial_conditions.beta = np.repeat(BSEDict['beta'], self.initial_conditions.kstar1.size)
        self.initial_conditions.xi = np.repeat(BSEDict['xi'], self.initial_conditions.kstar1.size)
        self.initial_conditions.acc2 = np.repeat(BSEDict['acc2'], self.initial_conditions.kstar1.size)
        self.initial_conditions.epsnov = np.repeat(BSEDict['epsnov'], self.initial_conditions.kstar1.size)
        self.initial_conditions.eddfac = np.repeat(BSEDict['eddfac'], self.initial_conditions.kstar1.size)
        self.initial_conditions.gamma = np.repeat(BSEDict['gamma'], self.initial_conditions.kstar1.size)
        self.initial_conditions.bconst = np.repeat(BSEDict['bconst'], self.initial_conditions.kstar1.size)
        self.initial_conditions.CK = np.repeat(BSEDict['CK'], self.initial_conditions.kstar1.size)
        self.initial_conditions.merger = np.repeat(BSEDict['merger'], self.initial_conditions.kstar1.size)
        self.initial_conditions.windflag = np.repeat(BSEDict['windflag'], self.initial_conditions.kstar1.size)
        self.initial_conditions.fbkickswitch = np.repeat(BSEDict['fbkickswitch'], self.initial_conditions.kstar1.size)

        initial_conditions = np.vstack([self.initial_conditions.kstar1, self.initial_conditions.kstar2, self.initial_conditions.mass1_binary, self.initial_conditions.mass2_binary, self.initial_conditions.porb, self.initial_conditions.ecc, self.initial_conditions.metallicity[0:self.initial_conditions.mass1_binary.size], self.initial_conditions.tphysf, self.initial_conditions.neta, self.initial_conditions.bwind, self.initial_conditions.hewind, self.initial_conditions.alpha1, self.initial_conditions.lambdaf, self.initial_conditions.ceflag, self.initial_conditions.tflag, self.initial_conditions.ifflag, self.initial_conditions.wdflag, self.initial_conditions.bhflag, self.initial_conditions.nsflag, self.initial_conditions.mxns, self.initial_conditions.pts1, self.initial_conditions.pts2, self.initial_conditions.pts3, self.initial_conditions.sigma, self.initial_conditions.beta, self.initial_conditions.xi, self.initial_conditions.acc2, self.initial_conditions.epsnov, self.initial_conditions.eddfac, self.initial_conditions.gamma, self.initial_conditions.bconst, self.initial_conditions.CK, self.initial_conditions.merger, self.initial_conditions.windflag, self.initial_conditions.fbkickswitch]).T

        # calculate maximum number of processes
        nproc = min(kwargs.pop('nproc', 1), len(initial_conditions))
        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                # kstar, mass, orbital period (days), eccentricity, metaliccity, evolution time (millions of years)
                print(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9],
                                        f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17], f[18], f[19],
                                        f[20], f[21], f[22], f[23], f[24], f[25], f[26], f[27], f[28], f[29],
                                        f[30], f[31], f[32], f[33], f[34])
                tmp = _evolvebin.evolv2(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9],
                                        f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17], f[18], f[19],
                                        f[20], f[21], f[22], f[23], f[24], f[25], f[26], f[27], f[28], f[29],
                                        f[30], f[31], f[32], f[33], f[34])
                return f, tmp[np.argwhere(tmp[:,0]>0),:]
            except Exception as e:
                if nproc == 1:
                    raise
                else:
                    return f, e

        # evolve sysyems
        output = mp_utils.multiprocess_with_queues(
            nproc, _evolve_single_system, initial_conditions, raise_exceptions=False)

        # raise exceptions (from multiprocessing, single process raises inline)
        for f, x in output:
            if isinstance(x, Exception):
                x.args = ('Failed to evolve %s: %s' % (f, str(x)),)
                raise x

        return output
