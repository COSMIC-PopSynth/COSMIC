# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of astro-traj.
#
# astro-traj is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# astro-traj is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with astro-traj.  If not, see <http://www.gnu.org/licenses/>.

"""`sample`
"""

import numpy as np
from gwpy.utils import mp as mp_utils
from aCOSMIC import _popbintd

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'Evolve'

class Evolve:
    def __init__(self):
        '''
        initialize Evolve
        '''


    def evolve(self, **kwargs):
        """After setting a number of initial conditions
        we evolve the system.

        Parameters
        ----------
        nproc : `int`, optional, default: 1
            number of CPUs to use for parallel file reading

        kwargs: Optional TrainingSet and LabelledSamples args
            that will download images in a specila way
            ./"Label"/"SampleType"/"image"

        Returns
        -------
        An evolve binary
        """

        # define multiprocessing method
        def _evolve_single_system(f):
            try:
                return f, _popbintd.evolv2(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14], f[15], f[16], f[17], f[18], f[19], f[20], f[21], f[22], f[23])
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
