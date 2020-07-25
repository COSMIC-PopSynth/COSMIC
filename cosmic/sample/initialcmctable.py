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

"""`InitialCMCTable`
"""

import pandas as pd
import numpy as np

__author__ = 'Scott Coughlin <scottcoughlin2014@u.northwestern.edu>'
__credits__ = 'Carl Rodriguez <carllouisrodriguez@gmail.com>'
__all__ = ['InitialCMCTable']


INITIAL_CONDITIONS_COLUMNS_CMC_SINGLES = ['id', 'k', 'm', 'Reff', 'r', 'vr', 'vt', 'binind']

INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES = ['index', 'id1', 'k1', 'm1', 'Reff1', 'id2', 'k2', 'm2', 'Reff2', 'a', 'e']

class InitialCMCTable():
    @classmethod
    def InitialCMCSingles(cls, id_idx, k, m, Reff, r, vr,vt, binind):
        """Create A Table of CMC Singles

        Parameters
        ----------
        m1 : float
            Primary mass [Msun]
        m2 : float
            Secondary mass [Msun]
        porb : float
            Orbital period [days]
        ecc : float
            Eccentricity
        kstar1 : array
            0-14 Initial stellar type of the larger object;
            main sequence stars are 0 if m < 0.7 Msun and 1 otherwise
        kstar2 : array
            0-14 Initial stellar type of the smaller object;
            main sequence stars are 0 if m < 0.7 Msun and 1 otherwise
        metallicity : float
            Metallicity of the binaries; Z_sun = 0.02

        **kwargs

            binfrac : float
                System-specific probability of the primary star being in a binary

        Returns
        -------
        InitialBinaries : DataFrame
            Single binary initial conditions

        """
        bin_dat = pd.DataFrame(np.vstack([
                                            id_idx, k, m, Reff, r, vr,vt, binind,
                                          ]).T,
                               columns = INITIAL_CONDITIONS_COLUMNS_CMC_SINGLES)

        return bin_dat

    @classmethod
    def InitialCMCBinaries(cls, index, id1, k1, m1, Reff1, id2, k2, m2, Reff2, a, e):
        """Create A Table of CMC Binaries

        Parameters
        ----------
        m1 : float
            Primary mass [Msun]
        m2 : float
            Secondary mass [Msun]
        porb : float
            Orbital period [days]
        ecc : float
            Eccentricity
        kstar1 : array
            0-14 Initial stellar type of the larger object;
            main sequence stars are 0 if m < 0.7 Msun and 1 otherwise
        kstar2 : array
            0-14 Initial stellar type of the smaller object;
            main sequence stars are 0 if m < 0.7 Msun and 1 otherwise
        metallicity : float
            Metallicity of the binaries; Z_sun = 0.02

        **kwargs

            binfrac : float
                System-specific probability of the primary star being in a binary

        Returns
        -------
        InitialBinaries : DataFrame
            Single binary initial conditions

        """
        bin_dat = pd.DataFrame(np.vstack([
                                          index, id1, k1, m1, Reff1, id2, k2, m2, Reff2, a, e
                                          ]).T,
                               columns = INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES)

        return bin_dat

    @classmethod
    def sampler(cls, format_, *args, **kwargs):
        """Fetch a method to generate an initial binary sample

        Parameters
        ----------
        format : str
            the method name; Choose from 'independent' or 'multidim'

        *args
            the arguments necessary for the registered sample
            method; see help(InitialCMCTable.sampler('independent')
            to see the arguments necessary for the independent sample
        """
        # standard registered fetch
        from .sampler.sampler import get_sampler
        sampler = get_sampler(format_, cls)
        return sampler(*args, **kwargs)
