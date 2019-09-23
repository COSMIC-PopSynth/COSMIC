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

"""`InitialBinaryTable`
"""

import numpy as np
import multiprocessing as mp
import math
import random
import scipy.integrate
import sys

from astropy.table import Table, Column
from astropy import units
import pandas as pd

from cosmic.utils import idl_tabulate, rndm

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['InitialBinaryTable']


G = 6.67384*math.pow(10, -11.0)
c = 2.99792458*math.pow(10, 8.0)
parsec = 3.08567758*math.pow(10, 16)
Rsun = 6.955*math.pow(10, 8)
Msun = 1.9891*math.pow(10,30)
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569*10**7.0
Tobs = 3.15569*10**7.0
geo_mass = G/c**2

INITIAL_CONDITIONS_COLUMNS = ['kstar_1', 'kstar_2', 'mass1_binary', 'mass2_binary', 'porb', 'ecc',
                             'metallicity', 'tphysf', 'mass0_1', 'mass0_2',
                             'rad1', 'rad2', 'lumin1', 'lumin2', 'massc1', 'massc2',
                             'radc1', 'radc2', 'menv1', 'menv2', 'renv1', 'renv2',
                             'ospin1', 'ospin2', 'b_0_1', 'b_0_2', 'bacc1', 'bacc2',
                             'tacc1', 'tacc2', 'epoch1', 'epoch2', 'tms1', 'tms2',
                             'bhspin1','bhspin2', 'tphys']

INITIAL_CONDITIONS_MISC_COLUMNS = ['binfrac']

if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_COLUMNS_ALL = INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_COLUMNS_ALL = INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_COLUMNS_ALL.extend(INITIAL_CONDITIONS_MISC_COLUMNS)


class InitialBinaryTable():
    @classmethod
    def InitialBinaries(cls, m1, m2, porb, ecc, tphysf, kstar1, kstar2, metallicity, **kwargs):
        """Create single binary

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
        tphysf : float
            Time to evolve the binary [Myr]
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

            mass0_1,mass0_2,rad1,rad2,lumin1,lumin2,
            massc1,massc2,radc1,radc2,menv1,menv2,renv1,renv2,
            ospin1,ospin2,b_0_1,b_0_2,bacc1,bacc2,
            tacc1,tacc2,epoch1,epoch2,tms1,tms2
            bhspin1,bhspin2

        Returns
        -------
        InitialBinaries : DataFrame
            Single binary initial conditions

        """
        # System-specific probability of the primary star being in a binary
        binfrac = kwargs.pop('binfrac', np.ones(np.array(m1).size))

        # The following are in general not needed for COSMIC, however
        # if you wish to start a binary at a very specific point in its evolution
        # these kwarg arguments can help you do this.
        # For instance the Globular Cluster code CMC requires this behavior.
        mass0_1 = kwargs.pop('mass0_1', m1)
        mass0_2 = kwargs.pop('mass0_2', m2)
        rad1 = kwargs.pop('rad1', np.zeros(np.array(m1).size))
        rad2 = kwargs.pop('rad2', np.zeros(np.array(m1).size))
        lumin1 = kwargs.pop('lumin1', np.zeros(np.array(m1).size))
        lumin2 = kwargs.pop('lumin2', np.zeros(np.array(m1).size))
        massc1 = kwargs.pop('massc1', np.zeros(np.array(m1).size))
        massc2 = kwargs.pop('massc2', np.zeros(np.array(m1).size))
        radc1 = kwargs.pop('radc1', np.zeros(np.array(m1).size))
        radc2 = kwargs.pop('radc2', np.zeros(np.array(m1).size))
        menv1 = kwargs.pop('menv1', np.zeros(np.array(m1).size))
        menv2 = kwargs.pop('menv2', np.zeros(np.array(m1).size))
        renv1 = kwargs.pop('renv1', np.zeros(np.array(m1).size))
        renv2 = kwargs.pop('renv2', np.zeros(np.array(m1).size))
        ospin1 = kwargs.pop('ospin1', np.zeros(np.array(m1).size))
        ospin2 = kwargs.pop('ospin2', np.zeros(np.array(m1).size))
        b_0_1 = kwargs.pop('b_0_1', np.zeros(np.array(m1).size))
        b_0_2 = kwargs.pop('b_0_2', np.zeros(np.array(m1).size))
        bacc1 = kwargs.pop('bacc1', np.zeros(np.array(m1).size))
        bacc2 = kwargs.pop('bacc2', np.zeros(np.array(m1).size))
        tacc1 = kwargs.pop('tacc1', np.zeros(np.array(m1).size))
        tacc2 = kwargs.pop('tacc2', np.zeros(np.array(m1).size))
        epoch1 = kwargs.pop('epoch1', np.zeros(np.array(m1).size))
        epoch2 = kwargs.pop('epoch2', np.zeros(np.array(m1).size))
        tms1 = kwargs.pop('tms1', np.zeros(np.array(m1).size))
        tms2 = kwargs.pop('tms2', np.zeros(np.array(m1).size))
        bhspin1 = kwargs.pop('bhspin1', np.zeros(np.array(m1).size))
        bhspin2 = kwargs.pop('bhspin2', np.zeros(np.array(m1).size))
        tphys = kwargs.pop('tphys', np.zeros(np.array(m1).size))

        bin_dat = pd.DataFrame(np.vstack([
                                            kstar1, kstar2,
                                            m1, m2, porb, ecc,
                                            metallicity, tphysf,
                                            mass0_1, mass0_2, rad1, rad2,
                                            lumin1, lumin2, massc1, massc2,
                                            radc1, radc2, menv1, menv2, renv1, renv2,
                                            ospin1, ospin2, b_0_1, b_0_2, bacc1, bacc2,
                                            tacc1, tacc2, epoch1, epoch2, tms1, tms2,
                                            bhspin1,bhspin2,tphys,binfrac,
                                          ]).T,
                               columns = INITIAL_CONDITIONS_COLUMNS_ALL)

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
            method; see help(InitialBinaryTable.sampler('independent')
            to see the arguments necessary for the independent sample
        """
        # standard registered fetch
        from .sampler.sampler import get_sampler
        sampler = get_sampler(format_, cls)
        return sampler(*args, **kwargs)
