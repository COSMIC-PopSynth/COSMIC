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
# along with astro-traj.  If not, see <http://www.gnu.org/licenses/>.

"""`sample`
"""

import numpy as np
import multiprocessing as mp
import math
import random
import scipy.integrate

from astropy.table import Table, Column
from astropy import units

from aCOSMIC.utils import idl_tabulate, rndm

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'Sample'


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


class InitialBinaryTable(Table):
    @classmethod
    def SingleBinary(cls, m1, m2, porb, ecc, tphysf, kstar1, kstar2, metallicity):
        """Create single binary

            Parameters:

                m1 (float) :
                    This is the mass in solar masses of the larger object

                m2 (float) :
                    This is the mass in solar masses of the smaller object

                porb (float) :
                    This is the orbital period in days

                ecc (float) :
                    Eccentricity This value is between 0 and 1

                tphysf (float) :
                    How long to evolve the binary in millions of years

                kstar1 (array) :
                    0-14 Initial stellar type of the larger object; 
                    main sequence stars are 0 if m < 0.7 Msun and 1 otherwise

                kstar2 (array) :
                    0-14 Initial stellar type of the smaller object; 
                    main sequence stars are 0 if m < 0.7 Msun and 1 otherwise

                metallicity (float):
                    metallicity of the galaxy where the binary lives.


            Returns:
                `SingleBinary`
        """
        tab = Table()
        tab['mass1_binary'] = Column([m1], unit=units.Msun,
                               description='Mass of larger object')
        tab['mass2_binary'] = Column([m2], unit=units.Msun,
                               description='Mass of smaller object') 
        tab['porb'] = Column([porb], unit=units.day,
                               description='Orbital Period')
        tab['ecc'] = Column([ecc], unit=None,
                               description='Eccentricity')
        tab['tphysf'] = Column([tphysf], unit=units.year*1000000,
                               description='Binary Evolution Time')
        tab['kstar1'] = Column([kstar1], unit=None,
                               description='Stellar type of larger object')
        tab['kstar2'] = Column([kstar2], unit=None,
                               description='Stelar type of smaller object')
        tab['metallicity'] = Column([metallicity], unit=None,
                               description='metallicity of the galaxy')

        return tab

    @classmethod
    def BinaryGrid(cls, m1, m2, porb, ecc, tphysf, kstar1, kstar2, metallicity):
        """Create a grid of binaries

            Parameters:

                m1 (array) :
                    This is the mass in solar masses of the larger object

                m2 (array) :
                    This is the mass in solar masses of the smaller object

                porb (array) :
                    This is the orbital period in days

                ecc (array) :
                    Eccentricity This value is between 0 and 1

                tphysf (array) :
                    How long to evolve the binary in millions of years

                kstar1 (array) :
                    0-14 Initial stellar type of the larger object; 
                    main sequence stars are 0 if m < 0.7 Msun and 1 otherwise

                kstar2 (array) :
                    0-14 Initial stellar type of the smaller object; 
                    main sequence stars are 0 if m < 0.7 Msun and 1 otherwise

                metallicity (array):
                    metallicity of the galaxy where the binary lives.


            Returns:
                `BinaryGrid`
        """
        tab = Table()
        tab['mass1_binary'] = Column(m1, unit=units.Msun,
                               description='Mass of larger object')
        tab['mass2_binary'] = Column(m2, unit=units.Msun,
                               description='Mass of smaller object')
        tab['porb'] = Column(porb, unit=units.day,
                               description='Orbital Period')
        tab['ecc'] = Column(ecc, unit=None,
                               description='Eccentricity')
        tab['tphysf'] = Column(tphysf, unit=units.year*1000000,
                               description='Binary Evolution Time')
        tab['kstar1'] = Column(kstar1, unit=None,
                               description='Stellar type of larger object')
        tab['kstar2'] = Column(kstar2, unit=None,
                               description='Stelar type of smaller object')
        tab['metallicity'] = Column(metallicity, unit=None,
                               description='metallicity of the galaxy')

        return tab

    @classmethod
    def sampler(cls, format_, *args, **kwargs):
        """Fetch a table of events from a database

        Parameters
        ----------
        format : `str`, `~sqlalchemy.engine.Engine`
            the format of the remote data, see _Notes_ for a list of
            registered formats, OR an SQL database `Engine` object

        *args
            all other positional arguments are specific to the
            data format, see below for basic usage
        Notes
        -----"""
        # standard registered fetch
        from .sampler.sampler import get_sampler
        sampler = get_sampler(format_, cls)
        return sampler(*args, **kwargs)
