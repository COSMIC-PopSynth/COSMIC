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

"""`InitialBinaryTable`
"""

import numpy as np
import math
import sys

import pandas as pd

__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = ["InitialBinaryTable"]


G = 6.67384 * math.pow(10, -11.0)
c = 2.99792458 * math.pow(10, 8.0)
parsec = 3.08567758 * math.pow(10, 16)
Rsun = 6.955 * math.pow(10, 8)
Msun = 1.9891 * math.pow(10, 30)
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569 * 10 ** 7.0
Tobs = 3.15569 * 10 ** 7.0
geo_mass = G / c ** 2

INITIAL_CONDITIONS_COLUMNS = []

INITIAL_CONDITIONS_COLUMNS_CORE = [
    "kstar_1",
    "kstar_2",
    "mass_1",
    "mass_2",
    "porb",
    "ecc",
    "metallicity",
    "tphysf",
]

INITIAL_CONDITIONS_COLUMNS_EXTRAS = [
    "mass0_1",
    "mass0_2",
    "rad_1",
    "rad_2",
    "lum_1",
    "lum_2",
    "massc_1",
    "massc_2",
    "radc_1",
    "radc_2",
    "menv_1",
    "menv_2",
    "renv_1",
    "renv_2",
    "omega_spin_1",
    "omega_spin_2",
    "B_1",
    "B_2",
    "bacc_1",
    "bacc_2",
    "tacc_1",
    "tacc_2",
    "epoch_1",
    "epoch_2",
    "tms_1",
    "tms_2",
    "bhspin_1",
    "bhspin_2",
    "tphys",
]

INITIAL_CONDITIONS_COLUMNS.extend(INITIAL_CONDITIONS_COLUMNS_CORE)
INITIAL_CONDITIONS_COLUMNS.extend(INITIAL_CONDITIONS_COLUMNS_EXTRAS)

INITIAL_CONDITIONS_MISC_COLUMNS = ["binfrac"]

if sys.version_info.major == 2 and sys.version_info.minor == 7:
    INITIAL_CONDITIONS_COLUMNS_ALL = INITIAL_CONDITIONS_COLUMNS[:]
else:
    INITIAL_CONDITIONS_COLUMNS_ALL = INITIAL_CONDITIONS_COLUMNS.copy()

INITIAL_CONDITIONS_COLUMNS_ALL.extend(INITIAL_CONDITIONS_MISC_COLUMNS)


class InitialBinaryTable:
    @classmethod
    def InitialBinaries(
        cls, m1, m2, porb, ecc, tphysf, kstar1, kstar2, metallicity, **kwargs
    ):
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
        binfrac = kwargs.pop("binfrac", np.ones(np.array(m1).size))

        # The following are in general not needed for COSMIC, however
        # if you wish to start a binary at a very specific point in its evolution
        # these kwarg arguments can help you do this.
        # For instance the Globular Cluster code CMC requires this behavior.
        mass0_1 = kwargs.pop("mass0_1", m1)
        mass0_2 = kwargs.pop("mass0_2", m2)
        rad1 = kwargs.pop("rad_1", np.zeros(np.array(m1).size))
        rad2 = kwargs.pop("rad_2", np.zeros(np.array(m1).size))
        lumin1 = kwargs.pop("lumin_1", np.zeros(np.array(m1).size))
        lumin2 = kwargs.pop("lumin_2", np.zeros(np.array(m1).size))
        massc1 = kwargs.pop("massc_1", np.zeros(np.array(m1).size))
        massc2 = kwargs.pop("massc_2", np.zeros(np.array(m1).size))
        radc1 = kwargs.pop("radc_1", np.zeros(np.array(m1).size))
        radc2 = kwargs.pop("radc_2", np.zeros(np.array(m1).size))
        menv1 = kwargs.pop("menv_1", np.zeros(np.array(m1).size))
        menv2 = kwargs.pop("menv_2", np.zeros(np.array(m1).size))
        renv1 = kwargs.pop("renv_1", np.zeros(np.array(m1).size))
        renv2 = kwargs.pop("renv_2", np.zeros(np.array(m1).size))
        ospin1 = kwargs.pop("ospin_1", np.zeros(np.array(m1).size))
        ospin2 = kwargs.pop("ospin_2", np.zeros(np.array(m1).size))
        b_0_1 = kwargs.pop("B_1", np.zeros(np.array(m1).size))
        b_0_2 = kwargs.pop("B_2", np.zeros(np.array(m1).size))
        bacc1 = kwargs.pop("bacc_1", np.zeros(np.array(m1).size))
        bacc2 = kwargs.pop("bacc_2", np.zeros(np.array(m1).size))
        tacc1 = kwargs.pop("tacc_1", np.zeros(np.array(m1).size))
        tacc2 = kwargs.pop("tacc_2", np.zeros(np.array(m1).size))
        epoch1 = kwargs.pop("epoch_1", np.zeros(np.array(m1).size))
        epoch2 = kwargs.pop("epoch_2", np.zeros(np.array(m1).size))
        tms1 = kwargs.pop("tms_1", np.zeros(np.array(m1).size))
        tms2 = kwargs.pop("tms_2", np.zeros(np.array(m1).size))
        bhspin1 = kwargs.pop("bhspin_1", np.zeros(np.array(m1).size))
        bhspin2 = kwargs.pop("bhspin_2", np.zeros(np.array(m1).size))
        tphys = kwargs.pop("tphys", np.zeros(np.array(m1).size))

        bin_dat = pd.DataFrame(
            np.vstack(
                [
                    kstar1,
                    kstar2,
                    m1,
                    m2,
                    porb,
                    ecc,
                    metallicity,
                    tphysf,
                    mass0_1,
                    mass0_2,
                    rad1,
                    rad2,
                    lumin1,
                    lumin2,
                    massc1,
                    massc2,
                    radc1,
                    radc2,
                    menv1,
                    menv2,
                    renv1,
                    renv2,
                    ospin1,
                    ospin2,
                    b_0_1,
                    b_0_2,
                    bacc1,
                    bacc2,
                    tacc1,
                    tacc2,
                    epoch1,
                    epoch2,
                    tms1,
                    tms2,
                    bhspin1,
                    bhspin2,
                    tphys,
                    binfrac,
                ]
            ).T,
            columns=INITIAL_CONDITIONS_COLUMNS_ALL,
        )

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
