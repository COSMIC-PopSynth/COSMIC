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
    def ScaleToNBodyUnits(cls, Singles, Binaries, virial_radius=1):
        """Rescale the single masses, radii, and velocities into N-body units
           i.e. \sum m = M = 1
                 Kinetic Energy   = 0.25
                 Potential Energy = -0.5

        Parameters
        ----------
        Singles : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        Binaries : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        virial_radius : float
            Virial radius of the cluster in parsec (default 1pc)

        Returns
        -------
        None: Pandas dataframes are modified in place

        """

        ## Normalize the masses to the total cluster mass
        M_total = sum(Singles['m'])
        Singles['m'] /= M_total
        Binaries['m1'] /= M_total
        Binaries['m2'] /= M_total

        ## Take the radii, and offset by one
        radius = np.array(Singles['r'])
        radius_p1 = np.append(radius[1:],[1e100])

        ## Masses and velocities
        mass = np.array(Singles['m'])
        cumul_mass = np.cumsum(mass)
        vr = np.array(Singles['vr'])
        vt = np.array(Singles['vt'])

        ## Then compute the total kinetic and potential energy
        ## There's probably a cleaner way to do the PE (this is a one-line version 
        ##  of the for loop we use in CMC; vectorized and pythonic, but sloppy)
        KE = 0.5*sum(mass*(vr**2 + vt**2))
        PE = 0.5*sum(mass[::-1]*np.cumsum((cumul_mass*(1.0/radius - 1.0/radius_p1))[::-1]))

        ## Compute the position and velocity scalings
        rfac = 2*PE
        vfac = 1.0/np.sqrt(4*KE)

        ## Scale the positions and velocities s.t. KE=0.25, PE=-0.5
        Singles['r'] *= rfac
        Singles['vr'] *= vfac
        Singles['vt'] *= vfac

        ## Finally, scale the radii and seperations from BSE into code units  
        PARSEC_TO_AU = 206264.806
        DistConv = 1/PARSEC_TO_AU/virial_radius

        Singles['Reff'] *= DistConv
        Binaries['a'] *= DistConv
        Binaries['Reff1'] *= DistConv
        Binaries['Reff2'] *= DistConv



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
