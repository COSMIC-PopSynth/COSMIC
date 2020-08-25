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
import h5py
from astropy.io import fits
from astropy.table import Table

__author__ = 'Scott Coughlin <scottcoughlin2014@u.northwestern.edu>'
__credits__ = 'Carl Rodriguez <carllouisrodriguez@gmail.com>'
__all__ = ['InitialCMCTable']


INITIAL_CONDITIONS_COLUMNS_CMC_SINGLES = ['id', 'k', 'm', 'Reff', 'r', 'vr', 'vt', 'binind']

INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES = ['index', 'id1', 'k1', 'm1', 'Reff1', 'id2', 'k2', 'm2', 'Reff2', 'a', 'e']

class InitialCMCTable(pd.DataFrame):
    scaled_to_nbody_units = False
    metallicity = None
    mass_of_cluster = None

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

        Singles.scaled_to_nbody_units = True
        Binaries.scaled_to_nbody_units = True
        return

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
        bin_dat = cls(np.vstack([
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
        bin_dat = cls(np.vstack([
                                index, id1, k1, m1, Reff1, id2, k2, m2, Reff2, a, e
                                ]).T,
                                columns = INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES)

        return bin_dat

    @classmethod
    def write(cls, Singles, Binaries, filename="input.hdf5", **kwargs):
        """Save Singles and Binaries to HDF5 or FITS file

        Parameters
        ----------
        Singles : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        Binaries : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        filename : (str)
            Must end in ".fits" or ".hdf5/h5"

        **kwargs

            virial_radius
            rtid
        Returns
        -------
            None:

        """
        virial_radius = kwargs.pop("virial_radius", 1)
        rtid = kwargs.pop("rtid", 1000000.0)

        # verify parameters
        if ( ('.hdf5' in filename) or ('.h5' in filename) ):
            savehdf5 = True
            savefits = False
        elif '.fits' in filename:
            savefits = True
            savehdf5 = False
        else:
            raise ValueError("File extension not recognized, valid file types are fits and hdf5")

        # If a user has not already scaled the units of the Singles and Binaries tables, and the attribute mass_of_cluster is None, then
        # we can calculate it now
        if (not Singles.scaled_to_nbody_units) and (Singles.mass_of_cluster is None):
            Singles.mass_of_cluster = np.sum(Singles['m'])
            InitialCMCTable.ScaleToNBodyUnits(Singles, Binaries, virial_radius=virial_radius)
        elif (Singles.scaled_to_nbody_units) and (Singles.mass_of_cluster is None):
            # we cannot get the pre-scaled mass of the cluster
            raise ValueError("In order to save the initial conditions correctly, "
                             "you cannot feed in Singles and Binaries which have already been scaled")

        if Singles.metallicity is None:
            raise ValueError("The user has not supplied a metallicity for the cluster. Please set the Singles.metallicity attribute")

        # Need to append special rows to the start and end of the Singles table
        singles = pd.DataFrame(np.zeros((1, Singles.shape[1])), index=[0], columns=Singles.columns)
        singles_bottom = pd.DataFrame(np.zeros((1, Singles.shape[1])), index=[0], columns=Singles.columns)
        singles = singles.append(Singles)
        singles = singles.append(singles_bottom)
        singles['r'].iloc[-1] = 1e40
        singles['r'].iloc[0] = 2.2250738585072014e-308

        # Add a special row to the end of Bianries table
        binaries = pd.DataFrame(np.zeros((1, Binaries.shape[1])), index=[0], columns=Binaries.columns)
        binaries = binaries.append(Binaries)

        if savehdf5:
            singles.to_hdf(filename, key="CLUS_OBJ_DATA", mode='w')
            binaries.to_hdf(filename, key="CLUS_BINARY_DATA")
            with h5py.File(filename, 'a') as f:
                f["CLUS_OBJ_DATA/block0_values"].attrs['EXTNAME'] = 'CLUS_OBJ_DATA'
                f["CLUS_OBJ_DATA/block0_values"].attrs['NOBJ'] = int(len(singles)) - 2
                f["CLUS_OBJ_DATA/block0_values"].attrs['NBINARY'] =  int(len(binaries)) - 1
                f["CLUS_OBJ_DATA/block0_values"].attrs['MCLUS'] = Singles.mass_of_cluster
                f["CLUS_OBJ_DATA/block0_values"].attrs['RVIR'] = virial_radius
                f["CLUS_OBJ_DATA/block0_values"].attrs['RTID'] = rtid
                f["CLUS_OBJ_DATA/block0_values"].attrs['Z'] = Singles.metallicity


        if savefits:
            Singles_fits = Table.from_pandas(singles)
            Binaries_fits = Table.from_pandas(binaries)

            hdu1 = fits.table_to_hdu(Singles_fits)
            hdu2 = fits.table_to_hdu(Binaries_fits)

            # create a header
            hdr = fits.Header()
            hdr['COMMENT'] = 'CMC Configured Initial Conditions'
            hdr['COMMENT'] = 'Produced by COSMIC'
            primary_hdu = fits.PrimaryHDU(header=hdr)

            hdu1.header['EXTNAME'] = 'CLUS_OBJ_DATA'
            hdu1.header['NOBJ'] = int(len(Singles_fits)) - 2
            hdu1.header['NBINARY'] =  int(len(Binaries_fits)) - 1
            hdu1.header['MCLUS'] = Singles.mass_of_cluster
            hdu1.header['RVIR'] = virial_radius
            hdu1.header['RTID'] = rtid
            hdu1.header['Z'] = Singles.metallicity


            # put all the HDUs together
            hdul = fits.HDUList([primary_hdu, hdu1, hdu2])

            # write it out
            hdul.writeto(filename, overwrite=True,)

        return


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
