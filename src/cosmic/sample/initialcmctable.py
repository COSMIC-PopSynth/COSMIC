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

"""`InitialCMCTable`
"""

import pandas as pd
import numpy as np
import h5py
from astropy.io import fits
from astropy.table import Table

__author__ = "Scott Coughlin <scottcoughlin2014@u.northwestern.edu>"
__credits__ = "Carl Rodriguez <carllouisrodriguez@gmail.com>"
__all__ = ["InitialCMCTable"]


INITIAL_CONDITIONS_COLUMNS_CMC_SINGLES = [
    "id",
    "k",
    "m",
    "Reff",
    "r",
    "vr",
    "vt",
    "binind",
]

INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES = [
    "index",
    "id1",
    "k1",
    "m1",
    "Reff1",
    "id2",
    "k2",
    "m2",
    "Reff2",
    "a",
    "e",
]


class InitialCMCTable(pd.DataFrame):
    scaled_to_nbody_units = False
    metallicity = None
    mass_of_cluster = None
    virial_radius = None
    tidal_radius = None
    central_bh = 0.0
    scale_with_central_bh = False

    def ScaleCentralBHMass(self, Mtotal):
        """Rescale the central BH mass; needed since this is a class attribute
            Parameters
            ----------

            Mtotal : float
                total mass of the cluster
        """
        self.central_bh /= Mtotal

    @classmethod
    def ScaleToNBodyUnits(cls, Singles, Binaries, virial_radius=1, central_bh=0, scale_with_central_bh=False):
        """Rescale the single masses, radii, and velocities into N-body units
           i.e. \sum m = M = 1
                 Kinetic Energy   = 0.25
                 Potential Energy = -0.5

        Note that this is already done for r, vr, and vt from the profile generators.
        However, after the stellar masses are assigned we need to redo it, and 
        the stellar radii and seperations need to be converted from RSUN to code units

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

        # Normalize the masses to the total cluster mass
        M_total = sum(Singles["m"])
        Singles["m"] /= M_total
        Binaries["m1"] /= M_total
        Binaries["m2"] /= M_total

        # Note if there's an central BH, we want to add it
        # to the mass of the already normalized stars (i.e. M_tot = 1+BH)
        Singles.ScaleCentralBHMass(M_total)

        # Take the radii, and offset by one
        radius = np.array(Singles["r"])
        radius_p1 = np.append(radius[1:], [1e100])

        # Masses and velocities
        mass = np.array(Singles["m"])
        cumul_mass = np.cumsum(mass)
        vr = np.array(Singles["vr"])
        vt = np.array(Singles["vt"])

        # Then compute the total kinetic and potential energy
        # There's probably a cleaner way to do the PE (this is a one-line version
        #  of the for loop we use in CMC; vectorized and pythonic, but sloppy)
        KE = 0.5 * np.sum(mass * (vr ** 2 + vt ** 2))
        PE = 0.5 * np.sum(
            mass[::-1]
            * np.cumsum((cumul_mass * (1.0 / radius - 1.0 / radius_p1))[::-1])
        )

        # If scaling with central BH, add the potential from it
        if scale_with_central_bh:
            PE += np.sum(Singles.central_bh*mass / radius) 

        # Compute the position and velocity scalings
        rfac = 2 * PE
        vfac = 1.0 / np.sqrt(4 * KE)

        # Scale the positions and velocities s.t. KE=0.25, PE=-0.5
        Singles["r"] *= rfac
        Singles["vr"] *= vfac
        Singles["vt"] *= vfac

        # Finally, scale the radii and seperations from BSE into code units
        PARSEC_PER_RSUN = 2.2546101516664447e-08 
        DistConv = PARSEC_PER_RSUN / virial_radius

        Singles["Reff"] *= DistConv
        Binaries["a"] *= DistConv
        Binaries["Reff1"] *= DistConv
        Binaries["Reff2"] *= DistConv

        Singles.scaled_to_nbody_units = True
        Binaries.scaled_to_nbody_units = True
        return

    @classmethod
    def InitialCMCSingles(cls, id_idx, k, m, Reff, r, vr, vt, binind):
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
        bin_dat = cls(
            np.vstack(
                [
                    id_idx,
                    k,
                    m,
                    Reff,
                    r,
                    vr,
                    vt,
                    binind,
                ]
            ).T,
            columns=INITIAL_CONDITIONS_COLUMNS_CMC_SINGLES,
        )

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
        bin_dat = cls(
            np.vstack([index, id1, k1, m1, Reff1, id2, k2, m2, Reff2, a, e]).T,
            columns=INITIAL_CONDITIONS_COLUMNS_CMC_BINARIES,
        )

        return bin_dat

    @classmethod
    def write(cls, Singles, Binaries, filename="input.hdf5", **kwargs):
        """Save Singles and Binaries to HDF5 or FITS file

        Parameters
        ----------
        Singles : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        Binaries : DataFrame
            Pandas DataFrame from the InitialCMCBinaries function
        filename : (str)
            Must end in ".fits" or ".hdf5/h5"

        Optional Parameteres
        --------------------
        These are automatically set in the Singles df, but can be OVERWRITTEN here 

        virial_radius : `float`
            the initial virial radius of the cluster, in parsecs

        tidal_radius : `float`
            the initial tidal radius of the cluster, in units of the virial_radius

        metallicity : `float`
            the stellar metallicity of the cluster


        Returns
        -------
            None:

        """

        # verify parameters
        if (".hdf5" in filename) or (".h5" in filename):
            savehdf5 = True
            savefits = False
        elif ".fits" in filename:
            savefits = True
            savehdf5 = False
        else:
            raise ValueError(
                "File extension not recognized, valid file types are fits and hdf5"
            )

        virial_radius = kwargs.pop('virial_radius',Singles.virial_radius)
        tidal_radius = kwargs.pop('tidal_radius',Singles.tidal_radius)
        metallicity = kwargs.pop('metallicity',Singles.metallicity)
        central_bh = kwargs.pop('central_bh',Singles.central_bh)
        scale_with_central_bh = kwargs.pop('scale_with_central_bh',Singles.scale_with_central_bh)

        # If a user has not already scaled the units of the Singles and Binaries tables,
        # and the attribute mass_of_cluster is None, then
        # we can calculate it now
        if (not Singles.scaled_to_nbody_units) and (Singles.mass_of_cluster is None):
            Singles.mass_of_cluster = np.sum(Singles["m"]) + central_bh
            InitialCMCTable.ScaleToNBodyUnits(
                Singles, Binaries, virial_radius=virial_radius, central_bh=central_bh, scale_with_central_bh=scale_with_central_bh
            )
        elif (not Singles.scaled_to_nbody_units) and (Singles.mass_of_cluster is not None):
            InitialCMCTable.ScaleToNBodyUnits(
                Singles, Binaries, virial_radius=virial_radius, central_bh=central_bh, scale_with_central_bh=scale_with_central_bh
            )
        elif (Singles.scaled_to_nbody_units) and (Singles.mass_of_cluster is None):
            # we cannot get the pre-scaled mass of the cluster
            raise ValueError(
                "In order to save the initial conditions correctly, "
                "you cannot feed in Singles and Binaries which have already been scaled"
            )

        if Singles.metallicity is None:
            raise ValueError(
                "The user has not supplied a metallicity for the cluster. Please set the Singles.metallicity attribute"
            )

        # Need to append special rows to the start and end of the Singles table
        singles = pd.DataFrame(
            np.zeros((1, Singles.shape[1])), index=[0], columns=Singles.columns
        )
        singles_bottom = pd.DataFrame(
            np.zeros((1, Singles.shape[1])), index=[0], columns=Singles.columns
        )
        singles = pd.concat([singles, Singles])
        singles = pd.concat([singles, singles_bottom])
        
        singles.iloc[-1, singles.columns.get_loc("r")] = 1e40
        singles.iloc[0, singles.columns.get_loc("r")] = 2.2250738585072014e-308
        singles.iloc[0, singles.columns.get_loc("m")] = Singles.central_bh

        # Add a special row to the end of Bianries table
        binaries = pd.DataFrame(
            np.zeros((1, Binaries.shape[1])), index=[0], columns=Binaries.columns
        )
        binaries = pd.concat([binaries, Binaries])

        if savehdf5:
            singles.to_hdf(filename, key="CLUS_OBJ_DATA", mode="w")
            binaries.to_hdf(filename, key="CLUS_BINARY_DATA")
            with h5py.File(filename, "a") as f:
                f["CLUS_OBJ_DATA/block0_values"].attrs["EXTNAME"] = "CLUS_OBJ_DATA"
                f["CLUS_OBJ_DATA/block0_values"].attrs["NOBJ"] = int(len(singles)) - 2
                f["CLUS_OBJ_DATA/block0_values"].attrs["NBINARY"] = (
                    int(len(binaries)) - 1
                )
                f["CLUS_OBJ_DATA/block0_values"].attrs[
                    "MCLUS"
                ] = Singles.mass_of_cluster
                f["CLUS_OBJ_DATA/block0_values"].attrs["RVIR"] = virial_radius 
                f["CLUS_OBJ_DATA/block0_values"].attrs["RTID"] = tidal_radius 
                f["CLUS_OBJ_DATA/block0_values"].attrs["Z"] = metallicity

        if savefits:
            Singles_fits = Table.from_pandas(singles)
            Binaries_fits = Table.from_pandas(binaries)

            hdu1 = fits.table_to_hdu(Singles_fits)
            hdu2 = fits.table_to_hdu(Binaries_fits)

            # create a header
            hdr = fits.Header()
            hdr["COMMENT"] = "CMC Configured Initial Conditions"
            hdr["COMMENT"] = "Produced by COSMIC"
            primary_hdu = fits.PrimaryHDU(header=hdr)

            hdu1.header["EXTNAME"] = "CLUS_OBJ_DATA"
            hdu1.header["NOBJ"] = int(len(Singles_fits)) - 2
            hdu1.header["NBINARY"] = int(len(Binaries_fits)) - 1
            hdu1.header["MCLUS"] = Singles.mass_of_cluster
            hdu1.header["RVIR"] = virial_radius
            hdu1.header["RTID"] = tidal_radius
            hdu1.header["Z"] = metallicity

            # put all the HDUs together
            hdul = fits.HDUList([primary_hdu, hdu1, hdu2])

            # write it out
            hdul.writeto(
                filename,
                overwrite=True,
            )

        return

    @classmethod
    def read(cls, filename):
        """Read Singles and Binaries to HDF5 or FITS file

        Parameters
        ----------
        filename : (str)
            Must end in ".fits" or ".hdf5/h5"

        Returns
        -------
        Singles : DataFrame
            Pandas DataFrame from the InitialCMCSingles function
        Binaries : DataFrame
            Pandas DataFrame from the InitialCMCBinaries function
        """
        # verify parameters
        if (".hdf5" in filename) or (".h5" in filename):
            savehdf5 = True
            savefits = False
        elif ".fits" in filename:
            savefits = True
            savehdf5 = False
        else:
            raise ValueError(
                "File extension not recognized, valid file types are fits and hdf5"
            )

        if savehdf5:
            Singles = pd.read_hdf(filename, "CLUS_OBJ_DATA")
            Binaries = pd.read_hdf(filename, "CLUS_BINARY_DATA")
        elif savefits:
            Singles = cls(Table.read(filename,hdu=1).to_pandas())
            Binaries = cls(Table.read(filename, hdu=2).to_pandas())
        return Singles, Binaries


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
