# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2021)
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

"""`MC_samp`
"""

import numpy as np
import math
import scipy.special as ss

__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = [
    "mass_weighted_number",
    "select_component_mass",
    "sample_sech_squared",
    "sample_exponential_radial",
    "sample_exponential_vertical",
    "sample_exponential_square_radial",
    "galactic_positions",
]


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

# solar coordinates in the galaxy: in parsecs from
# (Chaper 8 of Galactic Structure and stellar Pops book) Yoshii (2013)
############################################################################
x_sun = 8.0
y_sun = 0.0
z_sun = 25.0 / 1000.0


def mass_weighted_number(dat, total_sampled_mass, component_mass):
    """Compute the total number of systems in the synthetic catalog
    based on the total sampled mass of the simulated system and
    the total mass of a given galactic component

    Parameters
    ----------
    dat : DataFrame
        DataFrame containing the fixed population created from cosmic-pop
    total_sampled_mass : float
        total amount of mass sampled to generate the fixed population
        including single stars
    component_mass : float
        mass of the Galactic component we are simulating

    Returns
    -------
    nSystems : int
        number of systems in a Milky Way population for the selected
        Galactic population and fixed population
    """

    nSystems = int(len(dat.mass_1) * component_mass / total_sampled_mass)
    return nSystems


def select_component_mass(gx_component):
    """Select the Galactic component mass according to
    McMillan (2011)

    Parameters
    ----------
    gx_component : str
        Choose from: 'ThinDisk', 'Bulge', 'ThickDisk'

    Returns
    -------
    gx_component_mass : float
        Galactic component mass [Msun]
    """

    if gx_component == "ThinDisk":
        gx_component_mass = 4.32e10
    elif gx_component == "Bulge":
        gx_component_mass = 8.9e9
    elif gx_component == "ThickDisk":
        gx_component_mass = 1.44e10

    return gx_component_mass


def sample_sech_squared(size, scale_height=0.3):
    """Sample a collection of numbers of size=size distributed according
    to a sech sqaured function with a user specific scale height

    Parameters
    ----------
    size : int
        Size of the sample
    scale_height : float
        Scale height of the distribution; Default=0.3

    Returns
    -------
    distributed_nums : array
        Array of sampled values
    """

    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = np.arctanh((2 * rand_nums - 1)) * scale_height

    return distributed_nums


def sample_exponential_radial(size, scale_height):
    """Sample a collection of numbers of size=size distributed according
    to a radial exponential function with a user specific scale height

    Parameters
    ----------
    size : int
        Size of the sample
    scale_height : float
        Scale height of the distribution

    Returns
    -------
    distributed_nums : array
        Array of sampled values
    """

    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = -scale_height * np.log(1.0 - rand_nums)

    return distributed_nums


def sample_exponential_vertical(size, scale_height):
    """Sample a collection of numbers of size=size distributed according
    to a vertical exponential function with a user specific scale height

    Parameters
    ----------
    size : int
        Size of the sample
    scale_height : float
        Scale height of the distribution

    Returns
    -------
    distributed_nums : array
        Array of sampled values
    """

    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = -scale_height * np.log(1.0 - rand_nums)

    # CHOOSE A POSITION ABOVE AND BELOW THE DISK
    ###########################################################################
    pos_neg_choose = np.random.uniform(0, 1, size)
    (negInd,) = np.where(pos_neg_choose < 0.5)
    distributed_nums[negInd] = distributed_nums[negInd] * (-1.0)

    return distributed_nums


def sample_exponential_square_radial(size, scale_height):
    """Sample a collection of numbers of size=size distributed according
    to an exponential squared function with a user specific scale height

    Parameters
    ----------
    size : int
        Size of the sample
    scale_height : float
        Scale height of the distribution

    Returns
    -------
    distributed_nums : array
        Array of sampled values
    """

    rand_nums = np.random.uniform(0, 1, size)
    distributed_nums = scale_height * ss.erfinv(rand_nums)

    return distributed_nums


def galactic_positions(gx_component, size, model="McMillan"):
    """Sample a set of Galactic positions of size=size distributed
    according to the user specified model. X,Y,Z positions in [pc];
    Galactocentric distance in [kpc]

    Parameters
    ----------
    gx_component : str
        Choose from : 'ThinDisk', 'Bulge', 'ThickDisk'
    size : int
        Size of the sample
    model : str
        Current default model is 'McMillan'

    Returns
    -------
    xGx, yGx, zGx, inc, OMEGA, omega : array
        Array of sampled positions in Galactic cartesian coordinates
        centered on the Galactic center and orientations in radians
    """
    # SAMPLE FROM VERY SIMPLE DISTRIBUTION FUNCTIONS
    ###########################################################################
    if gx_component == "ThinDisk":
        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        #######################################################################
        # Assign the radial and vertical positions:
        if model == "double_exp":
            r = sample_exponential_radial(size, 2.5)
            z = sample_exponential_vertical(size, 0.3)

        if model == "sech_squared":
            r = sample_exponential_radial(size, 2.5)
            z = sample_sech_squared(size, 0.3)

        if model == "McMillan":
            r = sample_exponential_radial(size, 2.9)
            z = sample_exponential_vertical(size, 0.3)

        # Assign the azimuthal positions:
        phi = np.random.uniform(0, 2 * np.pi, size)

        # convert to cartesian:
        xGX = r * np.cos(phi)
        yGX = r * np.sin(phi)
        zGX = z

    elif gx_component == "Bulge":

        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        #######################################################################
        # Assign the radial positions:
        if model == "exp_squared":
            r = sample_exponential_square_radial(size, 0.5)
            # Assign the polar positions (uniform in cos(theta):
            theta = np.pi - np.arccos(np.random.uniform(-1, 1, size))

            # Assign the azimuthal positions:
            phi = np.random.uniform(0, 2 * np.pi, size)

            # convert to cartesian:
            xGX = r * np.cos(phi) * np.sin(theta)
            yGX = r * np.sin(phi) * np.sin(theta)
            zGX = r * np.cos(theta)

        elif model == "McMillan":
            r_save = []
            z_save = []
            # sample double exp func and then rejection sample
            while len(r_save) < size:
                rcut = 2.1
                q = 0.5
                r0 = 0.075
                alpha = -1.8
                r = np.random.uniform(0, 5, size * 10)
                z = np.random.uniform(0, 3, size * 10)
                prob = np.random.uniform(0, 1, size * 10)
                sample_func = np.exp(-(r ** 2 + (z / q) ** 2) / rcut ** 2)
                actual_func = (1 + np.sqrt((r ** 2 + (z / q) ** 2)) / r0) ** (
                    alpha
                ) * sample_func
                (indSave,) = np.where(prob < actual_func)
                for ii in indSave:
                    r_save.append(r[ii])
                    z_save.append(z[ii])
            r = np.array(r_save[:size])
            z = np.array(z_save[:size])
            ind_pos_neg = np.random.uniform(0, 1, len(z))
            (ind_negative,) = np.where(ind_pos_neg > 0.5)
            z[ind_negative] = -z[ind_negative]

            # Assign the azimuthal positions:
            phi = np.random.uniform(0, 2 * np.pi, size)

            # convert to cartesian:
            xGX = r * np.cos(phi)
            yGX = r * np.sin(phi)
            zGX = z

    elif gx_component == "ThickDisk":
        # COMPUTE THE POSITION AND ORIENTATION OF EACH BINARY
        ###############################################################################
        # Assign the radial and vertical positions:
        if model == "double_exp":
            r = sample_exponential_radial(size, 2.5)
            z = sample_exponential_vertical(size, 1.0)

        if model == "McMillan":
            r = sample_exponential_radial(size, 3.31)
            z = sample_exponential_vertical(size, 0.9)

        # Assign the azimuthal position of the star
        phi = np.random.uniform(0, 2 * np.pi, size)

        # convert to cartesian:
        xGX = r * np.cos(phi)
        yGX = r * np.sin(phi)
        zGX = z

    # assign an inclination, argument of periapsis, and longitude of ascending node
    inc = np.pi - np.arccos(np.random.uniform(-1, 0, size))
    OMEGA = np.random.uniform(0, 2 * np.pi, size)
    omega = np.random.uniform(0, 2 * np.pi, size)

    return xGX, yGX, zGX, inc, OMEGA, omega
