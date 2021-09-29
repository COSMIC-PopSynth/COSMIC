# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
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

"""`FIRE`
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import os

from cosmic.utils import epanechnikov


__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['SFH', 'positions']

def SFH(n_samp):
    """Samples evolution time and metallicity for sources according
    to Galaxy m12i in the Latte Suite. If you want to use this, you
    must have the file: FIRE.h5 which you need to access from Katie.
    The Latte suite of FIRE-2 cosmological zoom-in
    baryonic simulations of Milky Way-mass galaxies `(Wetzel et al 2016) <http://adsabs.harvard.edu/abs/2016ApJ...827L..23W>`_ 
    part of the Feedback In Realistic Environments (FIRE) simulation 
    project, were run using the Gizmo gravity plus hydrodynamics code in 
    meshless finite-mass (MFM) mode `(Hopkins 2015) <http://adsabs.harvard.edu/abs/2015MNRAS.450...53H>`_  
    and the FIRE-2 physics model `(Hopkins et al 2018) <http://adsabs.harvard.edu/abs/2018MNRAS.480..800H>`_.
    Parameters
    ----------
    n_samp : int
        number of binaries to assign evolution time and metallicity
    Returns
    -------
    timess : array
        array of evolution times (or age) in Myr
    met : array
        array of metallicities (0.02=solar)
    """
    cwd = os.getcwd()
    data_dir = cwd+'/data/'
    FIRE = pd.read_hdf(data_dir+'FIRE.h5')
   
    FIRE_cut = FIRE.loc[:,['age','met']]
     
    dat_sample = FIRE_cut.sample(n=n_samp,replace=True)
    # Convert the times into Myr
    times = np.array(dat_sample.age)*1e3
    # Convert the solar unit metallicity to BSE units
    mets = np.array(dat_sample.met)*0.02

    return times, mets

def positions(times, mets):
    """Assigns the Galactocentric cartesian position by finding
    the star particle in m12i with the closest age and metallicity
    to the times and metallicities provided. The cartesian position
    is sampled with an Epanechnikov kernel with kernel length supplied
    by the Ananke framework `(Sanderson+2019) <http://adsabs.harvard.edu/abs/2018arXiv180610564S>_`
    Parameters
    ----------
    timess : array
        array of evolution times (or age) in Myr
    met : array
        array of metallicities (0.02=solar)
    Returns
    -------
    xGx, yGx, zGx : array
        Cartesian Galactocentric coordinates in kpc
    """
    cwd = os.getcwd()
    data_dir = cwd+'/data/'
    FIRE = pd.read_hdf(data_dir+'FIRE.h5')

    part_data_ind = FIRE.age.searchsorted(times/1e3)
    part_data = FIRE.iloc[part_data_ind-1]
    
    part_data['rGx'] = (part_data.xGx**2 + part_data.yGx**2 + part_data.zGx**2)**0.5
    r_Gx = epanechnikov(part_data.kern_len, np.array(part_data.rGx))
    phi_Gx = np.random.uniform(0,2*np.pi, len(r_Gx))
    theta_Gx = np.pi - np.arccos(np.random.uniform(-1, 1, len(r_Gx)))

    x_Gx = r_Gx * np.sin(phi_Gx) * np.sin(theta_Gx)
    y_Gx = r_Gx * np.cos(phi_Gx) * np.sin(theta_Gx)
    z_Gx = r_Gx * np.cos(theta_Gx)

    return x_Gx, y_Gx, z_Gx
