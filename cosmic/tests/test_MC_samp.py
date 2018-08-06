"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import pandas as pd

import cosmic.MC_samp as MC_samp

x_sun = 8000.0
y_sun = 0.0
z_sun = 25

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
FIXED_POP = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'dat_ThinDisk_11_11_3.h5'), key='bcm')
MASS_FIXED = np.max(pd.read_hdf(os.path.join(TEST_DATA_DIR, 'dat_ThinDisk_11_11_3.h5'), key='totalMass'))
GX_COMPONENT_THIN = 'ThinDisk'
GX_COMPONENT_THICK = 'ThickDisk'
GX_COMPONENT_BULGE = 'Bulge'
GX_COMPONENT_MASS_THIN = 4.32e10
GX_COMPONENT_MASS_THICK = 1.44e10
GX_COMPONENT_MASS_BULGE = 8.9e9
N_GX = 135840169

class TestMC_samp(unittest2.TestCase):
    """`TestCase` for the MC samp method
    """

    def test_mass_weighted_number(self):
        # Check that the total population number is scaled properly by mass
        n_pop = MC_samp.mass_weighted_number(FIXED_POP, MASS_FIXED, GX_COMPONENT_MASS_THIN)
        self.assertEqual(n_pop, N_GX)

    def test_select_component_mass(self):
        # Check that the Galactic component masses are selected properly
        thin_mass = MC_samp.select_component_mass(GX_COMPONENT_THIN)
        self.assertEqual(thin_mass, GX_COMPONENT_MASS_THIN)

        thick_mass = MC_samp.select_component_mass(GX_COMPONENT_THICK)
        self.assertEqual(thick_mass, GX_COMPONENT_MASS_THICK)

        bulge_mass = MC_samp.select_component_mass(GX_COMPONENT_BULGE)
        self.assertEqual(bulge_mass, GX_COMPONENT_MASS_BULGE)

     
