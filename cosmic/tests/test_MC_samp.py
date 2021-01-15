"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import pandas as pd

from .. import MC_samp

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
SECH_SQUARED_SCALE_HEIGHT = 0.3/2
N_SAMP = 1000000
SECH_SQUARED_VAR_TEST = SECH_SQUARED_SCALE_HEIGHT**2*np.pi**2/3.0
SECH_SQUARED_MEAN_TEST = 0.0
EXP_RADIAL_SCALE_HEIGHT = 0.3
EXP_RADIAL_MEAN_TEST = EXP_RADIAL_SCALE_HEIGHT
EXP_RADIAL_VAR_TEST = EXP_RADIAL_SCALE_HEIGHT**2
EXP_VERT_SCALE_HEIGHT = 0.3
EXP_VERT_MEAN_TEST = 0.0
EXP_VERT_VAR_TEST = 2*EXP_RADIAL_SCALE_HEIGHT**2
EXP_SQUARE_RAD_SCALE_HEIGHT = 0.3
EXP_SQUARE_RAD_MEAN_TEST = EXP_SQUARE_RAD_SCALE_HEIGHT/np.sqrt(np.pi)
EXP_SQUARE_RAD_VAR_TEST = EXP_SQUARE_RAD_SCALE_HEIGHT*((2+np.pi)*np.abs(EXP_SQUARE_RAD_SCALE_HEIGHT) - 4*EXP_SQUARE_RAD_SCALE_HEIGHT)/(2*np.pi)
MCMILLAN_THINDISK_MEAN_TEST_XY = 2.9
MCMILLAN_THINDISK_MEAN_TEST_Z = 0.0
SECHSQUARE_THINDISK_MEAN_TEST_XY = 2.5
SECHSQUARE_THINDISK_MEAN_TEST_Z = 0.0
DOUBLEEXP_THINDISK_MEAN_TEST_XY = 2.5
DOUBLEEXP_THINDISK_MEAN_TEST_Z = 0.0
DOUBLEEXP_BULGE_MEAN_TEST_XYZ  = 0.5/np.sqrt(np.pi)
MCMILLAN_BULGE_MEAN_TEST_XY = 0.5
MCMILLAN_THICKDISK_MEAN_TEST_XY = 3.31
MCMILLAN_THICKDISK_MEAN_TEST_Z = 0.0
DOUBLEEXP_THICKDISK_MEAN_TEST_XY = 2.5
DOUBLEEXP_THICKDISK_MEAN_TEST_Z = 0.0

class TestMC_samp(unittest.TestCase):
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

    def test_sample_sech_squared(self):
        # Check that the mean & variance of the sampled distribution matches
        # the analytical mean & variance of the sech squared distribution
        sech_squared_sample = MC_samp.sample_sech_squared(N_SAMP, 2*SECH_SQUARED_SCALE_HEIGHT)
        dist_mean = np.mean(sech_squared_sample)
        self.assertTrue(np.abs(dist_mean - SECH_SQUARED_MEAN_TEST) < 1e-3)
        dist_var = np.var(sech_squared_sample)
        self.assertTrue(np.abs(dist_var - SECH_SQUARED_VAR_TEST) < 1e-3)

    def test_sample_exponential_radial(self):
        # Check that the mean & variance of the sampled distribution matches
        # the analytical mean & variance of the radial exponential distribution
        exp_radial_sample = MC_samp.sample_exponential_radial(N_SAMP, EXP_RADIAL_SCALE_HEIGHT)
        dist_mean = np.mean(exp_radial_sample)
        self.assertTrue(np.abs(dist_mean - (EXP_RADIAL_MEAN_TEST)) < 1e-3)
        dist_var = np.var(exp_radial_sample)
        self.assertTrue(np.abs(dist_var - EXP_RADIAL_VAR_TEST) < 1e-3)

    def test_sample_exponential_vertical(self):
        # Check that the mean & variance of the sampled distribution matches
        # the analytical mean & variance of the radial exponential distribution
        exp_vertical_sample = MC_samp.sample_exponential_vertical(N_SAMP, EXP_VERT_SCALE_HEIGHT)
        dist_mean = np.mean(exp_vertical_sample)
        self.assertTrue(np.abs(dist_mean - (EXP_VERT_MEAN_TEST)) < 1e-3)
        dist_var = np.var(exp_vertical_sample)
        self.assertTrue(np.abs(dist_var - EXP_VERT_VAR_TEST) < 1e-3)

    def test_sample_exponential_square_radial(self):
        # Check that the mean & variance of the sampled distribution matches
        # the analytical mean & variance of the radial exponential squared distribution
        exp_square_radial_sample = MC_samp.sample_exponential_square_radial(N_SAMP, EXP_SQUARE_RAD_SCALE_HEIGHT)
        dist_mean = np.mean(exp_square_radial_sample)
        self.assertTrue(np.abs(dist_mean - (EXP_SQUARE_RAD_MEAN_TEST)) < 1e-3)
        dist_var = np.var(exp_square_radial_sample)
        self.assertTrue(np.abs(dist_var - EXP_SQUARE_RAD_VAR_TEST) < 1e-3)

    def test_galactic_positions(self):
        # Check that the galactic positions are sampled properly
        np.random.seed(2)
        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('ThinDisk', N_SAMP, model='McMillan')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - MCMILLAN_THINDISK_MEAN_TEST_XY) < 1e-2)
        self.assertTrue(np.abs(np.mean(zGX) - MCMILLAN_THINDISK_MEAN_TEST_Z) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('ThinDisk', N_SAMP, model='sech_squared')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - SECHSQUARE_THINDISK_MEAN_TEST_XY) < 1e-2)
        self.assertTrue(np.abs(np.mean(zGX) - SECHSQUARE_THINDISK_MEAN_TEST_Z) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('ThinDisk', N_SAMP, model='double_exp')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - DOUBLEEXP_THINDISK_MEAN_TEST_XY) < 1e-2)
        self.assertTrue(np.abs(np.mean(zGX) - DOUBLEEXP_THINDISK_MEAN_TEST_Z) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('Bulge', N_SAMP, model='exp_squared')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2+zGX**2)**0.5) - DOUBLEEXP_BULGE_MEAN_TEST_XYZ) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('Bulge', 50000, model='McMillan')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - MCMILLAN_BULGE_MEAN_TEST_XY) < 1e-1)
        self.assertTrue(np.abs(np.mean(zGX) - 0.0 < 1e-2))
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('ThickDisk', N_SAMP, model='McMillan')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - MCMILLAN_THICKDISK_MEAN_TEST_XY) < 1e-2)
        self.assertTrue(np.abs(np.mean(zGX) - MCMILLAN_THICKDISK_MEAN_TEST_Z) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)

        xGX, yGX, zGX, inc, OMEGA, omega = MC_samp.galactic_positions('ThinDisk', N_SAMP, model='double_exp')
        self.assertTrue(np.abs(np.mean((xGX**2+yGX**2)**0.5) - DOUBLEEXP_THICKDISK_MEAN_TEST_XY) < 1e-2)
        self.assertTrue(np.abs(np.mean(zGX) - DOUBLEEXP_THICKDISK_MEAN_TEST_Z) < 1e-2)
        self.assertTrue(np.abs(np.mean(inc) - 1.0) < 1e-2)
        self.assertTrue(np.abs(np.mean(OMEGA) - np.pi) < 1e-2)
        self.assertTrue(np.abs(np.mean(omega) - np.pi) < 1e-2)
