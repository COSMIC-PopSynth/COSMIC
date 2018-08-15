"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd

import cosmic.gxreal as gxreal

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
GX_MODEL_MCMILLAN = 'McMillan'
DAT_LIST = ['mass_1', 'mass_2', 'porb', 'ecc']
N_SAMP = 10

class TestGxreal(unittest2.TestCase):
    """`TestCase` for the gxreal method
    """

    def test_compute_n_sample(self):
        # Check that the total population number is scaled properly by mass
        gx_real = gxreal.GxReal(FIXED_POP, MASS_FIXED, GX_MODEL_MCMILLAN, GX_COMPONENT_THIN, DAT_LIST)
        gx_real.n_samp = gx_real.compute_n_sample()

        self.assertEqual(gx_real.n_samp, N_GX)

    def test_sample_population(self):
        gx_real = gxreal.GxReal(FIXED_POP, MASS_FIXED, GX_MODEL_MCMILLAN, GX_COMPONENT_THIN, DAT_LIST)
        gx_real.n_samp = N_SAMP
        gx_real.gxrealization = gx_real.sample_population()

        self.assertTrue(np.min(gx_real.gxrealization.mass_1) >= 0.0)
        self.assertTrue(np.min(gx_real.gxrealization.mass_2) >= 0.0)
        self.assertTrue(np.min(gx_real.gxrealization.porb) >= 0.0) 
        self.assertTrue(np.min(gx_real.gxrealization.ecc) >= 0.0)
        
        self.assertTrue(np.max(gx_real.gxrealization.mass_1) <= np.max(FIXED_POP.mass_1) + np.max(FIXED_POP.mass_1)/20.0)
        self.assertTrue(np.max(gx_real.gxrealization.mass_2) <= np.max(FIXED_POP.mass_2) + np.max(FIXED_POP.mass_2)/20.0)
        self.assertTrue(np.max(gx_real.gxrealization.porb) <= np.max(FIXED_POP.porb) + np.max(FIXED_POP.porb)/50.0)
        self.assertTrue(np.max(gx_real.gxrealization.ecc) <= np.max(FIXED_POP.ecc) + np.max(FIXED_POP.ecc)/10.0) 
