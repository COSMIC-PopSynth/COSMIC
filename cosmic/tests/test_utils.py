"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd
import scipy.special as special

import cosmic.utils as utils

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
y = np.random.uniform(0.1,0.2,10)
kstar_single = [[10], [11], [12], [13], [14]]
kstar_double = [10, 14]
x_dat = pd.DataFrame(np.vstack([10*x, 10*f]).T, columns=['x_dat', 'f_dat'])
x_sample = np.vstack([np.random.uniform(0, 1, 10), np.random.uniform(0, 1, 10)]).T

IDL_TABULATE_ANSWER = 0.5
MASS_SUM_SINGLE = [41.0, 44.0, 50.0, 132.0, 320.0]
MASS_SUM_MULTIPLE = 301.0
X_TRANS_SUM = -2.7199038e-07  
BW_KNUTH = 0.333

class TestUtils(unittest2.TestCase):
    """`TestCase` for the utilities method
    """
    def test_idl_tabulate(self):
        # Give this custom integrator a simple integration
        # of a line from x = 0 to 1 and y= 0 to 1
        self.assertAlmostEqual(utils.idl_tabulate(x,f), IDL_TABULATE_ANSWER)
    
    def test_idl_tabulate_Err(self):
        # Force an error by sending in x = [0] 
        self.assertEqual(utils.idl_tabulate(np.array([0]),f), 0)

    def test_min_max_mass_single_kstar(self):
        # Send in a single typed binary with both components
        # being the same type
        for ii in range(len(kstar_single)):
            m_list = utils.mass_min_max_select(kstar_single[ii], kstar_single[ii])
            self.assertEqual(np.sum([m_list]), MASS_SUM_SINGLE[ii])

    def test_min_max_mass_multiple_kstar(self):
        # Send in a range of types for a binary for both components
        m_list = utils.mass_min_max_select(kstar_double, kstar_double)
        self.assertEqual(np.sum([m_list]), MASS_SUM_MULTIPLE)
    
    def test_param_transform(self):
        # Send in a range of numbers that should have a 
        # minimum of 0.0 and a maximum of 1.0 after transformation
        self.assertTrue(np.min(utils.param_transform(f)) >= 0.0)
        self.assertTrue(np.max(utils.param_transform(f)) <= 1.0)         

    def test_dat_transform(self):
        # send in DataFrame of 10*x (defined above)
        x_trans = utils.dat_transform(10*x_dat, ['x_dat', 'f_dat'])
        self.assertTrue(np.max(special.expit(x_trans[0])) < 1)
        self.assertTrue(np.min(special.expit(x_trans[0])) > 0)

    def test_dat_un_transform(self):
        # send in sampled data set, which will just be random
        # between -inf to +inf and should be transformed to be between
        # 0 and 10
        x_un_trans = utils.dat_un_transform(special.logit(x_sample), x_dat, ['x_dat', 'f_dat'])
        self.assertTrue(np.min(x_un_trans[0]) >= np.min(x_dat.x_dat))
        self.assertTrue(np.max(x_un_trans[0]) <= np.max(x_dat.x_dat))

    def test_binwidth_selector(self):
        # Check that the Knuth's bw is selected properly
        bw = utils.knuth_bw_selector(np.array([x]))
        self.assertTrue(bw.round(3) == BW_KNUTH)

