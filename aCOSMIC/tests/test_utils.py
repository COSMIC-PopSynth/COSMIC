"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate

import aCOSMIC.utils as utils

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
kstar_single = [[10], [11], [12], [13], [14]]
kstar_double = [10, 14]

IDL_TABULATE_ANSWER = 0.5
MASS_SUM_SINGLE = [11.0, 20.0, 34.0, 112.0, 330.0]
MASS_SUM_MULTIPLE = 301.0

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

