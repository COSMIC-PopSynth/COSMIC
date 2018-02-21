"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
import scipy.integrate



IDL_TABULATE_ANSWER = 0.5
class TestSample(unittest2.TestCase):
    """`TestCase` for the aCOSMIC
    """
    def test_idl_tabulate(self):
        # Give this custom integrator a simple integration
        # of a line from x = 0 to 1 and y= 0 to 1
        from aCOSMIC.sample import idl_tabulate
        self.assertAlmostEqual(idl_tabulate(x,f), IDL_TABULATE_ANSWER)
    def test_sample_masses(self):
        self.assertEqual(1, 1)
