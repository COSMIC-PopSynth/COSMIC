"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)

IDL_TABULATE_ANSWER = 0.5

class TestUtils(unittest2.TestCase):
    """`TestCase` for the aCOSMIC
    """
    def test_idl_tabulate(self):
        # Give this custom integrator a simple integration
        # of a line from x = 0 to 1 and y= 0 to 1
        from aCOSMIC.sample import idl_tabulate
        self.assertAlmostEqual(idl_tabulate(x,f), IDL_TABULATE_ANSWER)
