"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd

import aCOSMIC.lisa_sensitivity as lisa_sens

LISA_SENS_1mHz = 4.739480706393628e-19
LISA_ROOT_PSD_1mHz = 1.0666666666666666e-18

class TestSensitivity(unittest2.TestCase):
    """`TestCase` for the lisa_sensitivity method
    """
    def test_lisa_sensitivity(self):
        lisa_sensitivity_interp = lisa_sens.lisa_sensitivity()
        self.assertTrue(lisa_sensitivity_interp(1e-3) == LISA_SENS_1mHz)

    def test_lisa_root_psd(self):
        lisa_root_psd_interp = lisa_sens.lisa_root_psd()
        self.assertTrue(lisa_root_psd_interp(1e-3) == LISA_ROOT_PSD_1mHz)

