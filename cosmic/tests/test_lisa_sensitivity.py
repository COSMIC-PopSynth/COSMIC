"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd

import cosmic.lisa_sensitivity as lisa_sens

LISA_SENS_1mHz = 4.04241209e-21

class TestSensitivity(unittest2.TestCase):
    """`TestCase` for the lisa_sensitivity method
    """
    def test_lisa_characteristic_noise(self):
        lisa_sensitivity_interp = lisa_sens.lisa_characteristic_noise()
        self.assertAlmostEqual(lisa_sensitivity_interp(1e-3), LISA_SENS_1mHz)

