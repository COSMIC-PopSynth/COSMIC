"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd

import cosmic.GW_calcs as GW_calcs

G = 6.67384e-11
c = 2.99792458e8
parsec = 3.08567758e16
Rsun = 6.955e8
Msun = 1.9891e30
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569e7

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')

class TestGWcalcs(unittest2.TestCase):
    """`TestCase` for the GW calcs method
    """

    def test_m_chirp(self):
        # Test the chirp mass calc with GW150914
        m_chirp = GW_calcs.m_chirp(M1, M2)
        self.assertAlmostEqual(m_chirp, MCHIRP_TEST)

    def test_peak_gw_freq(self):
        # Test the peak gw freq with GW150914 
        # w/ ecc=0.5 at 1mHz GW freq
        f_gw_peak = GW_calcs.peak_gw_freq(M1, M2, PORB, ECC)
        self.assertAlmostEqual(f_gw_peak, F_PEAK_TEST)

    def test_peters_gfac(self):
        # Test the g factor with ecc=0.5 and 10 harmonics
        g_factor = GW_calcs.peters_gfac(ECC, N_HARMONIC)
        g_factor_sum = g_factor.sum()
        self.assertAlmostEqual(g_factor_sum, G_FAC_TEST)

 
