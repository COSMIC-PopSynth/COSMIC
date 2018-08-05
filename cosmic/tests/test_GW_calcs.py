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
BCM_DAT = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'GW_dat.h5'), key='gw_obs')
GW150914 = BCM_DAT.iloc[0]
N_HARMONIC = 10

MCHIRP_TEST = 28.25414674924588
F_PEAK_TEST = 0.028749555601321534
G_FAC_TEST = 0.313913837667
GW150914_SNR_TEST = 0.244493738333
GW150914_F_PEAK_TEST = 0.0287495556013
PSD_SUM_TEST = 1.2881240029e-40


class TestGWcalcs(unittest2.TestCase):
    """`TestCase` for the GW calcs method
    """

    def test_m_chirp(self):
        # Test the chirp mass calc with a single binary
        m_chirp = GW_calcs.m_chirp(GW150914.mass_1, GW150914.mass_2)
        self.assertAlmostEqual(m_chirp, MCHIRP_TEST)

    def test_peak_gw_freq(self):
        # Test the peak gw freq with a single binary
        f_gw_peak = GW_calcs.peak_gw_freq(m1=GW150914.mass_1*Msun, m2=GW150914.mass_2*Msun, porb=GW150914.porb, ecc=GW150914.ecc)
        self.assertAlmostEqual(f_gw_peak, F_PEAK_TEST)

    def test_peters_gfac(self):
        # Test the g factor with 10 harmonics
        g_factor = GW_calcs.peters_gfac(GW150914.ecc, N_HARMONIC)
        g_factor_sum = g_factor.sum()
        self.assertAlmostEqual(g_factor_sum, G_FAC_TEST)

    def test_SNR(self):
        # Test the SNR with observed GWs with dummy orbital periods and eccs
        SNR_dat = GW_calcs.LISA_SNR(BCM_DAT.mass_1, BCM_DAT.mass_2, BCM_DAT.porb, BCM_DAT.ecc, BCM_DAT.dist*1000*parsec, 10, 4*sec_in_year) 
        self.assertAlmostEqual(SNR_dat.SNR.iloc[0], GW150914_SNR_TEST)
        self.assertAlmostEqual(SNR_dat.gw_freq.iloc[0], GW150914_F_PEAK_TEST)

    def test_PSD(self):
        # Test the PSDs with observed GWs with dummy orbital periods and eccs
        PSD_dat = GW_calcs.LISA_PSD(BCM_DAT.mass_1, BCM_DAT.mass_2, BCM_DAT.porb, BCM_DAT.ecc, BCM_DAT.dist*1000*parsec, 10, 4*sec_in_year)
        self.assertAlmostEqual(np.sum(PSD_dat.PSD), PSD_SUM_TEST)
