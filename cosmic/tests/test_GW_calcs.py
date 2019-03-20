"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd

import cosmic.GW_calcs as GW_calcs
from cosmic.lisa_sensitivity import lisa_characteristic_noise as LISA_hc

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
m_in_au = 1.496e+11

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
BCM_DAT = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'GW_dat.h5'), key='gw_obs')
GW150914 = BCM_DAT.iloc[0]
N_HARMONIC = 10
Z = 1.0
D_LUM = 100.0
ECC_ARRAY = np.array([0.0, 0.1, 0.6, 0.75, 0.83, 0.9])

MCHIRP_TEST = 28.25414674924588
F_PEAK_TEST = 0.028749555601321534
G_FAC_TEST = 0.9316274283103428
GW150914_SNR_TEST = 0.2444938519579853
GW150914_F_PEAK_TEST = 0.0287495556013
PSD_SUM_TEST = 1.2881240029e-40
D_CO_TEST = 3353.2
D_LUM_TEST = 6706.4
Z_TEST = 0.0228
F_E_TEST = 11.730543071355376
N_MAX_TEST = np.array([2, 20, 50, 150, 300, 500])
HC_2_CIRC_TEST = 2.558119894716951e-40
HC_2_10_TEST = 4.317193898813392e-42
DA_DT_TEST = -0.009584709398702118
DE_DT_TEST = -1.497435006650899e-11
A_TEST = np.array([0.00677441, 0.00677329, 0.00677217, 0.00677105, 0.00676993, 0.00676882, 0.0067677 , 0.00676659, 0.00676548])
E_TEST = np.array([0.61592   , 0.61565748, 0.61539538, 0.61513372, 0.61487249, 0.61461168, 0.6143513 , 0.61409134, 0.6138318])
T_TEST = np.array([0.55555556, 1.11111111, 1.66666667, 2.22222222, 2.77777778, 3.33333333, 3.88888889, 4.44444444, 5.        ])
SNR_TEST = 0.41087030269445507


class TestGWcalcs(unittest2.TestCase):
    """`TestCase` for the GW calcs method
    """

    def test_sep_from_p(self):
        sep = GW_calcs.sep_from_p(1, 1, 0)
        self.assertEqual(sep, 1.0)

    def test_comoving_distance(self):
        # Test the comoving distance calculator
        d_co = GW_calcs.comoving_distance(Z)
        self.assertAlmostEqual(np.round(d_co, 1), D_CO_TEST)

    def test_comoving_distance(self):
        # Test the luminosity distance calculator
        d_lum = GW_calcs.luminosity_distance(Z)
        self.assertAlmostEqual(np.round(d_lum, 1), D_LUM_TEST)
  
    def test_z_from_lum_distance(self):
        # Test the calculator to get z from d_lum
        z = GW_calcs.z_from_lum_distance(D_LUM)
        self.assertAlmostEqual(np.round(z, 4), Z_TEST)

    def test_m_chirp(self):
        # Test the chirp mass calc with a single binary
        m_chirp = GW_calcs.m_chirp(GW150914.mass_1, GW150914.mass_2)
        self.assertAlmostEqual(m_chirp, MCHIRP_TEST)

    def test_peak_gw_freq(self):
        # Test the peak gw freq with a single binary
        f_gw_peak = GW_calcs.peak_gw_freq(m1=GW150914.mass_1*Msun, m2=GW150914.mass_2*Msun, porb=GW150914.porb, ecc=GW150914.ecc)
        self.assertAlmostEqual(f_gw_peak, F_PEAK_TEST)

    def test_peters_g(self):
        # Test the g factor at the n=10 harmonic
        g_factor = GW_calcs.peters_g(GW150914.ecc, N_HARMONIC)
        self.assertAlmostEqual(g_factor, G_FAC_TEST)

    def test_F_e(self):
        # Test the f(d) factor
        f_e = GW_calcs.F_e(GW150914.ecc)
        self.assertAlmostEqual(f_e, F_E_TEST)

    def test_n_max(self):
        n_max = GW_calcs.n_max(ECC_ARRAY)
        self.assertTrue(np.allclose(n_max, N_MAX_TEST))

    def test_hc2_circ(self):
        hc_2_circ = GW_calcs.hc2_circ(D=GW150914.dist/1000, f_orb=1/GW150914.porb, m1=GW150914.mass_1, m2=GW150914.mass_2)
        self.assertAlmostEqual(hc_2_circ, HC_2_CIRC_TEST)

    def test_hc2(self):
        hc_2_10 = GW_calcs.hc2(D=GW150914.dist/1000, f_orb=1/GW150914.porb, m1=GW150914.mass_1, m2=GW150914.mass_2, e=GW150914.ecc, n=10)
        self.assertAlmostEqual(hc_2_10, HC_2_10_TEST)

    def test_da_dt(self):
        da_dt = GW_calcs.da_dt(GW150914.mass_1*Msun, GW150914.mass_2*Msun, GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2)*m_in_au, GW150914.ecc)
        self.assertAlmostEqual(da_dt, DA_DT_TEST)

    def test_da_dt(self):
        de_dt = GW_calcs.de_dt(GW150914.mass_1*Msun, GW150914.mass_2*Msun, GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2)*m_in_au, GW150914.ecc)
        self.assertAlmostEqual(de_dt, DE_DT_TEST)

    def test_peters_evolution(self):
        a, e, t = GW_calcs.peters_evolution(GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2), GW150914.ecc, GW150914.mass_1, GW150914.mass_2, 5*sec_in_year, 10)
        self.assertTrue(np.allclose(a, A_TEST))
        self.assertTrue(np.allclose(e, E_TEST))
        self.assertTrue(np.allclose(t, T_TEST))
        self.assertTrue(np.min(e) > 0.0)
        self.assertTrue(np.min(a) > 0.0)
        self.assertTrue(np.max(e) < 1.0)
        self.assertTrue(np.max(t) <= 5.0)

    def test_snr_calc(self):
        LISA_noise = LISA_hc()
        t_obs = 5*sec_in_year
        z = GW_calcs.z_from_lum_distance(GW150914.dist/1000)
        a = GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2)
        e = GW150914.ecc
        a_evol, e_evol, t_evol = GW_calcs.peters_evolution(GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2), GW150914.ecc, GW150914.mass_1, GW150914.mass_2, 5*sec_in_year, 1e3)
        f_orb = 1/(GW_calcs.p_from_a(a, GW150914.mass_1, GW150914.mass_2)*sec_in_year)
        f_orb_evol = 1/(GW_calcs.p_from_a(a_evol, GW150914.mass_1, GW150914.mass_2) * sec_in_year)
        t_evol_log = np.logspace(-6, np.log10(t_obs), 5000)
        forb_evol_log = np.interp(t_evol_log, xp = t_evol * sec_in_year, fp = f_orb_evol)
        e_evol_log = np.interp(t_evol_log, xp = t_evol * sec_in_year, fp = e_evol)
        n_harm = int(GW_calcs.n_max(np.array([e]))[0])
        h_c_squared = []
        freqs = []
        for n in range(1,n_harm):
            h_c_squared.append(GW_calcs.hc2(GW150914.mass_1, GW150914.mass_2, forb_evol_log, GW150914.dist/1000, e_evol_log, n))
            freqs.append(forb_evol_log*n)
        h_c_squared = np.array(h_c_squared)
        freqs = np.array(freqs)
        snr = GW_calcs.snr_calc(freqs*(1+z), h_c_squared**0.5, LISA_noise(freqs*(1+z)), 1)
        self.assertAlmostEqual(snr, SNR_TEST)

    def test_snr_chirping(self):
        snr = GW_calcs.snr_chirping(GW150914.mass_1, GW150914.mass_2, GW_calcs.sep_from_p(GW150914.porb/sec_in_year, GW150914.mass_1, GW150914.mass_2), GW150914.ecc, GW150914.dist/1000, 5)
        self.assertAlmostEqual(snr, SNR_TEST)
         
