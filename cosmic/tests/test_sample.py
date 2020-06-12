"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import pandas as pd
from cosmic.sample.sampler.independent import Sample
from cosmic.sample.sampler.multidim import MultiDim
from scipy.optimize import curve_fit
from cosmic.utils import a_from_p

SAMPLECLASS = Sample()
MULTIDIMSAMPLECLASS = MultiDim()
TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')

#distribution slopes/power laws
KROUPA_93_HI = -2.7
KROUPA_93_MID = -2.2
KROUPA_93_LO = -1.3
KROUPA_01_HI = -2.3
KROUPA_01_LO = -1.3
SALPETER_55 = -2.35
SANA12_PORB_POWER_LAW = -0.55
FLAT_SLOPE = 0.0
THERMAL_SLOPE = 2.0
SANA12_ECC_POWER_LAW = -0.45

N_BINARY_SELECT = 85
VANHAAFTEN_BINFRAC_MAX = 0.9989087986493874
VANHAAFTEN_BINFRAC_MIN = 0.6192803136799157
MULTIDIM_BINFRAC_MAX = 0.6146916774140262
MULTIDIM_BINFRAC_MIN = 0.13786300908773025
CONST_SFR_SUM = 460028.2453521937
BURST_SFR_SUM = 946002.8245352195
KSTAR_SOLAR = 1.0
MOE_TOTAL_MASS = 20.27926225850954
METALLICITY_1000 = 0.02
METALLICITY_13000 = 0.02*0.15

def power_law_fit(data, n_bins=100):
    def line(x, a, b):
        return x*a + b
    def center_bins(bins):
        mid_bin = []
        for bin_lo, bin_hi in zip(bins[:-1], bins[1:]):
            mid_bin.append(bin_lo + (bin_hi-bin_lo)/2)
        return mid_bin

    hist, bins = np.histogram(data, bins=n_bins)
    bins = center_bins(bins)

    xdata = np.log10(bins)
    ydata = np.log10(hist)

    popt, pcov = curve_fit(line, xdata, ydata)

    slope, intercept = popt[0], popt[1]

    return slope

def linear_fit(data):
    def line(x, a, b):
        return x*a + b
    def center_bins(bins):
        mid_bin = []
        for bin_lo, bin_hi in zip(bins[:-1], bins[1:]):
            mid_bin.append(bin_lo + (bin_hi-bin_lo)/2)
        return mid_bin
    
    hist, bins = np.histogram(data, bins=50, density=True)
    bins = center_bins(bins)
    popt, pcov = curve_fit(line, bins, hist)

    slope, intercept = popt[0], popt[1]

    return slope

class TestSample(unittest.TestCase):
    """`TestCase` for the cosmic Sample class, which generates several
        independent initial parameters drawn from specified distributions
    """

    def test_sample_primary_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        mass, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa93', size=1000000)
        mass_hi = mass[mass > 1.0]
        # filter out highest masses because they kill us with the histograms
        mass_hi = mass_hi[mass_hi < 10.0]
        mass_mid = mass[(mass <= 1.0)]
        mass_mid = mass_mid[mass_mid > 0.5]
        mass_lo = mass[mass <= 0.5]

        #few bins for the most massive stars
        power_slope_hi = power_law_fit(mass_hi, n_bins=50)
        power_slope_mid = power_law_fit(mass_mid)
        power_slope_lo = power_law_fit(mass_lo)

        self.assertEqual(np.round(power_slope_hi, 1), KROUPA_93_HI)
        self.assertEqual(np.round(power_slope_mid, 1), KROUPA_93_MID)
        self.assertEqual(np.round(power_slope_lo, 1), KROUPA_93_LO)


    def test_sample_primary_kroupa01(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        mass, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=1000000)
        mass_hi = mass[mass > 0.5]
        # filter out highest masses because they kill us with the histograms
        mass_hi = mass_hi[mass_hi < 10.0]
        mass_lo = mass[mass <= 0.5]

        #few bins for the most massive stars
        power_slope_hi = power_law_fit(mass_hi, n_bins=50)
        power_slope_lo = power_law_fit(mass_lo)

        self.assertEqual(np.round(power_slope_hi, 1), KROUPA_01_HI)
        self.assertEqual(np.round(power_slope_lo, 1), KROUPA_01_LO)

    def test_sample_primary_salpeter55(self):
        np.random.seed(3)
        # Check that the sample_primary function samples mass correctly
        mass, total_mass = SAMPLECLASS.sample_primary(primary_model='salpeter55', size=10000000)
        #filter on mass to get better statistics
        power_slope = power_law_fit(mass[mass < 1.0], n_bins=50)
        self.assertEqual(np.round(power_slope, 2), SALPETER_55)

    def test_sample_secondary(self):
        np.random.seed(2)
        # Check that the sample_secondary function samples secondary mass correctly
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='salpeter55', size=1000000)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1)
        ind_massive, = np.where(mass1 > 5.0)
        q = mass2[ind_massive]/mass1[ind_massive]
        slope = linear_fit(q)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

    def test_binary_select(self):
        np.random.seed(2)
        # Check that the binary select function chooses binarity properly
        m1_b, m1_s, binfrac = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=0.0)
        self.assertEqual(len(m1_b), 0)
        m1_b, m1_s, binfrac = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=1.0)
        self.assertEqual(len(m1_b), 99)
        m1_b, m1_s, binfrac = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model='vanHaaften')
        self.assertEqual(len(m1_b), N_BINARY_SELECT)

    def test_binary_fraction(self):
        np.random.seed(2)
        # Check that the binary fraction tracking is correct
        m1_b, m1_s, binfrac = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=0.5)
        self.assertEqual(binfrac.max(), 0.5)
        m1_b, m1_s, binfrac = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model='vanHaaften')
        self.assertEqual(binfrac.max(), VANHAAFTEN_BINFRAC_MAX)
        self.assertEqual(binfrac.min(), VANHAAFTEN_BINFRAC_MIN)

    def test_sample_ecc(self):
        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='thermal', size=100000)
        slope = linear_fit(ecc)
        self.assertEqual(np.round(slope, 1), THERMAL_SLOPE)

        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='uniform', size=100000)
        slope = linear_fit(ecc)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

        np.random.seed(4)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='sana12', size=1000000)
        power_slope = power_law_fit(ecc)
        self.assertEqual(np.round(power_slope, 2), SANA12_ECC_POWER_LAW)

        np.random.seed(4)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='circular', size=1000000)
        self.assertEqual(np.mean(ecc), 0.0)

    def test_sample_porb(self):
        # next do Sana12
        np.random.seed(4)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=1000000)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1)
        ecc = SAMPLECLASS.sample_ecc(ecc_model='sana12', size=len(mass1))
        porb = SAMPLECLASS.sample_porb(mass1, mass2, ecc, 'sana12', size=len(mass1))
        power_slope = power_law_fit(np.log10(porb))
        self.assertAlmostEqual(np.round(power_slope, 2), SANA12_PORB_POWER_LAW)

        # next do Renzo19
        np.random.seed(4)
        porb_power = SAMPLECLASS.sample_porb(mass1+15, mass2, ecc, 'renzo19', size=len(mass1))
        power_slope = power_law_fit(np.log10(porb_power))
        self.assertEqual(np.round(power_slope, 2), SANA12_PORB_POWER_LAW)
        porb_log_uniform = SAMPLECLASS.sample_porb(mass1, mass2, ecc, 'renzo19', size=len(mass1))
        ind_log_uniform, = np.where(mass1 <= 15)
        porb_log_uniform = porb_log_uniform[ind_log_uniform]
        uniform_slope = linear_fit(np.log10(porb_log_uniform))
        self.assertEqual(np.round(uniform_slope, 1), FLAT_SLOPE)

        np.random.seed(2)
        # Check that the sample_porb function samples porb properly
        porb = SAMPLECLASS.sample_porb(mass1, mass2, ecc, 'log_uniform', size=len(mass1))
        # filter out the weird cuts for RLO and convert to sep
        sep = a_from_p(porb, mass1, mass2)
        sep = sep[sep > 10]
        uniform_slope = linear_fit(np.log10(sep))
        self.assertEqual(np.round(uniform_slope, 1), FLAT_SLOPE)

    def test_sample_SFH(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH='const' correctly
        times, met = SAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                            SF_duration=10000.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        np.random.seed(2)
        # Check that the sample SFH function samples SFH='burst' correctly
        times, met = SAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                            SF_duration=1000.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        times, met = SAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                            SF_duration=0.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), 100*10000.0)
        self.assertAlmostEqual(np.mean(met), 0.02)

    def test_set_kstar(self):
        # Check that the kstar is selected properly
        kstar = SAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)

    def test_Moe_sample(self):
        # Test the multidim sampler and system-by-system binary fraction
        m1, m2, porb, ecc, mass_singles, mass_binaries, n_singles, n_binaries, binfrac = MULTIDIMSAMPLECLASS.initial_sample(rand_seed = 2, size=10, nproc=1)
        self.assertEqual(np.sum(mass_singles), MOE_TOTAL_MASS)
        self.assertAlmostEqual(binfrac.max(), MULTIDIM_BINFRAC_MAX)
        self.assertAlmostEqual(binfrac.min(), MULTIDIM_BINFRAC_MIN)

    def test_sample_MultiDim_SFH(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH='const' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                                    SF_duration=10000.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        np.random.seed(2)
        # Check that the sample SFH function samples SFH='burst' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                                    SF_duration=1000.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SF_start=10000.0,\
                                                    SF_duration=0.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), 100*10000.0)
        self.assertAlmostEqual(np.mean(met), 0.02)

    def test_set_kstar_MultiDim(self):
        # Check that the kstar is selected properly
        kstar = MULTIDIMSAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)


