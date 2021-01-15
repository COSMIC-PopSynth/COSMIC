"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import pandas as pd
from ..sample.sampler.independent import Sample
from ..sample.sampler.multidim import MultiDim
from ..sample.sampler.cmc import CMCSample
from ..sample.cmc import elson
from ..sample.initialcmctable import InitialCMCTable
from scipy.optimize import curve_fit
from ..utils import a_from_p

SAMPLECLASS = Sample()
MULTIDIMSAMPLECLASS = MultiDim()
CMCSAMPLECLASS = CMCSample()
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

VR_TEST_ARRAY = np.array([-0.18594816,  0.1400368 ,  0.1523644 ,  0.22576684,  0.08073556,
        0.55498596, -0.67882622,  0.40616153,  0.14757831, -0.08412516,
        0.45161456, -0.8976193 ,  0.44370247,  0.1465378 , -0.4286771 ,
        0.16186843, -0.1133108 ,  0.96282722,  0.67939811,  0.16191544,
       -0.2639312 , -0.21315978, -0.11847003, -0.20368144, -0.32500842,
       -0.21170882, -0.09862926, -0.18438057, -0.24658641, -0.32690077,
        0.89204725, -0.51284319, -0.35446582,  0.14927286, -0.11127869,
        0.14568106,  0.27395772,  0.43082869,  0.25500385, -0.33574652,
        0.37472843, -0.59610097, -0.03829008, -0.16599029, -0.35312626,
        0.19659834, -0.06067012,  0.33754019, -0.48883645,  0.27897062,
        0.16407262,  0.04158386, -0.45989258, -0.08093133, -0.49408438,
       -0.66423923,  1.04430277,  0.3084299 ,  0.56908086, -0.13619488,
        0.21995691,  0.21899305,  0.18831186,  0.00114122, -0.07080843,
       -0.63989728,  0.18849505,  0.1295082 ,  0.33902113, -0.08302271,
        0.52236542, -0.13373999, -0.37063424,  0.27660045, -0.13614144,
       -0.36378358,  0.03581035,  0.21635155,  0.06900961,  0.13764414,
        1.09899499,  0.21926952,  0.39769433, -0.50881124,  0.092561  ,
        0.19426666, -0.09103391, -0.21798363,  0.09434753, -0.07313366,
       -0.13701777,  0.26382516,  0.29047077,  0.06467314, -0.07786409,
        0.08246067,  0.12688664,  0.09449071, -0.19330743, -0.03493555])

VT_TEST_ARRAY = np.array([0.48747186, 0.55155613, 0.56147967, 0.41743075, 0.88614991,
       0.50633149, 0.52669552, 0.71334914, 0.46877172, 0.28577335,
       0.40308517, 0.19699487, 0.55134965, 0.27718976, 0.41963609,
       0.61626828, 0.33090972, 0.40562583, 0.38143321, 0.31186274,
       0.49214035, 0.23397304, 0.73316965, 0.60170198, 0.59117713,
       0.51571931, 0.38856547, 0.36578666, 0.21618999, 0.26064677,
       0.22461443, 0.58505281, 0.16842403, 1.31870369, 0.05625215,
       0.71926501, 0.20665954, 0.49372543, 0.60643738, 0.13477519,
       0.42816799, 0.56503676, 0.27173246, 0.68237931, 0.256907  ,
       0.23324984, 0.5000218 , 0.19305519, 0.86587864, 0.60529773,
       0.59690595, 0.02606873, 0.44770205, 0.65147474, 0.20353073,
       0.16486576, 0.51679633, 0.22066847, 0.28146501, 0.67812263,
       0.16349706, 0.19833249, 0.20864183, 0.16841779, 0.50967491,
       0.10819983, 0.85968885, 0.32615345, 0.43991402, 0.10173283,
       0.42292425, 0.98021201, 0.36943507, 0.41638183, 0.25181444,
       0.29835775, 0.18222665, 0.29207768, 0.26636605, 0.48151551,
       0.07546877, 0.49540014, 0.4531577 , 0.17817517, 0.39682114,
       0.16825046, 0.25753726, 0.12097214, 0.25279382, 0.33493362,
       0.23479711, 0.58643927, 0.37307201, 0.32841265, 0.36744431,
       0.04668889, 0.19594559, 0.05090213, 0.20597884, 0.17267217])

R_TEST_ARRAY = np.array([ 0.24202934, 0.3082052,  0.30985062, 0.31534635, 0.3700414,  0.38275938,
        0.43989968, 0.44484695, 0.47661871, 0.48257492, 0.51624937, 0.53928952,
        0.55181529, 0.58177015, 0.59677683, 0.65081737, 0.66393235, 0.66628891,
        0.69233502, 0.72411303, 0.72940649, 0.74752975, 0.75777753, 0.76806457,
        0.80608474, 0.8286205,  0.84132768, 0.84592651, 0.88850017, 0.9005115,
        0.92527742, 0.94413073, 0.95665949, 0.9776045,  0.9929805,  0.99494626,
        1.00579509, 1.02426239, 1.04900808, 1.06394796, 1.06395337, 1.06518658,
        1.10229627, 1.11968891, 1.12343952, 1.1308054,  1.14555028, 1.14674108,
        1.16170566, 1.16311196, 1.16468701, 1.17437419, 1.2314099,  1.26563713,
        1.29129198, 1.29790148, 1.30477888, 1.31420214, 1.3171364,  1.31715867,
        1.33713795, 1.36281262, 1.37545575, 1.39177066, 1.41388615, 1.42814808,
        1.43161668, 1.47521102, 1.51124964, 1.51892413, 1.56015319, 1.57544702,
        1.62993147, 1.63592578, 1.70015741, 1.71205554, 1.76087382, 1.77551052,
        1.77819171, 1.93336049, 2.10396789, 2.11254154, 2.26123724, 2.33356325,
        2.39180499, 2.45001404, 2.47751094, 2.49881394, 2.768836, 2.91788118,
        2.96571776, 3.00182127, 3.3786418,  3.50975295, 5.58299447, 6.40743644,
        7.05155775, 7.35588129, 9.36791099, 15.57253768])

REFF_TEST_ARRAY = np.array([3.94190562, 5.99895482])

SINGLES_CMC_FITS, BINARIES_CMC_FITS = InitialCMCTable.read(filename=os.path.join(TEST_DATA_DIR, "input_cmc.fits"))

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
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1, qmin = 0.1)
        ind_massive, = np.where(mass1 > 5.0)
        q = mass2[ind_massive]/mass1[ind_massive]
        slope = linear_fit(q)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

    def test_binary_select(self):
        np.random.seed(2)
        # Check that the binary select function chooses binarity properly
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=0.0)
        self.assertEqual(len(m1_b), 0)
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=1.0)
        self.assertEqual(len(m1_b), 99)
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model='vanHaaften')
        self.assertEqual(len(m1_b), N_BINARY_SELECT)

    def test_binary_fraction(self):
        np.random.seed(2)
        # Check that the binary fraction tracking is correct
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model=0.5)
        self.assertEqual(binfrac.max(), 0.5)
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100), binfrac_model='vanHaaften')
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
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1, qmin=0.1)
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
        m1, m2, porb, ecc, mass_singles, mass_binaries, n_singles, n_binaries, binfrac = MULTIDIMSAMPLECLASS.initial_sample(rand_seed = 2, size=10, nproc=1, mp_seeds=[0])
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

class TestCMCSample(unittest.TestCase):
    def test_elson_profile(self):
        np.random.seed(2)
        vr, vt, r = CMCSAMPLECLASS.set_vr_vt_r(N=100, r_max=300, gamma=4)
        np.testing.assert_allclose(VR_TEST_ARRAY, vr, rtol=1e-5)
        np.testing.assert_allclose(VT_TEST_ARRAY, vt, rtol=1e-5)
        np.testing.assert_allclose(R_TEST_ARRAY, r, rtol=1e-5)

    def test_set_reff(self):
        reff = CMCSAMPLECLASS.set_reff(mass=np.array([10.0, 20.0]), metallicity=0.02, params=os.path.join(TEST_DATA_DIR,'Params.ini'))
        np.testing.assert_allclose(REFF_TEST_ARRAY, reff)

    def test_cmc_sampler(self):
        np.random.seed(2)
        # Test generating CMC initial conditions and test saving the output to files
        Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.2, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', cluster_profile='plummer', met=0.014, size=20, params=os.path.join(TEST_DATA_DIR,'Params.ini'), gamma=4, r_max=100, qmin=0.1)
        InitialCMCTable.write(Singles, Binaries, filename="input.hdf5")
        InitialCMCTable.write(Singles, Binaries, filename="input.fits")
        Singles, Binaries = InitialCMCTable.read(filename="input.hdf5")
        Singles, Binaries = InitialCMCTable.read(filename="input.fits")
        # read the test files and compare to the static unit tests files
        pd.testing.assert_frame_equal(Singles, SINGLES_CMC_FITS)
        pd.testing.assert_frame_equal(Binaries, BINARIES_CMC_FITS)
