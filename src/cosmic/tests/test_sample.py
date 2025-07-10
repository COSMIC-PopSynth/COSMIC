"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import pandas as pd
from cosmic.sample import InitialBinaryTable
from cosmic.sample.sampler.independent import Sample
from cosmic.sample.sampler.multidim import MultiDim
from cosmic.sample.sampler.cmc import CMCSample
from cosmic.sample.cmc import elson
from cosmic.sample.initialcmctable import InitialCMCTable
from scipy.optimize import curve_fit
from cosmic.utils import a_from_p, get_porb_norm

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
# Since we integrate up to 0.7 in ecc, the slope is 4.0 rounded at the 0th decimal
THERMAL_SLOPE = 4.0
SANA12_ECC_POWER_LAW = -0.45
SANA12_ECC_POWER_LAW_ROUND = -0.5
MEAN_RAGHAVAN = 4.9
SIGMA_RAGHAVAN = 2.3
CLOSE_BINARY_FRAC = 0.42


N_BINARY_SELECT = 85
VANHAAFTEN_BINFRAC_MAX = 0.9989087986493874
VANHAAFTEN_BINFRAC_MIN = 0.6192803136799157
OFFNER_MASS_RANGES = [(0.075,0.15), (0.15,0.30), (0.3,0.6), (0.75,1.00), (0.85,1.25), (1.00,1.25), (1.6,2.4), (3,5), (5,8), (8,17), (17,50)]
OFFNER_DATA = [0.19, 0.23, 0.30, 0.42, 0.47, 0.50, 0.68, 0.81, 0.89, 0.93, 0.96]
OFFNER_ERRORS = [0.03, 0.02, 0.02, 0.03, 0.03, 0.04, 0.07, 0.06, 0.05, 0.04, 0.04]
MULTIDIM_BINFRAC_MAX = 0.6146916774140262
MULTIDIM_BINFRAC_MIN = 0.13786300908773025
CONST_SFR_SUM = 460028.2453521937
BURST_SFR_SUM = 946002.8245352195
KSTAR_SOLAR = 1.0
MOE_TOTAL_MASS = 20.27926225850954
METALLICITY_1000 = 0.02
METALLICITY_13000 = 0.02*0.15

QMIN_LOWMASS = 0.1
QMIN_HIGHMASS = 0.7
M2MIN_LOWMASS = 0.08
M2MIN_HIGHMASS = 12.0
BINFRAC_LOWMASS = 499139
BINFRAC_HIGHMASS = 2205

KING_TEST_DATA = np.load(os.path.join(TEST_DATA_DIR, "cmc_king_test.npz"))
ELSON_TEST_DATA = np.load(os.path.join(TEST_DATA_DIR, "cmc_elson_test.npz"))
PLUMMER_TEST_DATA = np.load(os.path.join(TEST_DATA_DIR, "cmc_plummer_test.npz"))
R_PLUMMER_TEST_ARRAY, VR_PLUMMER_TEST_ARRAY, VT_PLUMMER_TEST_ARRAY = PLUMMER_TEST_DATA["arr_0"], PLUMMER_TEST_DATA["arr_1"], PLUMMER_TEST_DATA["arr_2"] 
R_ELSON_TEST_ARRAY, VR_ELSON_TEST_ARRAY, VT_ELSON_TEST_ARRAY = ELSON_TEST_DATA["arr_0"], ELSON_TEST_DATA["arr_1"], ELSON_TEST_DATA["arr_2"] 
R_KING_TEST_ARRAY, VR_KING_TEST_ARRAY, VT_KING_TEST_ARRAY = KING_TEST_DATA["arr_0"], KING_TEST_DATA["arr_1"], KING_TEST_DATA["arr_2"]

REFF_TEST_ARRAY = np.array([3.94190562, 5.99895482])

SINGLES_CMC_FITS, BINARIES_CMC_FITS = InitialCMCTable.read(filename=os.path.join(TEST_DATA_DIR, "input_cmc.fits"))
SINGLES_CMC_HDF5, BINARIES_CMC_HDF5 = InitialCMCTable.read(filename=os.path.join(TEST_DATA_DIR, "input_cmc.hdf5"))

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

    hist, bins = np.histogram(data, bins=100, density=True)
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
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='salpeter55', size=10000000)
        ind_massive, = np.where(mass1 > 15.0)
        ind_not_massive, = np.where(mass1 < 15.0)

        np.random.seed(2)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass=mass1, qmin=0.1)
        q = mass2[ind_massive] / mass1[ind_massive]
        slope = linear_fit(q)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

        np.random.seed(2)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass=mass1, qmin=-1)
        q = mass2[ind_not_massive] / mass1[ind_not_massive]
        slope = linear_fit(q)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

        np.random.seed(2)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass=mass1, m2_min=0.08)
        q = mass2[ind_massive] / mass1[ind_massive]
        slope = linear_fit(q)
        self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

        # This is now redundant since you should only sample with qmin or m2_min
        #np.random.seed(2)
        #mass2 = SAMPLECLASS.sample_secondary(primary_mass=mass1, qmin=0.1, m2_min=0.08)
        #q = mass2[ind_massive] / mass1[ind_massive]
        #slope = linear_fit(q)
        #self.assertEqual(np.round(slope, 1), FLAT_SLOPE)

    def test_sample_q(self):
        """Test you can sample different mass ratio distributions"""
        np.random.seed(2)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=10000000)
        for slope in [0, 1, 2]:
            mass2 = SAMPLECLASS.sample_secondary(primary_mass=mass1, q_power_law=slope, qmin=0.0)
            q = mass2 / mass1
            fit_slope = power_law_fit(q)
            self.assertEqual(np.round(fit_slope, 1), slope)

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

        test_fracs = []
        test_errs = []
        primary_mass = np.array([float(x) for x in np.logspace(np.log10(0.08), np.log10(150), num=100000)])
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=primary_mass, binfrac_model='offner22')
        for i in range(len(OFFNER_MASS_RANGES)):
            low, high = OFFNER_MASS_RANGES[i][0], OFFNER_MASS_RANGES[i][1]
            offner_value = OFFNER_DATA[i]
            offner_error = OFFNER_ERRORS[i]
            bins_count = len(m1_b[(m1_b >= low) & (m1_b <= high)])
            singles_count = len(m1_s[(m1_s >= low) & (m1_s <= high)])
            bin_frac = bins_count / (bins_count + singles_count)
            error = abs(offner_value - bin_frac)
            self.assertLess(error, offner_error)

    def test_msort(self):
        np.random.seed(2)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=1000000)
        # Check that qmin_msort and m2_min_msort are workings as expected
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1, qmin=0.1, msort=15, qmin_msort=0.7, m2_min_msort=12)
        ind_light, = np.where(mass1 < 15.0)
        ind_massive, = np.where(mass1 >= 15.0)
        m2_light = mass2[ind_light]
        m2_massive = mass2[ind_massive]
        q_light = mass2[ind_light]/mass1[ind_light]
        q_massive = mass2[ind_massive]/mass1[ind_massive]
        assert m2_light.min() > np.min(mass1[ind_light]) * 0.1
        assert m2_massive.min() > np.min(mass1[ind_massive]) * 0.7
        assert q_light.min() > QMIN_LOWMASS
        assert q_massive.min() > QMIN_HIGHMASS
        # Check that the binary fraction tracking is correct when using msort
        m1_b, m1_s, binfrac, bin_index = SAMPLECLASS.binary_select(primary_mass=mass1, binfrac_model=0.5, msort=15, binfrac_model_msort=0.8)
        assert np.sum(binfrac==0.5) == BINFRAC_LOWMASS
        assert np.sum(binfrac==0.8) == BINFRAC_HIGHMASS

#
    def test_sample_porb(self):
        # next do Sana12
        np.random.seed(4)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=100000)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1, qmin=0.1)
        rad1 = SAMPLECLASS.set_reff(mass=mass1, metallicity=0.02)
        rad2 = SAMPLECLASS.set_reff(mass=mass2, metallicity=0.02)
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1, mass2, rad1, rad2, 'sana12', size=mass1.size
        )
        power_slope = power_law_fit(np.log10(porb))
        self.assertEqual(np.round(power_slope, 2), SANA12_PORB_POWER_LAW)

        # now some custom power laws
        for slope in [-0.5, 0.5, 1]:
            porb,aRL_over_a = SAMPLECLASS.sample_porb(
                mass1, mass2, rad1, rad2, porb_model={
                    "min": 0.15, "max": 5, "slope": slope
                }, size=mass1.size
            )
            power_slope = power_law_fit(np.log10(porb))
            self.assertEqual(np.round(power_slope, 1), slope)

        np.random.seed(5)
        # next do Renzo+19
        m1_high = mass1+15
        rad1_high = SAMPLECLASS.set_reff(mass=m1_high, metallicity=0.02)
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            m1_high, mass2, rad1_high, rad2, 'renzo19', size=m1_high.size
        )
        porb_cut = porb[porb > 2.5]
        power_slope = power_law_fit(np.log10(porb_cut))
        self.assertAlmostEqual(np.round(power_slope, 2), SANA12_PORB_POWER_LAW)

        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1, mass2, rad1, rad2, 'renzo19', size=mass1.size
        )
        ind_log_uniform, = np.where(mass1 <= 15)
        porb_log_uniform = porb[ind_log_uniform]
        uniform_slope = linear_fit(np.log10(porb_log_uniform))
        self.assertEqual(np.round(uniform_slope, 1), FLAT_SLOPE)

        # next check the log uniform
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1, mass2, rad1, rad2, 'log_uniform', size=mass1.size
        )
        power_slope = linear_fit(np.log10(porb))
        sep = a_from_p(porb, mass1, mass2)
        sep = sep[sep > 10]
        uniform_slope = linear_fit(np.log10(sep))
        self.assertEqual(np.round(uniform_slope, 1), FLAT_SLOPE)

        # next check raghavan10
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1, mass2, rad1, rad2, 'raghavan10', size=mass1.size
        )
        log_porb_mean = np.mean(np.log10(porb))
        log_porb_sigma = np.std(np.log10(porb))
        self.assertTrue(np.round(log_porb_mean, 1) >= MEAN_RAGHAVAN-0.15)
        self.assertEqual(np.round(log_porb_sigma, 0), np.round(SIGMA_RAGHAVAN, 0))

        # next check moe19
        from cosmic.utils import get_met_dep_binfrac
        from scipy.interpolate import interp1d
        from scipy.stats import kstest
        metallicity = 0.001
        # this is a metallicity dependent population:
        binfrac = get_met_dep_binfrac(metallicity)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=100000)
        (mass1_binaries, mass_single, binfrac_binaries, binary_index,) = SAMPLECLASS.binary_select(
            mass1, binfrac_model=binfrac,
            )
        mass2_binaries = SAMPLECLASS.sample_secondary(
            primary_mass=mass1_binaries, qmin=0.1
            )
        rad1 = SAMPLECLASS.set_reff(mass=mass1_binaries, metallicity=metallicity)
        rad2 = SAMPLECLASS.set_reff(mass=mass2_binaries, metallicity=metallicity)
        
        
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1_binaries, mass2_binaries, rad1, rad2, 'moe19', size=mass1_binaries.size, met=metallicity
        )
        
        binary_frac = len(porb) / (len(mass1))
        self.assertAlmostEqual(np.round(binary_frac, 2), binfrac)
        
                

    def test_sample_ecc(self):
        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly

        # first sample orbital periods
        np.random.seed(4)
        mass1, total_mass = SAMPLECLASS.sample_primary(primary_model='kroupa01', size=100000)
        mass2 = SAMPLECLASS.sample_secondary(primary_mass = mass1, qmin=0.1)
        rad1 = SAMPLECLASS.set_reff(mass=mass1, metallicity=0.02)
        rad2 = SAMPLECLASS.set_reff(mass=mass2, metallicity=0.02)
        porb,aRL_over_a = SAMPLECLASS.sample_porb(
            mass1, mass2, rad1, rad2, 'sana12', size=mass1.size
        )

        # now we feed aRL_over_a into sample_ecc
        ecc = SAMPLECLASS.sample_ecc(aRL_over_a, ecc_model='thermal', size=mass1.size)
        ecc_cut = ecc[ecc < 0.7]
        slope = linear_fit(ecc_cut)
        self.assertEqual(np.round(slope, 0), THERMAL_SLOPE)

        ecc = SAMPLECLASS.sample_ecc(aRL_over_a, ecc_model='sana12', size=mass1.size)
        ecc_cut = ecc[ecc < 0.91]
        power_slope = power_law_fit(ecc_cut)
        self.assertEqual(np.round(power_slope, 1), SANA12_ECC_POWER_LAW_ROUND)

        ecc = SAMPLECLASS.sample_ecc(aRL_over_a, ecc_model='circular', size=mass1.size)
        self.assertEqual(np.mean(ecc), 0.0)


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
        m1, m2, porb, ecc, single_mass_list, mass_singles, mass_binaries, n_singles, n_binaries, binfrac = MULTIDIMSAMPLECLASS.initial_sample(rand_seed = 2, size=10, nproc=1, mp_seeds=[0])
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

    def test_sampling_binfrac_zero(self):
        # check that you can't sample based on size with a binary fraction of 0.0
        it_fails = False
        try:
            InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                       primary_model='kroupa01', ecc_model='thermal',
                                       porb_model='sana12', binfrac_model=0.0,
                                       SF_start=10.0, SF_duration=0.0, met=0.02,
                                       size=1000)
        except ValueError:
            it_fails = True
        self.assertTrue(it_fails)

    def test_sampling_targets_bad_input(self):
        # check that you get an error if you don't supply a size or total_mass
        it_fails = False
        try:
            InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                       primary_model='kroupa01', ecc_model='thermal',
                                       porb_model='sana12', binfrac_model=0.5,
                                       SF_start=10.0, SF_duration=0.0, met=0.02,
                                       sampling_target="total_mass")
        except ValueError:
            it_fails = True
        self.assertTrue(it_fails)
    
        it_fails = False
        try:
            InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                       primary_model='kroupa01', ecc_model='thermal',
                                       porb_model='sana12', binfrac_model=0.5,
                                       SF_start=10.0, SF_duration=0.0, met=0.02,
                                       total_mass=None, sampling_target="total_mass")
        except ValueError:
            it_fails = True
        self.assertTrue(it_fails)
    
        it_fails = False
        try:
            InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                       primary_model='kroupa01', ecc_model='thermal',
                                       porb_model='sana12', binfrac_model=0.5,
                                       SF_start=10.0, SF_duration=0.0, met=0.02,
                                       size=None, total_mass=None, sampling_target="size")
        except ValueError:
            it_fails = True
        self.assertTrue(it_fails)

    def test_sampling_targets_size(self):
        # check that you can sample based on size
        for size in np.random.randint(100, 1000, size=100):
            samples = InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                                primary_model='kroupa01', ecc_model='thermal',
                                                porb_model='sana12', binfrac_model=0.5,
                                                SF_start=10.0, SF_duration=0.0, met=0.02,
                                                size=size)
            self.assertGreaterEqual(len(samples[0]), size)

    def test_sampling_targets_mass(self):
        # check that you can sample based on total mass
        for mass in np.random.randint(100, 1000, size=100):
            samples = InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                                primary_model='kroupa01', ecc_model='thermal',
                                                porb_model='sana12', binfrac_model=0.5,
                                                SF_start=10.0, SF_duration=0.0, met=0.02,
                                                size=None, total_mass=mass, sampling_target="total_mass")
            self.assertGreaterEqual(samples[1] + samples[2], mass)

    def test_sampling_targets_mass_trimmed(self):
        # check that you can sample based on total mass and that the samples are trimmed properly
        for mass in np.random.randint(100, 1000, size=100):
            samples = InitialBinaryTable.sampler('independent', np.arange(16), np.arange(16),
                                                primary_model='kroupa01', ecc_model='thermal',
                                                porb_model='sana12', binfrac_model=0.5,
                                                SF_start=10.0, SF_duration=0.0, met=0.02,
                                                size=None, total_mass=mass, sampling_target="total_mass",
                                                trim_extra_samples=True)
            self.assertLessEqual(abs(samples[1] + samples[2] - mass), 300)

class TestCMCSample(unittest.TestCase):
    def test_plummer_profile(self):
        np.random.seed(2)
        r, vr, vt = CMCSAMPLECLASS.set_r_vr_vt('plummer',N=100, r_max=300)
        np.testing.assert_allclose(VR_PLUMMER_TEST_ARRAY, vr, rtol=1e-5)
        np.testing.assert_allclose(VT_PLUMMER_TEST_ARRAY, vt, rtol=1e-5)
        np.testing.assert_allclose(R_PLUMMER_TEST_ARRAY, r, rtol=1e-5)

    def test_elson_profile(self):
        np.random.seed(2)
        r, vr, vt = CMCSAMPLECLASS.set_r_vr_vt('elson',N=100, r_max=300, gamma=3)
        np.testing.assert_allclose(VR_ELSON_TEST_ARRAY, vr, rtol=1e-5)
        np.testing.assert_allclose(VT_ELSON_TEST_ARRAY, vt, rtol=1e-5)
        np.testing.assert_allclose(R_ELSON_TEST_ARRAY, r, rtol=1e-5)

    def test_king_profile(self):
        np.random.seed(2)
        r, vr, vt = CMCSAMPLECLASS.set_r_vr_vt('king',N=100, w_0=5)
        np.testing.assert_allclose(VR_KING_TEST_ARRAY, vr, rtol=1e-5)
        np.testing.assert_allclose(VT_KING_TEST_ARRAY, vt, rtol=1e-5)
        np.testing.assert_allclose(R_KING_TEST_ARRAY, r, rtol=1e-5)

    def test_set_reff(self):
        reff = CMCSAMPLECLASS.set_reff(mass=np.array([10.0, 20.0]), metallicity=0.02)
        np.testing.assert_allclose(REFF_TEST_ARRAY, reff)

    def test_cmc_sampler(self):
        np.random.seed(2)
        # Test generating CMC initial conditions and test saving the output to files
        Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.2, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', cluster_profile='plummer', met=0.014, size=20, params=os.path.join(TEST_DATA_DIR,'Params.ini'), gamma=4, r_max=100, qmin=0.1)
        InitialCMCTable.write(Singles, Binaries, filename="input.hdf5")
        InitialCMCTable.write(Singles, Binaries, filename="input.fits")
        Singles, Binaries = InitialCMCTable.read(filename="input.fits")
        # read the test files and compare to the static unit tests files
        pd.testing.assert_frame_equal(Singles, SINGLES_CMC_FITS)
        pd.testing.assert_frame_equal(Binaries, BINARIES_CMC_FITS)
        Singles, Binaries = InitialCMCTable.read(filename="input.hdf5")
        pd.testing.assert_frame_equal(Singles, SINGLES_CMC_HDF5)
        pd.testing.assert_frame_equal(Binaries, BINARIES_CMC_HDF5)
