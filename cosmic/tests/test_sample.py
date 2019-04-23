"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import pandas as pd
from cosmic.sample.sampler.independent import Sample
from cosmic.sample.sampler.multidim import MultiDim


SAMPLECLASS = Sample()
MULTIDIMSAMPLECLASS = MultiDim()
TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
TOTAL_SAMPLED_MASS_150_KROUPA93 = 25169.078513680262
TOTAL_SAMPLED_MASS_50_KROUPA93 = 2375.7308462270503
TOTAL_SAMPLED_MASS_KROUPA93 = 41.585324999945854
TOTAL_SAMPLED_MASS_150_SALPETER55 = 17593.2866338113
TOTAL_SAMPLED_MASS_50_SALPETER55 = 1615.8593881185334
TOTAL_SAMPLED_MASS_SALPETER55 = 28.182278678238717
TOTAL_SECONDARY_MASS = 16.15470927770034
N_BINARY_SELECT = 92
PORB = 855.304181220944
THERMAL_ECC_SUM = 5.7488819291685695
UNIFORM_ECC_SUM = 3.58801017672414
CONST_SFR_SUM = 460028.2453521937
BURST_SFR_SUM = 953997.1754647805
BP99_SFR_SUM = 660745.2648859287
BP99_BULGE_SFR_SUM = 811437.0508828465
KSTAR_SOLAR = 1.0
MOE_TOTAL_MASS = 31.134712126322306
METALLICITY_1000 = 0.02
METALLICITY_13000 = 0.02*0.15


class TestSample(unittest2.TestCase):
    """`TestCase` for the cosmic Sample class, which generates several
        independent initial parameters drawn from specified distributions
    """

    def test_sample_primary_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=50.0, primary_max=150.0, primary_model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_150_KROUPA93)

        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=10.0, primary_max=50.0, primary_model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_50_KROUPA93)

        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=0.08, primary_max=5.0, primary_model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_KROUPA93)

    def test_sample_primary_salpeter55(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=50.0, primary_max=150.0, primary_model='salpeter55', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_150_SALPETER55)

        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=10.0, primary_max=50.0, primary_model='salpeter55', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_50_SALPETER55)

        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=0.08, primary_max=5.0, primary_model='salpeter55', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_SALPETER55)

    def test_sample_secondary(self):
        np.random.seed(2)
        # Check that the sample_secondary function samples secondary mass correctly
        a_0 = SAMPLECLASS.sample_secondary(primary_mass = np.arange(10))
        self.assertEqual(np.sum(a_0), TOTAL_SECONDARY_MASS)

    def test_binary_select(self):
        np.random.seed(2)
        # Check that the binary select function chooses binarity properly
        m1_b, m1_s = SAMPLECLASS.binary_select(primary_mass=np.arange(1,100))
        self.assertEqual(len(m1_b), N_BINARY_SELECT)

    def test_sample_porb_abt(self):
        np.random.seed(2)
        # Check that the sample_porb function samples porb properly
        porb = SAMPLECLASS.sample_porb(1.0, 1.0, size=1)
        self.assertAlmostEqual(porb[0], PORB)

    def test_sample_ecc(self):
        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='thermal', size=10)
        self.assertEqual(ecc.sum(), THERMAL_ECC_SUM)

        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(ecc_model='uniform', size=10)
        self.assertEqual(ecc.sum(), UNIFORM_ECC_SUM)

    def test_sample_SFH(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH='const' correctly
        times, met = SAMPLECLASS.sample_SFH(SFH_model='const', component_age=10000.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        np.random.seed(2)
        # Check that the sample SFH function samples SFH='burst' correctly
        times, met = SAMPLECLASS.sample_SFH(SFH_model='burst', component_age=10000.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        times, met = SAMPLECLASS.sample_SFH(SFH_model='delta_burst',\
                                            component_age=10000.0,\
                                            met = 0.02, size=100)
        self.assertEqual(times.sum(), 100*10000.0)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        np.random.seed(2)
        times, met = SAMPLECLASS.sample_SFH(SFH_model='BP99',\
                                            component_age=10000.0,\
                                            met = 0.02, size=100)
        self.assertTrue(np.max(times) <= 13566.0)
        self.assertAlmostEqual(times.sum(), BP99_SFR_SUM)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        np.random.seed(2)
        times, met = SAMPLECLASS.sample_SFH(SFH_model='BP99_bulge',\
                                            component_age=10000.0,\
                                            met = 0.02, size=100)
        self.assertTrue(np.max(times) <= 13566.0)
        self.assertAlmostEqual(times.sum(), BP99_BULGE_SFR_SUM)


    def test_set_kstar(self):
        # Check that the kstar is selected properly
        kstar = SAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)

    def test_Moe_sample(self):
        m1, m2, porb, ecc, total_mass, n_samp = MULTIDIMSAMPLECLASS.initial_sample(rand_seed = 2, size=10, nproc=1)
        self.assertEqual(total_mass, MOE_TOTAL_MASS)

    def test_sample_MultiDim_SFH(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH='const' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SFH_model='const',\
                                                    component_age=10000.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        np.random.seed(2)
        # Check that the sample SFH function samples SFH='burst' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SFH_model='burst',\
                                                    component_age=10000.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SFH_model='delta_burst',\
                                                    component_age=10000.0,\
                                                    met = 0.02, size=100)
        self.assertEqual(times.sum(), 100*10000.0)
        self.assertAlmostEqual(np.mean(met), 0.02)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        np.random.seed(2)
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SFH_model='BP99',\
                                                    component_age=10000.0,\
                                                    met = 0.02, size=100)
        self.assertTrue(np.max(times) <= 13566.0)
        self.assertAlmostEqual(times.sum(), BP99_SFR_SUM)

        # Check that the sample SFH function samples SFH='delta_burst' correctly
        np.random.seed(2)
        times, met = MULTIDIMSAMPLECLASS.sample_SFH(SFH_model='BP99_bulge',\
                                                    component_age=10000.0,\
                                                    met = 0.02, size=100)
        self.assertTrue(np.max(times) <= 13566.0)
        self.assertAlmostEqual(times.sum(), BP99_BULGE_SFR_SUM)

    def test_set_kstar_MultiDim(self):
        # Check that the kstar is selected properly
        kstar = MULTIDIMSAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)


