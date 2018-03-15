"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import pandas as pd
from aCOSMIC.sample import Sample
from aCOSMIC.sample import MultiDimSample

SAMPLECLASS = Sample(0.02, size=10)
MULTIDIMSAMPLECLASS = MultiDimSample(0.02, size=10)
TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
TOTAL_SAMPLED_MASS_150_KROUPA93 = 25169.078513680262
TOTAL_SAMPLED_MASS_50_KROUPA93 = 2375.7308462270503
TOTAL_SAMPLED_MASS_KROUPA93 = 41.585324999945854 
TOTAL_SAMPLED_MASS_150_SALPETER55 = 17593.2866338113
TOTAL_SAMPLED_MASS_50_SALPETER55 = 1615.8593881185334
TOTAL_SAMPLED_MASS_SALPETER55 = 28.182278678238717 
TOTAL_SECONDARY_MASS = 16.15470927770034
N_BINARY_SELECT = 92
HAN_PORB = 29551.837204674266
THERMAL_ECC_SUM = 5.7488819291685695
UNIFORM_ECC_SUM = 3.58801017672414
CONST_SFR_SUM = 460028.2453521937
BURST_SFR_SUM = 953997.1754647805
KSTAR_SOLAR = 1.0
MOE_TOTAL_MASS = 62.735163577852106

class TestSample(unittest2.TestCase):
    """`TestCase` for the aCOSMIC Sample class, which generates several 
        independent initial parameters drawn from specified distributions
    """

    def test_sample_primary_150_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=50.0, primary_max=150.0, model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_150_KROUPA93)

    def test_sample_primary_50_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=10.0, primary_max=50.0, model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_50_KROUPA93)

    def test_sample_primary_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=0.08, primary_max=5.0, model='kroupa93', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_KROUPA93)

    def test_sample_primary_150_salpeter55(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=50.0, primary_max=150.0, model='salpeter55', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_150_SALPETER55)

    def test_sample_primary_50_salpeter55(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=10.0, primary_max=50.0, model='salpeter55', size=100)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_50_SALPETER55)

    def test_sample_primary_salpeter55(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(primary_min=0.08, primary_max=5.0, model='salpeter55', size=100)
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
 
    def test_sample_porb_han(self):
        np.random.seed(2)
        # Check that the sample_porb function samples porb properly
        porb = SAMPLECLASS.sample_porb(1.0, 1.0, model='Han', size=1)
        self.assertEqual(porb[0], HAN_PORB)

    def test_sample_ecc_thermal(self):
        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(model='thermal', size=10)
        self.assertEqual(ecc.sum(), THERMAL_ECC_SUM)

    def test_sample_ecc_uniform(self):
        np.random.seed(2)
        # Check that the sample_ecc function samples ecc properly
        ecc = SAMPLECLASS.sample_ecc(model='uniform', size=10)
        self.assertEqual(ecc.sum(), UNIFORM_ECC_SUM)

    def test_sample_SFH_const(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH correctly
        times = SAMPLECLASS.sample_SFH(model='const', component_age=10000.0, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)

    def test_sample_SFH_burst(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH correctly
        times = SAMPLECLASS.sample_SFH(model='burst', component_age=10000.0, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)

    def test_set_kstar(self):
        # Check that the kstar is selected properly
        kstar = SAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)


    def test_Moe_sample(self):
        np.random.seed(2)
        m1, m2, porb, ecc, total_mass = MULTIDIMSAMPLECLASS.initial_sample(size=10, nproc=1)
        self.assertEqual(total_mass, MOE_TOTAL_MASS)

    def test_sample_SFH_const(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH correctly
        times = MULTIDIMSAMPLECLASS.sample_SFH(model='const', component_age=10000.0, size=100)
        self.assertEqual(times.sum(), CONST_SFR_SUM)

    def test_sample_SFH_burst(self):
        np.random.seed(2)
        # Check that the sample SFH function samples SFH correctly
        times = MULTIDIMSAMPLECLASS.sample_SFH(model='burst', component_age=10000.0, size=100)
        self.assertEqual(times.sum(), BURST_SFR_SUM)

    def test_set_kstar(self):
        # Check that the kstar is selected properly
        kstar = MULTIDIMSAMPLECLASS.set_kstar(pd.DataFrame([1.0, 1.0, 1.0, 1.0, 1.0]))
        self.assertEqual(np.mean(kstar), KSTAR_SOLAR)
