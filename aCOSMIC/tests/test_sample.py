"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
from aCOSMIC.sample import Sample

SAMPLECLASS = Sample(0.02, size=10)
TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
TOTAL_SAMPLED_MASS_KROUPA93 = 2.138821630578966 
TOTAL_SAMPLED_MASS_SALPETER55 = 1.4513692733524741

class TestSample(unittest2.TestCase):
    """`TestCase` for the aCOSMIC
    """

    def test_sample_primary_kroupa93(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(10, size=10)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_KROUPA93)

    def test_sample_primary_salpeter55(self):
        np.random.seed(2)
        # Check that the sample_primary function samples mass correctly
        a_0, total_sampled_mass = SAMPLECLASS.sample_primary(kstar1_final=10, model='salpeter55', size=10)
        self.assertEqual(total_sampled_mass, TOTAL_SAMPLED_MASS_SALPETER55)
