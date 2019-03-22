"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf'), key='initC')
BPP_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bpp')
BCM_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bcm')
BSEDict = {}

class TestEvolve(unittest2.TestCase):
    """`TestCase` for the cosmic
    """
    def test_single_evolve(self):
        # Check that the sample_primary function samples mass correctly
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP.round(4), BPP_DF.round(4), check_dtype=False, check_exact=False, check_less_precise=True)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM.round(4), BCM_DF.round(4), check_dtype=False, check_exact=False, check_less_precise=True)
