"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample import initialbinarytable
from cosmic import evolve

import warnings
warnings.filterwarnings("ignore")

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
PARAMS_INI = os.path.join(TEST_DATA_DIR,'Params.ini')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf5'), key='initC')
KICK_INITC = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'kick_initial_conditions.h5'), key='initC')

init_conds_columns = initialbinarytable.INITIAL_CONDITIONS_COLUMNS_ALL

INIT_CONDITIONS_NO_BSE_COLUMNS = INIT_CONDITIONS[init_conds_columns]
BPP_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bpp')
BCM_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bcm')
BSEFlag_columns = list(set(evolve.INITIAL_BINARY_TABLE_SAVE_COLUMNS) - set(initialbinarytable.INITIAL_CONDITIONS_COLUMNS_ALL)) 
BSEDict = INIT_CONDITIONS[BSEFlag_columns].to_dict(orient='index')[0]
BSEDict['qcrit_array'] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
BSEDict['natal_kick_array'] = [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]]
BSEDict['fprimc_array'] = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]
BSEDict['grflag'] = 1
BSEDict['don_lim'] = -1
BSEDict['acc_lim'] = -1
BSEDict['wd_mass_lim'] = 0
BSEDict['kick_flag'] = -1

class TestEvolve(unittest.TestCase):
    """`TestCase` for the cosmic
    """
    def test_single_evolve_with_table(self):

        # Check that the sample_primary function samples mass correctly
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_dict(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, BSEDict=BSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_inifile(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, params=PARAMS_INI, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_dict_and_table(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_table(self):
        # Check that the sample_primary function samples mass correctly
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_dict(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, BSEDict=BSEDict, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_inifile(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, params=PARAMS_INI, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_dict_and_table(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_ejection_velocity_pfahl(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=KICK_INITC)

        self.assertAlmostEqual(kick_info['vsys_2_total'].iloc[0], 482.346136, places=5)
