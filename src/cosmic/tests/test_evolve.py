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
BSEDict['acc_lim'] = [-1, -1]
BSEDict['wd_mass_lim'] = 0
SSEDict = {'stellar_engine': 'sse'}
BSEDict['maltsev_mode'] = 0
BSEDict['maltsev_fallback'] = 0.5
BSEDict['maltsev_pf_prob'] = 0.1
BSEDict['kick_flag'] = -1
BSEDict['alpha1'] = [1.0, 1.0]
BSEDict['fryer_fmix'] = 1.0
BSEDict['fryer_mcrit_nsbh'] = 5.75



class TestEvolve(unittest.TestCase):
    """`TestCase` for the cosmic
    """
    #def test_single_evolve_with_table(self):

        # Check that the sample_primary function samples mass correctly
    #    EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
    #        initialbinarytable=INIT_CONDITIONS, randomseed=523574)

    #    pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
    #    pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_dict(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, BSEDict=BSEDict, SSEDict=SSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_inifile(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, params=PARAMS_INI, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_single_evolve_with_dict_and_table(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict, SSEDict=SSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    #def test_multi_evolve_with_table(self):
        # Check that the sample_primary function samples mass correctly
    #    EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
    #        initialbinarytable=INIT_CONDITIONS, n_per_block=100)

    #    pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
    #    pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_dict(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, BSEDict=BSEDict, SSEDict=SSEDict, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_inifile(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, params=PARAMS_INI, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_multi_evolve_with_dict_and_table(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict, SSEDict=SSEDict, randomseed=523574, n_per_block=100)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False)

    def test_fill_in_defaults(self):
        from cosmic.utils import get_default_BSE_settings
        defaults = get_default_BSE_settings()
        missing = ['pts1', 'windflag', 'qcrit_array', 'natal_kick_array', 'alpha1', 'acc_lim']

        # reference: explicitly use the defaults for the missing settings
        full_dict = defaults.copy()
        bpp_ref, bcm_ref, initC_ref, _ = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(), BSEDict=full_dict,
            SSEDict=SSEDict, randomseed=523574)

        # a partial dict should fail without filling in defaults
        partial_dict = {k: v for k, v in full_dict.items() if k not in missing}
        with self.assertRaises(Exception):
            Evolve.evolve(initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(), BSEDict=partial_dict,
                          SSEDict=SSEDict, randomseed=523574)

        # but should match the reference when filling them in
        bpp, bcm, initC, _ = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(), BSEDict=partial_dict,
            SSEDict=SSEDict, randomseed=523574, fill_in_defaults=True)
        print(bpp, bpp_ref)
        pd.testing.assert_frame_equal(bpp, bpp_ref)
        pd.testing.assert_frame_equal(bcm, bcm_ref)
        pd.testing.assert_frame_equal(initC, initC_ref)

        # with no BSEDict at all, every setting should be filled in from the defaults
        _, _, initC, _ = Evolve.evolve(initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(),
                                       randomseed=523574, fill_in_defaults=True)
        _, _, initC_ref, _ = Evolve.evolve(initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(),
                                           BSEDict=defaults, randomseed=523574)
        pd.testing.assert_frame_equal(initC, initC_ref)

    def test_mixed_array_and_flattened_settings(self):
        # the test BSEDict contains both array settings (e.g. alpha1) and their
        # flattened versions (e.g. alpha1_0), which shouldn't lead to duplicate columns
        assert 'alpha1' in BSEDict and 'alpha1_0' in BSEDict
        _, _, initC, _ = Evolve.evolve(initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(),
                                       BSEDict=BSEDict, SSEDict=SSEDict, randomseed=523574)
        self.assertFalse(initC.columns.duplicated().any())

        # conflicting values use the array value
        mixed = {**BSEDict, 'alpha1': [2.0, 3.0], 'alpha1_0': 1.0}
        _, _, initC, _ = Evolve.evolve(initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS.copy(),
                                        BSEDict=mixed, SSEDict=SSEDict, randomseed=523574)
        self.assertFalse(initC.columns.duplicated().any())
        self.assertEqual(initC['alpha1_0'].iloc[0], 2.0)
        self.assertEqual(initC['alpha1_1'].iloc[0], 3.0)

    def test_ejection_velocity_pfahl(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond, kick_info = Evolve.evolve(
            initialbinarytable=KICK_INITC)
        self.assertAlmostEqual(kick_info['vsys_2_total'].iloc[0], 15.813456, places=5)        
