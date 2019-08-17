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
PARAMS_INI = os.path.join(TEST_DATA_DIR,'Params.ini')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf5'), key='initC')

init_conds_columns = ['kstar_1', 'kstar_2', 'mass1_binary', 'mass2_binary', 'porb', 'ecc', 'metallicity', 'binfrac', 'tphysf']

INIT_CONDITIONS_NO_BSE_COLUMNS = INIT_CONDITIONS[init_conds_columns]
BPP_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bpp')
BCM_DF = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'unit_tests_results.hdf5'), key='bcm')
BSEFlag_columns = ['neta', 'bwind', 'hewind', 'alpha1', 'lambdaf',
'ceflag', 'tflag', 'ifflag', 'wdflag', 'pisn', 'bhflag', 'nsflag',
'cekickflag', 'cemergeflag', 'cehestarflag',
'mxns', 'pts1', 'pts2', 'pts3',
'ecsn', 'ecsn_mlow', 'aic', 'ussn', 'sigma', 'sigmadiv', 'bhsigmafrac', 'polar_kick_angle',
'beta', 'xi', 'acc2', 'epsnov',
'eddfac', 'gamma', 'bconst', 'ck', 'windflag', 'qcflag', 'eddlimflag']
BSEDict = INIT_CONDITIONS[BSEFlag_columns].to_dict(orient='index')[0]
BSEDict['qcrit_array'] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
BSEDict['natal_kick_array'] = [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]
BSEDict['fprimc_array'] = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]

class TestEvolve(unittest2.TestCase):
    """`TestCase` for the cosmic
    """
    def test_single_evolve_with_table(self):
        # Check that the sample_primary function samples mass correctly
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS,)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False, check_less_precise=True)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False, check_less_precise=True)

    def test_single_evolve_with_dict(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, BSEDict=BSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False, check_less_precise=True)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False, check_less_precise=True)

    def test_single_evolve_with_inifile(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS_NO_BSE_COLUMNS, params=PARAMS_INI, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False, check_less_precise=True)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False, check_less_precise=True)

    def test_single_evolve_with_dict_and_table(self):
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=INIT_CONDITIONS, BSEDict=BSEDict, randomseed=523574)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, BPP_DF, check_dtype=False, check_exact=False, check_less_precise=True)
        pd.testing.assert_frame_equal(EvolvedBinaryBCM, BCM_DF, check_dtype=False, check_exact=False, check_less_precise=True)
