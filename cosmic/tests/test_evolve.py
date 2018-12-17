"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable

bpp_columns = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2' , 'sep', 'ecc', 'RROL_1', 'RROL_2', 'evol_type']

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
BPP_ARRAY = np.load(os.path.join(TEST_DATA_DIR, 'bpp_array_ind_sampling.npy'))
INIT_CONDITIONS = np.load(os.path.join(TEST_DATA_DIR, 'init_conditions_ind_sampling.npy'))

bppDF = pd.DataFrame(BPP_ARRAY, columns=bpp_columns, index=[1] * len(BPP_ARRAY))
bppDF['bin_num'] = 1

SINGLE_BINARY = InitialBinaryTable.SingleBinary(m1=INIT_CONDITIONS[2],
                                                m2=INIT_CONDITIONS[3],
                                                porb=INIT_CONDITIONS[4],
                                                ecc=INIT_CONDITIONS[5],
                                                tphysf=INIT_CONDITIONS[7],
                                                kstar1=INIT_CONDITIONS[0],
                                                kstar2=INIT_CONDITIONS[1],
                                                metallicity=0.02)

BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'CK': -1000, 'bwind': 0.0, 'lambdaf': -1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 2, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'ppsn' : 1}

class TestEvolve(unittest2.TestCase):
    """`TestCase` for the cosmic
    """
    def test_single_evolve(self):
        # Check that the sample_primary function samples mass correctly
        EvolvedBinaryBPP, EvolvedBinaryBCM, initCond = Evolve.evolve(
            initialbinarytable=SINGLE_BINARY, BSEDict=BSEDict,
            idx=1)

        pd.testing.assert_frame_equal(EvolvedBinaryBPP, bppDF, check_dtype=False)
