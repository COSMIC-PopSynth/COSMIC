"""Unit test for cosmic kick routine
"""

__author__ = 'Tom Wagg <tomjwagg@gmail.com>'

import os
import unittest
import numpy as np
from scipy.stats import maxwell, norm
import pandas as pd

from cosmic.sample.initialbinarytable import InitialBinaryTable, INITIAL_CONDITIONS_COLUMNS_ALL
from cosmic.evolve import Evolve, INITIAL_BINARY_TABLE_SAVE_COLUMNS


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf5'), key='initC')
BSEFlag_columns = list(set(INITIAL_BINARY_TABLE_SAVE_COLUMNS) - set(INITIAL_CONDITIONS_COLUMNS_ALL))
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


class TestKick(unittest.TestCase):
    """`TestCase` for the cosmic kick routine
    """

    def test_disberg(self):
        """Test the Disberg kick routine
        """
        # define a binary that will result in a kick
        single_binary = InitialBinaryTable.InitialBinaries(
            m1=20, m2=0.1, porb=1e20,
            ecc=0.0, tphysf=200.0, kstar1=1, kstar2=1,
            metallicity=0.02
        )
        N = 10000
        ibt = single_binary.loc[single_binary.index.repeat(N)].reset_index()

        # avoid using the same randomseed
        del BSEDict["bin_num"], BSEDict["randomseed"]

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using disberg
        BSEDict["kickflag"] = 5
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, nproc=1)
        natal_kicks_disberg = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a lognormal distribution and ensure it matches the expected values
        mu_d, sigma_d = norm.fit(np.log(natal_kicks_disberg))
        self.assertTrue(np.round(mu_d, 2) == 5.61)
        self.assertTrue(np.round(sigma_d, 2) == 0.69)

    def test_hobbs(self):
        """Test the HObbs kick routine
        """
        # define a binary that will result in a kick
        single_binary = InitialBinaryTable.InitialBinaries(
            m1=20, m2=0.1, porb=1e20,
            ecc=0.0, tphysf=200.0, kstar1=1, kstar2=1,
            metallicity=0.02
        )
        N = 10000
        ibt = single_binary.loc[single_binary.index.repeat(N)].reset_index()

        # avoid using the same randomseed
        del BSEDict["bin_num"], BSEDict["randomseed"]

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using hobbs
        BSEDict["kickflag"] = 1
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, nproc=1)
        natal_kicks = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a maxwellian to the hobbs natal kicks and ensure it matches the expected values
        _, s_hobbs = maxwell.fit(natal_kicks, floc=0.0)
        self.assertTrue(np.round(s_hobbs, 0) == 265)
