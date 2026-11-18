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

import warnings
warnings.filterwarnings("ignore")


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
INIT_CONDITIONS = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'initial_conditions_for_testing.hdf5'), key='initC')
SSEDict = {"stellar_engine": "sse", "path_to_tracks": "", "path_to_he_tracks": ""}
BSEFlag_columns = list(set(INITIAL_BINARY_TABLE_SAVE_COLUMNS) - set(INITIAL_CONDITIONS_COLUMNS_ALL))
BSEDict = INIT_CONDITIONS[BSEFlag_columns].to_dict(orient='index')[0]
BSEDict['qcrit_array'] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
BSEDict['natal_kick_array'] = [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]]
BSEDict['fprimc_array'] = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                           2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]

# avoid using the same randomseed
del BSEDict["bin_num"], BSEDict["randomseed"]


class TestKick(unittest.TestCase):
    """`TestCase` for the cosmic kick routine
    """

    def test_disberg(self):
        """Test the Disberg kick routine
        """
        # define a simple single star that will result in a kick
        single_star = InitialBinaryTable.InitialBinaries(
            m1=20 + np.random.rand() * 0.01, m2=0.0, porb=0,
            ecc=-1, tphysf=100.0, kstar1=1, kstar2=15,
            metallicity=0.02
        )
        N = 10000
        ibt = single_star.loc[single_star.index.repeat(N)].reset_index()

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using disberg
        BSEDict["kickflag"] = 5
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, SSEDict=SSEDict, nproc=1)
        natal_kicks_disberg = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a lognormal distribution and ensure it matches the expected values
        mu_d, sigma_d = norm.fit(np.log(natal_kicks_disberg))
        self.assertTrue(np.round(mu_d, 1) == 5.6)
        self.assertTrue(np.round(sigma_d, 2) == 0.68)

    def test_hobbs(self):
        """Test the Hobbs kick routine
        """
        # define a simple single star that will result in a kick
        single_star = InitialBinaryTable.InitialBinaries(
            m1=20 + np.random.rand() * 0.01, m2=0.0, porb=0,
            ecc=-1, tphysf=100.0, kstar1=1, kstar2=15,
            metallicity=0.02
        )
        N = 10000
        ibt = single_star.loc[single_star.index.repeat(N)].reset_index()

        # turn off ECSN and bh fallback
        BSEDict["ecsn"] = 0
        BSEDict["ecsn_mlow"] = 0
        BSEDict["bhflag"] = 3

        # evolve using hobbs
        BSEDict["kickflag"] = 1
        _, _, _, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict, SSEDict=SSEDict, nproc=1)
        natal_kicks = kick_info['natal_kick'][kick_info['natal_kick'] != 0.0]

        # fit a maxwellian to the hobbs natal kicks and ensure it matches the expected values
        _, s_hobbs = maxwell.fit(natal_kicks, floc=0.0)
        self.assertTrue(np.round(s_hobbs) - 265 < 5)

    def test_kick_collision(self):
        """Test that a natal kick aimed at the companion destroys it

        A He star + MS binary goes through a SN from a circular orbit. We read
        the pre-SN state from the bpp, then re-evolve with a natal kick chosen
        so that the post-SN relative velocity points almost straight at the
        companion: the stars are approaching and the pericentre of the
        (unbound) post-SN orbit is ~0.5 Rsun, well inside the companion's
        ~4 Rsun radius, so this must end in a collision. The mirrored kick,
        which sends the remnant directly away from the companion, must instead
        disrupt the binary and leave the companion intact.
        """
        G = 1.3271244e11  # km^3 / (Msun s^2), as in kick.f
        rsunkm = 6.96e5
        binary = InitialBinaryTable.InitialBinaries(
            m1=6.0, m2=10.0, porb=5.0, ecc=0.0, tphysf=2.0,
            kstar1=7, kstar2=1, metallicity=0.02
        )
        bse = dict(BSEDict, kickflag=1)

        # first pass: get the state of the binary right before the SN
        bpp, _, _, _ = Evolve.evolve(initialbinarytable=binary, BSEDict=bse,
                                     SSEDict=SSEDict, nproc=1, randomseed=42)
        pre_sn = bpp[bpp["evol_type"] == 15].iloc[0]
        self.assertEqual(pre_sn["kstar_2"], 1)
        self.assertGreater(pre_sn["sep"], 0.0)
        self.assertEqual(pre_sn["ecc"], 0.0)

        # with a mean anomaly of 0, the exploding star sits at (a, 0, 0)
        # relative to its companion and moves with velocity (0, v_orb, 0)
        a = pre_sn["sep"] * rsunkm
        v_orb = np.sqrt(G * (pre_sn["mass_1"] + pre_sn["mass_2"]) / a)

        # target post-SN relative velocities (in-plane, km/s)
        for v_rel, expect_collision in [((-1500.0, 50.0), True),
                                        ((1500.0, 50.0), False)]:
            with self.subTest(v_rel=v_rel):
                kick = np.array([v_rel[0], v_rel[1] - v_orb])
                natal_kick_array = [
                    [np.linalg.norm(kick), 0.0,
                     np.degrees(np.arctan2(kick[1], kick[0])) % 360.0, 0.0, 0.0],
                    [-100.0, -100.0, -100.0, -100.0, 0.0]
                ]
                bpp, _, _, kick_info = Evolve.evolve(
                    initialbinarytable=binary,
                    BSEDict=dict(bse, natal_kick_array=natal_kick_array),
                    SSEDict=SSEDict, nproc=1, randomseed=42
                )
                post_sn = bpp.iloc[bpp["evol_type"].tolist().index(15) + 1]
                self.assertEqual(kick_info["disrupted"].iloc[0], 1)
                if expect_collision:
                    # coalescence, companion destroyed
                    self.assertEqual(post_sn["evol_type"], 6)
                    self.assertEqual(post_sn["kstar_2"], 15)
                else:
                    # disrupted, companion survives
                    self.assertEqual(post_sn["evol_type"], 11)
                    self.assertEqual(post_sn["kstar_2"], 1)
