"""Unit tests for cosmic.checkstate
"""

import unittest
import numpy as np

from cosmic import _evolvebin
from cosmic.checkstate import set_checkstates, CHECKSTATE_COLUMNS, DEFAULT_CONDITIONS
from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable

from test_evolve import BSEDict, SSEDict

import warnings
warnings.filterwarnings("ignore")

KSTAR_1 = 3 * int(np.flatnonzero(CHECKSTATE_COLUMNS == "kstar_1")[0])


def get_arrays():
    """Return copies of the condition arrays currently stored in the Fortran"""
    return (_evolvebin.checkstate_array.checkstate_array.copy(),
            _evolvebin.checkstate_exclude.checkstate_exclude.copy(),
            _evolvebin.checkstate_params.dtp_state.copy())


def evolve_binary(timestep_conditions, dtp=100.0):
    """Evolve a wide 5 + 0.5 Msun binary, which is effectively a single 5 Msun star"""
    binary = InitialBinaryTable.InitialBinaries(
        m1=5.0, m2=0.5, porb=1e5, ecc=0.0, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.014
    )
    # the test BSEDict includes a dtp, which takes priority over the table
    _, bcm, _, _ = Evolve.evolve(initialbinarytable=binary, BSEDict={**BSEDict, "dtp": dtp},
                                 SSEDict=SSEDict, timestep_conditions=timestep_conditions)
    return bcm


class TestSetCheckstates(unittest.TestCase):
    """`TestCase` for storing timestep_conditions in the Fortran arrays
    """
    def tearDown(self):
        set_checkstates([])

    def test_array_shapes(self):
        conditions, exclude, dtp_state = get_arrays()
        self.assertEqual(conditions.shape, DEFAULT_CONDITIONS.shape)
        self.assertEqual(exclude.shape, (dtp_state.size, CHECKSTATE_COLUMNS.size))

    def test_ranges(self):
        # inclusive and strict ends of a chained comparison
        set_checkstates([["4<=kstar_1<10", "dtp=0.5"], ["4<kstar_1<=10", "dtp=0.5"]])
        conditions, _, _ = get_arrays()
        np.testing.assert_array_equal(conditions[0, KSTAR_1:KSTAR_1 + 3], [4, 2, 10])
        np.testing.assert_array_equal(conditions[1, KSTAR_1:KSTAR_1 + 3], [4, 1, 10])

        # the same range given as separate conditions
        set_checkstates([["kstar_1>=4", "kstar_1<10", "dtp=0.5"]])
        conditions, _, _ = get_arrays()
        np.testing.assert_array_equal(conditions[0, KSTAR_1:KSTAR_1 + 3], [4, 2, 10])

    def test_equal_and_not_equal(self):
        set_checkstates([["kstar_1=4", "dtp=0.5"], ["kstar_1<10", "kstar_1!=4", "dtp=0.5"]])
        conditions, exclude, _ = get_arrays()
        np.testing.assert_array_equal(conditions[0, KSTAR_1:KSTAR_1 + 3], [4, 0, 4])
        np.testing.assert_array_equal(conditions[1, KSTAR_1:KSTAR_1 + 3], [-10e30, 6, 10])
        self.assertEqual(exclude[1, KSTAR_1 // 3], 4)

    def test_negative_values(self):
        epoch_1 = 3 * int(np.flatnonzero(CHECKSTATE_COLUMNS == "epoch_1")[0])
        bhspin_1 = 3 * int(np.flatnonzero(CHECKSTATE_COLUMNS == "bhspin_1")[0])

        # negative values in every position of a chained comparison, with and without spaces
        for condition in ["-10<epoch_1<-5", "-10 < epoch_1 < -5", "- 10<epoch_1<- 5"]:
            set_checkstates([[condition, "dtp=0.5"]])
            conditions, _, _ = get_arrays()
            np.testing.assert_array_equal(conditions[0, epoch_1:epoch_1 + 3], [-10, 3, -5])

        # negative values in simple comparisons, including with the value first
        set_checkstates([["epoch_1>=-10", "-5>epoch_1", "bhspin_1!=-0.5", "dtp=0.5"],
                         ["bhspin_1=-0.5", "dtp=0.5"]])
        conditions, exclude, _ = get_arrays()
        np.testing.assert_array_equal(conditions[0, epoch_1:epoch_1 + 3], [-10, 2, -5])
        np.testing.assert_array_equal(conditions[0, bhspin_1:bhspin_1 + 3], [-10e30, 4, 10e30])
        self.assertEqual(exclude[0, bhspin_1 // 3], -0.5)
        np.testing.assert_array_equal(conditions[1, bhspin_1:bhspin_1 + 3], [-0.5, 0, -0.5])

    def test_value_first(self):
        set_checkstates([["4=kstar_1", "dtp=0.5"], ["1<kstar_1", "10>=kstar_1", "dtp=0.5"]])
        conditions, _, _ = get_arrays()
        np.testing.assert_array_equal(conditions[0, KSTAR_1:KSTAR_1 + 3], [4, 0, 4])
        np.testing.assert_array_equal(conditions[1, KSTAR_1:KSTAR_1 + 3], [1, 1, 10])

    def test_dtp(self):
        set_checkstates([["dtp=10"], ["kstar_1=4", "DTP=None"]])
        conditions, _, dtp_state = get_arrays()
        np.testing.assert_array_equal(dtp_state[:3], [10, 13700, -1])
        np.testing.assert_array_equal(conditions[0], DEFAULT_CONDITIONS[0])
        self.assertEqual(_evolvebin.check_dtp.check_dtp, 1)

        set_checkstates([])
        self.assertEqual(_evolvebin.check_dtp.check_dtp, 0)

    def test_invalid_conditions(self):
        invalid = [
            [["kstar_1>1", "kstar_1>2", "dtp=0.5"]],        # lower bound set twice
            [["kstar_1=4", "kstar_1<10", "dtp=0.5"]],       # upper bound set twice
            [["kstar_1!=4", "kstar_1!=5", "dtp=0.5"]],      # two exclusions
            [["kstar_1=4"]],                                # no dtp
            [["kstar_1=4", "dtp=0.5", "dtp=1.0"]],          # dtp set twice
            [["kstar_1=4", "dtp>0.5"]],                     # dtp not set with =
            [["kstar_1=4", "dtp=-1"]],                      # negative dtp
            [["kstar_1=4", "dtp=never"]],                   # non-numeric dtp
            [["not_a_column=4", "dtp=0.5"]],                # invalid column
            [["dtp=0.5"]] * 16,                             # too many conditions
        ]
        for timestep_conditions in invalid:
            with self.assertRaises(ValueError, msg=str(timestep_conditions)):
                set_checkstates(timestep_conditions)

    def test_invalid_conditions_leave_arrays_unchanged(self):
        set_checkstates([["4<=kstar_1<10", "dtp=0.5"]])
        before = get_arrays()
        with self.assertRaises(ValueError):
            set_checkstates([["kstar_1=1", "dtp=0.1"], ["kstar_1=4"]])
        for array_before, array_after in zip(before, get_arrays()):
            np.testing.assert_array_equal(array_before, array_after)


class TestTimestepConditions(unittest.TestCase):
    """`TestCase` for the bcm resolution produced by timestep_conditions
    """
    def tearDown(self):
        set_checkstates([])

    def test_ranges(self):
        bcm = evolve_binary([["kstar_1<2", "dtp=5.0"], ["1<kstar_1<4", "dtp=0.01"],
                             ["kstar_1=4", "dtp=0.0"], ["4<kstar_1<10", "dtp=0.5"]])
        ms = bcm[bcm.kstar_1 == 1].tphys.values
        giant = bcm[bcm.kstar_1.isin([2, 3])].tphys.values
        self.assertAlmostEqual(np.median(np.diff(ms)), 5.0)
        self.assertAlmostEqual(np.median(np.diff(giant)), 0.01)

    def test_last_column(self):
        # teff_2 can never be negative so this condition should never apply
        bcm = evolve_binary([["kstar_1<10", "teff_2<0", "dtp=0.5"]])
        ms = bcm[bcm.kstar_1 == 1].tphys.values
        self.assertAlmostEqual(np.median(np.diff(ms)), 100.0)

    def test_not_equal(self):
        bcm = evolve_binary([["kstar_1<10", "kstar_1!=4", "dtp=0.5"]])
        ms = bcm[bcm.kstar_1 == 1].tphys.values
        self.assertAlmostEqual(np.median(np.diff(ms)), 0.5)
        # core helium burning lasts ~15 Myr so would have ~30 rows at 0.5 Myr
        self.assertLess((bcm.kstar_1 == 4).sum(), 3)

    def test_negative_bounds(self):
        # bhspin_1 is zero for these stars, so only the first condition applies
        bcm = evolve_binary([["-1<bhspin_1<1", "dtp=10.0"], ["bhspin_1<-0.5", "dtp=0.5"]])
        ms = bcm[bcm.kstar_1 == 1].tphys.values
        self.assertAlmostEqual(np.median(np.diff(ms)), 10.0)

    def test_dtp_only_condition(self):
        bcm = evolve_binary([["dtp=10.0"], ["kstar_1=4", "dtp=0.5"]])
        ms = bcm[bcm.kstar_1 == 1].tphys.values
        core_he = bcm[bcm.kstar_1 == 4].tphys.values
        self.assertAlmostEqual(np.median(np.diff(ms)), 10.0)
        self.assertAlmostEqual(np.median(np.diff(core_he)), 0.5)

    def test_switch_to_finer_dtp(self):
        # the giant phases are much shorter than 50 Myr so would be missed
        # if the next output time from the main sequence was kept
        bcm = evolve_binary([["kstar_1<2", "dtp=50.0"], ["kstar_1>=2", "mass_1>2", "dtp=0.01"]])
        giant = bcm[bcm.kstar_1.isin([2, 3, 4])].tphys.values
        self.assertGreater(len(giant), 100)
        self.assertAlmostEqual(np.median(np.diff(giant)), 0.01)

    def test_resume_after_dtp_zero(self):
        bcm = evolve_binary([["kstar_1=4", "dtp=0.0"]], dtp=10.0)
        wd = bcm[bcm.kstar_1 == 11].tphys.values
        self.assertGreater(len(wd), 100)
        self.assertAlmostEqual(np.median(np.diff(wd)), 10.0)
