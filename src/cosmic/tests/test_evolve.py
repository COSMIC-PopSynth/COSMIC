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
        self.assertAlmostEqual(kick_info['vsys_2_total'].iloc[0], 15.842697, places=5)        

    def test_vink_wind_metallicity_scaling(self):
        """Vink+2001 OB winds scale as Mdot ~ (Z/zsun)^0.85 on both sides of the bi-stability jump
        but a previous bug in the code had an incorrect scaling. This test checks this is correct always.
        Evolve single stars across their main sequence at a range
        # of metallicities, strip out the L, M and Teff terms of the Vink formula and
        # check that what is left scales with metallicity with a 0.85 exponent"""
        metallicities = np.logspace(-3.5, np.log10(0.03), 10)
        masses = {"cool": 6.0, "hot": 40.0}  # MS Teff below/above the 25 kK boundary

        n_Z = len(metallicities)
        m1 = np.repeat(list(masses.values()), n_Z)
        n_stars = len(m1)
        ibt = InitialBinaryTable.InitialBinaries(
            m1=m1, m2=np.zeros(n_stars), porb=np.zeros(n_stars), ecc=np.zeros(n_stars),
            tphysf=np.full(n_stars, 200.0), kstar1=np.ones(n_stars), kstar2=np.zeros(n_stars),
            metallicity=np.tile(metallicities, len(masses)))

        # Vink winds without LBV winds, no Eddington-dependent exponent, every (short) timestep in bcm
        # (drop bin_num from the shared test dict so that each star keeps its own)
        wind_dict = {k: v for k, v in BSEDict.items() if k != 'bin_num'}
        wind_dict.update(windflag=2, eddlimflag=0, dtp=0.0, pts1=0.001)
        zsun = wind_dict['zsun']
        _, bcm, _, _ = Evolve.evolve(initialbinarytable=ibt, BSEDict=wind_dict, SSEDict=SSEDict)

        ms = bcm[(bcm['kstar_1'] <= 1) & (bcm['deltam_1'] < 0)]
        L, M = ms['lum_1'], ms['mass_1']
        teff = 1000.0 * (1130.0 * L / ms['rad_1']**2)**0.25

        # Vink+2001 log10(Mdot) without the metallicity term (see SSE_mlwind.f)
        x = np.log10(teff / 4.0e4)
        log_mdot_noZ = np.where(
            teff > 25000,
            -6.697 + 2.194 * np.log10(L / 1e5) - 1.313 * np.log10(M / 30.0)
            - 1.226 * np.log10(2.6 / 2.0) + 0.933 * x - 10.92 * x**2,
            -6.688 + 2.210 * np.log10(L / 1e5) - 1.339 * np.log10(M / 30.0)
            - 1.601 * np.log10(1.3 / 2.0) + 1.07 * np.log10(teff / 2.0e4))
        residual = np.log10(-ms['deltam_1']) - log_mdot_noZ

        for i, (side, mass) in enumerate(masses.items()):
            in_regime = (teff > 25000) if side == "hot" else ((teff >= 12500) & (teff <= 25000))
            mean_residual = np.zeros(n_Z)
            for j in range(n_Z):
                this_star = ms['bin_num'] == i * n_Z + j
                use = this_star & in_regime
                # the star should spend (most of) its main sequence on the intended side
                self.assertGreater(use.sum(), 0.5 * this_star.sum(),
                                   msg=f"{mass} Msun star at Z={metallicities[j]:.4g} is not mostly {side}")
                mean_residual[j] = residual[use].mean()

            slope, intercept = np.polyfit(np.log10(metallicities / zsun), mean_residual, 1)
            self.assertAlmostEqual(slope, 0.85, delta=0.01,
                                   msg=f"Vink metallicity exponent on the {side} side is {slope:.3f}")
            self.assertAlmostEqual(intercept, 0.0, delta=0.01)
