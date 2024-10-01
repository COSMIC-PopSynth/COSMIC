"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

from cosmic.sample import InitialBinaryTable
from cosmic.evolve import Evolve

import warnings

import os
import unittest
import numpy as np
import scipy.integrate
import pandas as pd
import scipy.special as special
import pytest

from cosmic import utils

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
y = np.random.uniform(0.1,0.2,10)
kstar_single = [[10], [11], [12], [13], [14]]
kstar_double = [10, 14]
k1_range = [10, 11, 12]
k2_range = [10, 11, 12]
k1_range_false = np.arange(0,12)
k2_range_false = np.arange(0,12)
x_dat = pd.DataFrame(np.vstack([10*x, 10*f]).T, columns=['x_dat', 'f_dat'])
x_sample = np.vstack([np.random.uniform(0, 1, 10), np.random.uniform(0, 1, 10)]).T
wrong_dict = {'test_wrong_dict' : False}
alive_dict = {'binary_state' : [0]}
noLISA_dict = {'binary_state' : [0]}
false_dict = {'binary_state' : [0,1,2]}
conv_dict_formation = {'pop_select' : 'formation'}
conv_dict_1_SN = {'pop_select' : '1_SN'}
conv_dict_2_SN = {'pop_select' : '2_SN'}
conv_dict_disruption = {'pop_select' : 'disruption'}
conv_dict_final_state = {'pop_select' : 'final_state'}
conv_dict_XRB_form = {'pop_select' : 'XRB_form'}
conv_dict_false = {'pop_select' : 'wrong'}
conv_lim_dict = {"sep" : [10, 5000]}

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
BPP_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bpp')
BCM_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bcm')

IBT = InitialBinaryTable.InitialBinaries(m1=[100.0, 11.8,10**1.5], m2=[85.0, 11.1,21], porb=[10000.0,2211.0,0.1], ecc=[0.65,0.55,0.0], tphysf=[13700.0,13700.0,13700.0], kstar1=[1,1,1], kstar2=[1,1,14], metallicity=[0.005,0.02,0.002], binfrac=[0.5,0.5,0.5])

IDL_TABULATE_ANSWER = 0.5
MASS_SUM_SINGLE = [41.0, 41.6, 42.0, 126.0, 316.0]
MASS_SUM_MULTIPLE = 301.0
X_TRANS_SUM = -2.7199038e-07
BW_KNUTH = 0.333
_KNOWN_METHODS = ['select_final_state',
                  'binary_state']

class TestUtils(unittest.TestCase):
    """`TestCase` for the utilities method
    """
    def test_filter_bin_state(self):
        self.assertRaises(ValueError, utils.filter_bin_state, BCM_TEST, BPP_TEST, wrong_dict, kstar_double, kstar_double)

        bcm_true, bin_state_fraction = utils.filter_bin_state(BCM_TEST, BPP_TEST, alive_dict, k1_range, k2_range)

        self.assertTrue(bcm_true.tphys.all() >= 1.0)

        bcm_false, bin_state_fraction = utils.filter_bin_state(BCM_TEST, BPP_TEST, false_dict, k1_range_false, k2_range_false)


    def test_conv_select(self):
        self.assertRaises(ValueError, utils.conv_select, BCM_TEST, BPP_TEST, [11], [11], wrong_dict, {})

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [11], [11], conv_dict_formation['pop_select'], {})
        self.assertTrue(np.all(conv.sep >= 0))

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [13,14], range(0,15), conv_dict_1_SN['pop_select'], {})
        self.assertEqual(len(conv), 0)

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [13,14], range(0,15), conv_dict_2_SN['pop_select'], {})
        self.assertEqual(len(conv), 0)

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [13,14], range(0,15), conv_dict_disruption['pop_select'], {})
        self.assertEqual(len(conv), 4)

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [11], [11], conv_dict_final_state['pop_select'], {})
        self.assertEqual(len(conv), int(len(BCM_TEST)))

        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [13,14], range(0,15), conv_dict_XRB_form['pop_select'], {})
        self.assertEqual(len(conv), 0)

        self.assertRaises(ValueError, utils.conv_select, BCM_TEST, BPP_TEST, [11], [11], false_dict, {})

    def test_conv_lims(self):
        conv, bin_nums = utils.conv_select(BCM_TEST, BPP_TEST, [11], [11], conv_dict_formation['pop_select'], conv_lim_dict)
        self.assertTrue(conv.loc[conv.bin_num.isin(bin_nums)].sep.max() < conv_lim_dict["sep"][1])
        self.assertTrue(conv.loc[conv.bin_num.isin(bin_nums)].sep.min() > conv_lim_dict["sep"][0])

    def test_idl_tabulate(self):
        # Give this custom integrator a simple integration
        # of a line from x = 0 to 1 and y= 0 to 1
        self.assertAlmostEqual(utils.idl_tabulate(x,f), IDL_TABULATE_ANSWER)

    def test_idl_tabulate_Err(self):
        # Force an error by sending in x = [0]
        self.assertEqual(utils.idl_tabulate(np.array([0]),f), 0)

    def test_min_max_mass_single_kstar(self):
        # Send in a single typed binary with both components
        # being the same type
        for ii in range(len(kstar_single)):
            m_list = utils.mass_min_max_select(kstar_single[ii], kstar_single[ii])
            self.assertEqual(np.sum([m_list]), MASS_SUM_SINGLE[ii])

    def test_min_max_mass_multiple_kstar(self):
        # Send in a range of types for a binary for both components
        m_list = utils.mass_min_max_select(kstar_double, kstar_double)
        self.assertEqual(np.sum([m_list]), MASS_SUM_MULTIPLE)

    def test_param_transform(self):
        # Send in a range of numbers that should have a
        # minimum of 0.0 and a maximum of 1.0 after transformation
        self.assertTrue(np.min(utils.param_transform(f)) >= 0.0)
        self.assertTrue(np.max(utils.param_transform(f)) <= 1.0)

    def test_dat_transform(self):
        # send in DataFrame of 10*x (defined above)
        x_trans = utils.dat_transform(10*x_dat, ['x_dat', 'f_dat'])
        self.assertTrue(np.max(special.expit(x_trans[0])) < 1)
        self.assertTrue(np.min(special.expit(x_trans[0])) > 0)

    def test_dat_un_transform(self):
        # send in sampled data set, which will just be random
        # between -inf to +inf and should be transformed to be between
        # 0 and 10
        x_un_trans = utils.dat_un_transform(special.logit(x_sample), x_dat, ['x_dat', 'f_dat'])
        self.assertTrue(np.min(x_un_trans[0]) >= np.min(x_dat.x_dat))
        self.assertTrue(np.max(x_un_trans[0]) <= np.max(x_dat.x_dat))

    def test_binwidth_selector(self):
        # Check that the Knuth's bw is selected properly
        bw = utils.knuth_bw_selector(np.array([x]))
        self.assertTrue(bw.round(3) == BW_KNUTH)

    def test_error_check(self):
        BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'rtmsflag' : 0, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.05, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.5, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' :[[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.4, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 3, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'rembar_massloss' : 0.5, 'zsun' : 0.02, 'kickflag' : 1, 'grflag' : 1, 'acc_lim' : -1, 'don_lim' : -1}
        filters = {'binary_state': [0], 'timestep_conditions' : 'dtp=None'}
        convergence = {'convergence_params': ['mass_1', 'mass_2', 'sep', 'ecc'], 'pop_select': 'formation',\
                       'match': -5.0, 'convergence_limits' : {"sep" : [0,1000]}, 'match' : -3.0,\
                       'apply_convergence_limits' : True}
        sampling = {'sampling_method': 'multidim', 'SF_start': '13700.0', 'SF_duration' : 0.0, 'metallicity': 0.02, 'keep_singles': False}
        utils.error_check(BSEDict,filters,convergence,sampling)
        utils.error_check(BSEDict)
        assert 1==1

    def test_warning_check(self):
        with pytest.warns(UserWarning, match='At least one of your initial binaries is starting in Roche Lobe Overflow'):
            utils.check_initial_conditions(IBT)

    def test_no_RL_check_for_singles(self):
        """Make sure you don't get a divide by zero error when checking for Roche Lobe Overflow"""
        BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0,
                   'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5,
                   'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1,
                   'acc2': 1.5, 'grflag': 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0,
                   'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0,
                   'natal_kick_array': [[-100.0, -100.0, -100.0, -100.0, 0.0],
                                        [-100.0, -100.0, -100.0, -100.0, 0.0]], 'bhsigmafrac': 1.0,
                   'polar_kick_angle': 90, 'qcrit_array': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                   'cekickflag': 2, 'cehestarflag': 0, 'cemergeflag': 0, 'ecsn': 2.25,
                   'ecsn_mlow': 1.6, 'aic': 1, 'ussn': 0, 'sigmadiv': -20.0, 'qcflag': 5,
                   'eddlimflag': 0, 'fprimc_array': [2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                     2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                     2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0,
                                                     2.0/21.0, 2.0/21.0, 2.0/21.0, 2.0/21.0],
                   'bhspinflag': 0, 'bhspinmag': 0.0, 'rejuv_fac': 1.0, 'rejuvflag': 0, 'htpmb': 1,
                   'ST_cr': 1, 'ST_tide': 1, 'bdecayfac': 1, 'rembar_massloss': 0.5, 'kickflag' : 1,
                   'zsun': 0.014, 'bhms_coll_flag': 0, 'don_lim': -1, 'acc_lim': -1,
                   'rtmsflag': 0, 'wd_mass_lim': 1}

        initial_binaries = InitialBinaryTable.sampler('independent', np.linspace(0, 15, 16), np.linspace(0, 15, 16),
                                                    binfrac_model=0.5, SF_start=10.0,
                                                    SF_duration=0.0, met=0.02, size=10,
                                                    primary_model='kroupa01', ecc_model='sana12', porb_model='sana12',
                                                    keep_singles=True)[0]

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            Evolve.evolve(initialbinarytable=initial_binaries, BSEDict=BSEDict)

    def test_convert_kstar_evol_type(self):
        # convert to string
        bpp = utils.convert_kstar_evol_type(BPP_TEST)
        # convert back and then make sure that it is the same
        bpp = utils.convert_kstar_evol_type(BPP_TEST)
        pd.testing.assert_frame_equal(bpp, BPP_TEST)
