"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

from cosmic.sample import InitialBinaryTable

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd
import scipy.special as special
import pytest

import cosmic.utils as utils

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
alive_dict = {'mass_transfer_white_dwarf_to_co' : True,
             'select_final_state' : True,
             'binary_state' : [0],
             'lisa_sources' : True}
noLISA_dict = {'mass_transfer_white_dwarf_to_co' : True,
               'select_final_state' : True,
               'binary_state' : [0],
               'lisa_sources' : False}
false_dict = {'mass_transfer_white_dwarf_to_co' : False,
             'select_final_state' : False,
             'binary_state' : [0,1,2],
             'lisa_sources' : False}
conv_dict_true = {'lisa_convergence' : True}
conv_dict_false = {'lisa_convergence' : False}


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
BPP_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bpp')
BCM_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bcm')

IBT = InitialBinaryTable.MultipleBinary(m1=[100.0, 11.8,10**1.5], m2=[85.0, 11.1,21], porb=[10000.0,2211.0,0.1], ecc=[0.65,0.55,0.0], tphysf=[13700.0,13700.0,13700.0], kstar1=[1,1,1], kstar2=[1,1,14], metallicity=[0.005,0.02,0.002])

IDL_TABULATE_ANSWER = 0.5
MASS_SUM_SINGLE = [41.0, 41.6, 50.0, 132.0, 320.0]
MASS_SUM_MULTIPLE = 301.0
X_TRANS_SUM = -2.7199038e-07
BW_KNUTH = 0.333
_KNOWN_METHODS = ['mass_transfer_white_dwarf_to_co',
                  'select_final_state',
                  'binary_state',
                  'lisa_sources']

class TestUtils(unittest2.TestCase):
    """`TestCase` for the utilities method
    """
    def test_filter_bpp_bcm(self):
        self.assertRaises(ValueError, utils.filter_bpp_bcm, BCM_TEST, BPP_TEST, wrong_dict, kstar_double, kstar_double)

        bcm_true, bin_state_fraction = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, alive_dict, k1_range, k2_range)

        self.assertTrue(bcm_true.tphys.all() >= 1.0)
        self.assertTrue(len(bcm_true.loc[bcm_true.sep > 0.0]) >= 1)
        self.assertTrue(len(bcm_true.loc[(bcm_true.RROL_2 > 1)]) >= 1)
        self.assertTrue(bcm_true.porb.all() < 5.0)

        bcm_false, bin_state_fraction = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, false_dict, k1_range_false, k2_range_false)

        self.assertTrue(len(bcm_false.loc[bcm_false.tphys <= 1.0]) > 1)
        self.assertTrue(len(bcm_false.loc[bcm_false.sep == 0.0]) > 1)
        self.assertTrue(bcm_false.loc[(bcm_false.RROL_2 > 1)].kstar_2.all()<10)

        bcm_no_LISA, bin_state_fraction = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, noLISA_dict, k1_range, k2_range)
        self.assertTrue(len(bcm_no_LISA.loc[bcm_no_LISA.porb < 4.0]) > 1)

    def test_bcm_conv_select(self):
        self.assertRaises(ValueError, utils.bcm_conv_select, BCM_TEST, BPP_TEST, wrong_dict)

        bcm_1, bcm_2 = utils.bcm_conv_select(BCM_TEST, BCM_TEST, conv_dict_true)
        self.assertEqual(len(bcm_2), int(len(bcm_1)/2))
        self.assertTrue(bcm_1.porb.all() < np.log10(5000))
        self.assertTrue(bcm_2.porb.all() < np.log10(5000))

        bcm_1_F, bcm_2_F = utils.bcm_conv_select(BCM_TEST[:len(BCM_TEST)-10],\
                                                 BCM_TEST,\
                                                 conv_dict_false)
        self.assertEqual(len(BCM_TEST[len(BCM_TEST)-10:]), len(bcm_1_F)-len(bcm_2_F))
        self.assertTrue(len(bcm_1_F.porb.loc[bcm_1_F.porb > np.log10(5000)]) > 1)
        self.assertTrue(len(bcm_2_F.porb.loc[bcm_2_F.porb > 3.0]) > 1)

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
        BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0}
        filters = {'mass_transfer_white_dwarf_to_co': False, 'select_final_state': True, 'binary_state': [0], 'lisa_sources': False}
        convergence = {'lisa_convergence': False}
        utils.error_check(BSEDict,filters,convergence)
        utils.error_check(BSEDict)
        assert 1==1

    def test_warning_check(self):
        with pytest.warns(UserWarning, match='At least one of your initial binaries is starting in Roche Lobe Overflow'):
            utils.check_initial_conditions(IBT)

    def test_convert_kstar_evol_type(self):
        # convert to string
        bpp = utils.convert_kstar_evol_type(BPP_TEST)
        # convert back and then make sure that it is the same
        bpp = utils.convert_kstar_evol_type(BPP_TEST)
        pd.testing.assert_frame_equal(bpp, BPP_TEST)
