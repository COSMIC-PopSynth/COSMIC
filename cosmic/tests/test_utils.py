"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
import pandas as pd
import scipy.special as special

import cosmic.utils as utils

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
y = np.random.uniform(0.1,0.2,10)
kstar_single = [[10], [11], [12], [13], [14]]
kstar_double = [10, 14]
x_dat = pd.DataFrame(np.vstack([10*x, 10*f]).T, columns=['x_dat', 'f_dat'])
x_sample = np.vstack([np.random.uniform(0, 1, 10), np.random.uniform(0, 1, 10)]).T
wrong_dict = {'test_wrong_dict' : False}
alive_dict = {'mass_transfer_white_dwarf_to_co' : True, 
             'select_final_state' : True,
             'binary_state' : [0],
             'merger_type' : [-1],
             'LISA_sources' : True}
noLISA_dict = {'mass_transfer_white_dwarf_to_co' : True,
               'select_final_state' : True,
               'binary_state' : [0],
               'merger_type' : [-1],
               'LISA_sources' : False}
false_dict = {'mass_transfer_white_dwarf_to_co' : False,
             'select_final_state' : False,
             'binary_state' : [0, 1, 2, 3],
             'merger_type' : [100],
             'LISA_sources' : False}
conv_dict_true = {'LISA_convergence' : True}
conv_dict_false = {'LISA_convergence' : False}


TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
BPP_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bpp')
BCM_TEST = pd.read_hdf(os.path.join(TEST_DATA_DIR, 'utils_test.hdf'), key='bcm')


IDL_TABULATE_ANSWER = 0.5
MASS_SUM_SINGLE = [41.0, 44.0, 50.0, 132.0, 320.0]
MASS_SUM_MULTIPLE = 301.0
X_TRANS_SUM = -2.7199038e-07  
BW_KNUTH = 0.333
_KNOWN_METHODS = ['mass_transfer_white_dwarf_to_co',
                  'select_final_state',
                  'binary_state',
                  'merger_type',
                  'LISA_sources']


class TestUtils(unittest2.TestCase):
    """`TestCase` for the utilities method
    """
    def test_filter_bpp_bcm(self):
        self.assertRaises(ValueError, utils.filter_bpp_bcm, BCM_TEST, BPP_TEST, wrong_dict)

        bcm_true = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, alive_dict)
        
        self.assertTrue(bcm_true.tphys.all() >= 1.0)
        self.assertTrue(len(bcm_true.loc[bcm_true.sep > 0.0]) >= 1)
        self.assertTrue(len(bcm_true.loc[(bcm_true.RROL_2 > 1)]) >= 1)
        self.assertTrue(bcm_true.porb.all() < 4.0)

        bcm_false = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, false_dict)
        
        self.assertTrue(len(bcm_false.loc[bcm_false.tphys <= 1.0]) > 1)
        self.assertTrue(len(bcm_false.loc[bcm_false.sep == 0.0] > 1))
        self.assertTrue(bcm_false.loc[(bcm_false.RROL_2 > 1)].kstar_2.all()<10)

        bcm_no_LISA = utils.filter_bpp_bcm(BCM_TEST, BPP_TEST, noLISA_dict)
        self.assertTrue(len(bcm_no_LISA.loc[bcm_no_LISA.porb < 4.0]) > 1)
 
    def test_bcm_conv_select(self):
        self.assertRaises(ValueError, utils.filter_bpp_bcm, BCM_TEST, BPP_TEST, wrong_dict)

        bcm_1, bcm_2 = utils.bcm_conv_select(BCM_TEST, BCM_TEST, conv_dict_true)
        self.assertEqual(len(bcm_2), int(len(bcm_1)/2))
        self.assertTrue(bcm_1.porb.all() < 3.0)
        self.assertTrue(bcm_2.porb.all() < 3.0)

        bcm_1_F, bcm_2_F = utils.bcm_conv_select(BCM_TEST[:len(BCM_TEST)-10],\
                                                 BCM_TEST[len(BCM_TEST)-10:],\
                                                 conv_dict_false)
        self.assertEqual(len(BCM_TEST[len(BCM_TEST)-10:]), len(bcm_2_F))
        self.assertTrue(len(bcm_1_F.porb.loc[bcm_1_F.porb > 3.0]) > 1)
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

