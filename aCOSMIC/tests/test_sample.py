"""Unit test for aCOSMIC
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest2
import numpy as np
import scipy.integrate
from aCOSMIC.sample import Sample
tmp = Sample(0.02, size=10)

TEST_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')

f  = np.linspace(0,1,10)
x = np.linspace(0,1,10)
np.random.seed(2)


IDL_TABULATE_ANSWER = 0.5
SAMPLE_PRIMARY_ANSWER = np.load(os.path.join(TEST_DATA_DIR, 'Kroupa93.npy'))
print SAMPLE_PRIMARY_ANSWER

class TestSample(unittest2.TestCase):
    """`TestCase` for the aCOSMIC
    """
    from aCOSMIC.sample import Sample
    tmp = Sample(0.02, size=10)

    def test_idl_tabulate(self):
        # Give this custom integrator a simple integration
        # of a line from x = 0 to 1 and y= 0 to 1
        from aCOSMIC.sample import idl_tabulate
        self.assertAlmostEqual(idl_tabulate(x,f), IDL_TABULATE_ANSWER)

    def test_sample_primary(self):
        # Check that the sample_primary function samples mass correctly

        self.assertAlmostEqual(tmp.sample_primary(10, size=10), SAMPLE_PRIMARY_ANSWER)
