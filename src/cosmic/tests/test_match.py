"""Unit test for cosmic
"""

__author__ = 'Katie Breivik <katie.breivik@gmail.com>'

import os
import unittest
import numpy as np
import scipy.integrate
import pandas as pd

from cosmic import Match

np.random.seed(2)
sample = np.random.uniform(0,1,500)

MATCH_TEST = np.log10(1-0.99636185870762237007)

class TestMatch(unittest.TestCase):
    """`TestCase` for the match method
    """
    def test_Match(self):
        # test the match by dividing the sample into to sub samples
        # one containing half of the set and one containing the full set

        dataCm = [sample[:int(len(sample)/2)], sample]
        match, bin_width = Match.match(dataCm)
        self.assertAlmostEqual(match,MATCH_TEST)


