# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of astro-traj.
#
# astro-traj is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# astro-traj is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with astro-traj.  If not, see <http://www.gnu.org/licenses/>.

"""`sample`
"""

import numpy as np
import math
import random

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'Sample'

class Sample:
    def __init__(self):   
        '''
        initialize samples
        '''

    # sample primary and secondary masses
    def sample_kroupa93(self, size=None):
        '''
        Primary mass follows Kroupa (1993), normalization comes from
        `Hurley (2002) <https://arxiv.org/abs/astro-ph/0201220>`_ between 0.1 and 100 Msun
        '''
        a_0 = np.random.uniform(0.0, 0.9999797, size)
        low_cutoff = 0.740074, high_cutoff=0.908422


        lowIdx, = np.where(a_0 <= low_cutoff)
        midIdx, = np.where((a_0 > low_cutoff) & (a_0 < high_cutoff)) 
        highIdx, = np.where(a_0 >= high_cutoff)

        a_0[lowIdx] = ((0.1) ** (-3.0/10.0) - (a_0[lowIdx] / 0.968533)) ** (-10.0/3.0)
        a_0[midIdx] = ((0.5) ** (-6.0/5.0) - ((a_0[midIdx] - low_cutoff) / 0.129758)) ** (-5.0/6.0)
        a_0[highIdx] = (1 - ((a_0[highIdx] - high_cutoff) / 0.0915941)) ** (-10.0/17.0)
        
        return a_0

    def sample_massRatio(self, size=None):
        '''
        Secondary mass is computed from uniform mass ratio distribution draws motivated by
        `Mazeh et al. (1992) <http://adsabs.harvard.edu/abs/1992ApJ...401..265M>`_
        and `Goldberg & Mazeh (1994) <http://adsabs.harvard.edu/abs/1994ApJ...429..362G>`_
        '''
        
        a_0 = np.random.uniform(0.001, 1, size)
        
        return a_0 

    def binary_select(self, primary_mass):
        '''
        Binary fraction is set by `van Haaften et al.(2009)<http://adsabs.harvard.edu/abs/2013A%26A...552A..69V>`_ in appdx
        '''
        binary_fraction = 1/2.0 + 1/4.0 * np.log10(primary_mass)
        binary_choose =  np.random.uniform(0, 1.0, len(binaryfraction)

        binaryIdx, = np.where(binary_fraction > binary_choose)
        singleIdx, = np.where(binary_fraction < binary_choose)

        return primary_mass[binaryIdx], primary_mass[singleIdx]

    
     
