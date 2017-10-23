# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
#
# This file is part of aCOSMIC.
#
# aCOSMIC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# aCOSMIC is distributed in the hope that it will be useful,
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
    def __init__(self, B_0, bacc, bkick, epoch, massc, ospin, tacc, tphys, size):
        '''
        initialize samples
        '''
        self.B_0 = np.asarray(B_0).repeat(size).reshape(size, 2)
        self.bacc = np.asarray(bacc).repeat(size).reshape(size, 2)
        self.bkick = np.asarray(bkick).repeat(size).reshape(size, 12)
        self.epoch = np.asarray(epoch).repeat(size).reshape(size, 2)
        self.massc = np.asarray(massc).repeat(size).reshape(size, 2)
        self.ospin = np.asarray(ospin).repeat(size).reshape(size, 2)
        self.tacc = np.asarray(tacc).repeat(size).reshape(size, 2)
        self.tphys = np.asarray(tphys).repeat(size)


    # sample primary masses
    def sample_kroupa93(self, size=None):
        '''
        Primary mass follows Kroupa (1993), normalization comes from
        `Hurley (2002) <https://arxiv.org/abs/astro-ph/0201220>`_ between 0.1 and 100 Msun
        '''
        a_0 = np.random.uniform(0.0, 0.9999797, size)
        low_cutoff = 0.740074
        high_cutoff=0.908422


        lowIdx, = np.where(a_0 <= low_cutoff)
        midIdx, = np.where((a_0 > low_cutoff) & (a_0 < high_cutoff)) 
        highIdx, = np.where(a_0 >= high_cutoff)

        a_0[lowIdx] = ((0.1) ** (-3.0/10.0) - (a_0[lowIdx] / 0.968533)) ** (-10.0/3.0)
        a_0[midIdx] = ((0.5) ** (-6.0/5.0) - ((a_0[midIdx] - low_cutoff) / 0.129758)) ** (-5.0/6.0)
        a_0[highIdx] = (1 - ((a_0[highIdx] - high_cutoff) / 0.0915941)) ** (-10.0/17.0)
        
        return a_0


    # sample secondary mass
    def sample_secondary(self, primary_mass):
        '''
        Secondary mass is computed from uniform mass ratio distribution draws motivated by
        `Mazeh et al. (1992) <http://adsabs.harvard.edu/abs/1992ApJ...401..265M>`_
        and `Goldberg & Mazeh (1994) <http://adsabs.harvard.edu/abs/1994ApJ...429..362G>`_
        '''
        
        a_0 = np.random.uniform(0.001, 1, primary_mass.size)
        secondary_mass = primary_mass*a_0        

        return secondary_mass


    def binary_select(self, primary_mass):
        '''
        Binary fraction is set by `van Haaften et al.(2009)<http://adsabs.harvard.edu/abs/2013A%26A...552A..69V>`_ in appdx
        '''
        binary_fraction = 1/2.0 + 1/4.0 * np.log10(primary_mass)
        binary_choose =  np.random.uniform(0, 1.0, binary_fraction.size)

        binaryIdx, = np.where(binary_fraction > binary_choose)
        singleIdx, = np.where(binary_fraction < binary_choose)

        return primary_mass[binaryIdx], primary_mass[singleIdx]


    def sample_separation(self, size=None):
        '''
        Separation is sampled according to `Han (1998)<http://adsabs.harvard.edu/abs/1998MNRAS.296.1019H>`_
        '''

        a_0 = np.random.uniform(0, 1, size)
        low_cutoff = 0.0583333
        lowIdx, = np.where(a_0 <= low_cutoff)
        hiIdx, = np.where(a_0 > low_cutoff)
        
        a_0[lowIdx] = (a_0[lowIdx]/0.00368058)**(5/6.0)
        a_0[hiIdx] = np.exp(a_0[hiIdx]/0.07+math.log(10.0))
        
        # convert to meters
        a_0 = a_0*Rsun
        
        return a_0


    def sep_to_porb(mass1, mass2, sep):
        '''
        Use KEPLER III to convert from separation in meters to orbital period in seconds
        with masses given in solar masses
        '''
        porb_sec = (4*np.pi**2.0/(G*(mass1+mass2)*Msun)*(sep**3.0))**0.5
     
        return porb_sec


    def sample_ecc(self, size=None):
        '''
        Thermal eccentricity distribution following `Heggie (1975)<http://adsabs.harvard.edu/abs/1975MNRAS.173..729H>`_ 
        '''
        a_0 = np.random.uniform(0.0, 1.0, size)
 
        return a_0**0.5


    def sample_constant_SFH(self, t_component, size=None):
        '''
        Assign an evolution time assuming a constant star formation rate over the age: t_component which is specified in [Myr]
        '''

        tphys = np.uniform.random(0, t_component, size)
        return tphys

     
    def sample_gyr_burst(self, t_component, size=None):
        '''
        Assign an evolution time assuming constant star formation rate for 1Gyr starting at 't_component' Myr in the past
        '''

        tphys = t_component - np.uniform.random(0,1000, size)
        return tphys

     
    def set_kstar(self, mass):
        '''
        Initialize all stars according to: kstar=1 if M>=0.7 Msun; kstar=0 if M<0.7
        '''
         
        kstar = np.zeros(size)
        low_cutoff = 0.7
        lowIdx, = np.where(mass < low_cutoff)
        hiIdx, = np.where(mass >= low_cutoff)

        kstar[lowIdx] = 0
        kstar[hiIdx] = 1

        return kstar
