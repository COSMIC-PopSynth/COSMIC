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


G = 6.67384*math.pow(10, -11.0)
c = 2.99792458*math.pow(10, 8.0)
parsec = 3.08567758*math.pow(10, 16)
Rsun = 6.955*math.pow(10, 8)
Msun = 1.9891*math.pow(10,30)
day = 86400.0
rsun_in_au = 215.0954
day_in_year = 365.242
sec_in_day = 86400.0
sec_in_hour = 3600.0
hrs_in_day = 24.0
sec_in_year = 3.15569*10**7.0
Tobs = 3.15569*10**7.0
geo_mass = G/c**2


class Sample:
    def __init__(self, metallicity, size=None):
        '''
        initialize samples
        '''
        self.metallicity = np.asarray(metallicity).repeat(size)


    # sample primary masses
    def sample_primary(self, kstar1_final, model='kroupa93', size=None):
        '''
        kroupa93 follows Kroupa (1993), normalization comes from
        `Hurley (2002) <https://arxiv.org/abs/astro-ph/0201220>`_ between 0.1 and 100 Msun
        salpter55 follows Salpeter (1955) <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_ between 0.1 and 100 Msun
        '''

        if model=='kroupa93':
            # If the final binary contains a compact object (BH or NS),
            # we want to evolve 'size' binaries that could form a compact
            # object so we over sample the initial population
            if kstar1_final > 14.0:
                a_0 = np.random.uniform(0.0, 0.9999797, size*500)
            elif kstar1_final > 12.0:
                a_0 = np.random.uniform(0.0, 0.9999797, size*50)
            else:   
                a_0 = np.random.uniform(0.0, 0.9999797, size)
            
            low_cutoff = 0.740074
            high_cutoff=0.908422
     
            lowIdx, = np.where(a_0 <= low_cutoff)
            midIdx, = np.where((a_0 > low_cutoff) & (a_0 < high_cutoff)) 
            highIdx, = np.where(a_0 >= high_cutoff)
     
            a_0[lowIdx] = ((0.1) ** (-3.0/10.0) - (a_0[lowIdx] / 0.968533)) ** (-10.0/3.0)
            a_0[midIdx] = ((0.5) ** (-6.0/5.0) - ((a_0[midIdx] - low_cutoff) / 0.129758)) ** (-5.0/6.0)
            a_0[highIdx] = (1 - ((a_0[highIdx] - high_cutoff) / 0.0915941)) ** (-10.0/17.0)
            
            total_sampled_mass = np.sum(a_0)
            if kstar1_final > 13.0:
                a_0 = a_0[a_0 > 15.0]
            elif kstar1_final > 12.0: 
                a_0 = a_0[a_0 > 8.0]
            return a_0, total_sampled_mass
        
        elif model=='salpeter55':
            # If the final binary contains a compact object (BH or NS),
            # we want to evolve 'size' binaries that could form a compact
            # object so we over sample the initial population
            if kstar1_final > 14.0:
                mSamp = np.random.power(-1.3, size*500)
            elif kstar1_final > 12.0:
                mSamp = np.random.power(-1.3, size*50)
            else:
                mSamp = np.random.power(-1.3, size*500)
            
            # Transform sample from 0 to 1 to be between 0.1 and 100 Msun
            a_0 = mSamp*(100.0-0.1)+0.1

            total_sampled_mass = np.sum(a_0)
            if kstar1_final > 13.0:
                a_0 = a_0[a_0 > 15.0]
            elif kstar1_final > 12.0:
                a_0 = a_0[a_0 > 8.0]
            return a_0, total_sampled_mass 
                   
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


    def binary_select(self, primary_mass, model='vanHaaften'):
        '''
        Binary fraction is set by `van Haaften et al.(2009)<http://adsabs.harvard.edu/abs/2013A%26A...552A..69V>`_ in appdx
        '''
        
        if model=='vanHaaften':
            binary_fraction = 1/2.0 + 1/4.0 * np.log10(primary_mass)
            binary_choose =  np.random.uniform(0, 1.0, binary_fraction.size)
     
            binaryIdx, = np.where(binary_fraction > binary_choose)
            singleIdx, = np.where(binary_fraction < binary_choose)
     
            return primary_mass[binaryIdx], primary_mass[singleIdx]

    def sep_to_porb(self, mass1, mass2, sep):
        '''
        Use KEPLER III to convert from separation in meters to orbital period in seconds
        with masses given in solar masses
        '''
        
        porb_sec = (4*np.pi**2.0/(G*(mass1+mass2)*Msun)*(sep**3.0))**0.5

        return porb_sec

    def sample_porb(mass1, mass2, model='Han', size=None):
        '''
        If model='Han', separation is sampled according to `Han (1998)<http://adsabs.harvard.edu/abs/1998MNRAS.296.1019H>`_
        If model='log_normal', Use (ref from Aaron Geller) to sample a log normal distributed orbital period
        with mu=1e5.03 days and sigma=1e2.28 days; converts period from days to seconds.
        Separation is then converted to orbital period in seconds
        '''
        
        if model=='Han': 
            a_0 = np.random.uniform(0, 1, size)
            low_cutoff = 0.0583333
            lowIdx, = np.where(a_0 <= low_cutoff)
            hiIdx, = np.where(a_0 > low_cutoff)
            
            a_0[lowIdx] = (a_0[lowIdx]/0.00368058)**(5/6.0)
            a_0[hiIdx] = np.exp(a_0[hiIdx]/0.07+math.log(10.0))
            
            # convert to meters
            a_0 = a_0*Rsun
            
            porb_sec = sep_to_porb(mass1, mass2, a_0)
            return porb_sec
     
        if model=='log_normal':
            porb = np.random.lognormal(mean=1e5.03, sigma=1e2.28, size=size)

            return porb*sec_in_day


    def sample_ecc(self, model='thermal', size=None):
        '''
        If model=='thermal', thermal eccentricity distribution following `Heggie (1975)<http://adsabs.harvard.edu/abs/1975MNRAS.173..729H>`_ 
        
        If model=='uniform', Sample eccentricities uniformly between 0 and 1 following ref ref from Aaron Geller '''
        
        if model=='thermal':
            a_0 = np.random.uniform(0.0, 1.0, size)
 
            return a_0**0.5

        if model=='uniform':
            ecc = np.random.uniform(0.0, 1.0, size)

            return ecc


    def sample_SFH(self, model='const', size=None):
        '''
        'const': Assign an evolution time assuming a constant star formation rate over the age of the MW disk: 10,000 [Myr]
        'burst': Assign an evolution time assuming constant star formation rate for 1Gyr starting at 't_component' 10,000 [Myr] in the past
        '''

        if model='const':

            tphys = np.random.uniform(0, 10000.0, size)
            return tphys
 
        if model='burst':
            tphys = 10000.0 - np.random.uniform(0, 1000, size)
            return tphys
     
     
    def set_kstar(self, mass):
        '''
        Initialize all stars according to: kstar=1 if M>=0.7 Msun; kstar=0 if M<0.7
        '''
         
        kstar = np.zeros(mass.size)
        low_cutoff = 0.7
        lowIdx, = np.where(mass < low_cutoff)
        hiIdx, = np.where(mass >= low_cutoff)

        kstar[lowIdx] = 0
        kstar[hiIdx] = 1

        return kstar
