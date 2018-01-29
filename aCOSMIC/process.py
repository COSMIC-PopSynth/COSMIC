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
# along with aCOSMIC.  If not, see <http://www.gnu.org/licenses/>.

"""`gx_realization`
"""

import numpy as np
import math
import random
import scipy.special as ss
import scipy.stats as stats

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = 'gx_realization'


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


def param_transform(dat):
    '''
    Transforms a data set to limits between zero and one
    '''
        
    datMin = min(dat)-0.000001
    datMax = max(dat)+0.000001
    datZeroed = dat-datMin

    datTransformed = datZeroed/((datMax-datMin)*1.000001)
    return datTransformed

def mass_weighted_number(dat, total_sampled_mass, component_mass): 
    '''
    Compute the total number of systems in the synthetic catalog
    based on the total sampled mass of the simulated system and
    the total mass of a given galactic component
    '''
        
    nSystems = int(len(dat)*component_mass/total_sampled_mass)

    return nSystems

def dat_transform_gw(dat):
    '''
    Transform a data set to have limits between zero and one using 
    param_transform, then transform to logit space
    '''

    m1_trans = param_transform(dat.mass_1)
    m2_trans = param_transform(dat.mass_2)
    porb_trans = param_transform(dat.porb)
    #note: ecc is between zero and one
    ecc_trans = param_transform(dat.ecc)

    dat_trans = ss.logit(np.vstack([m1_trans, m2_trans, porb_trans, ecc_trans]))

    return dat_trans

def dat_un_transform_gw(dat_trans, dat_set):
    dat_min = min(dat_set)
    dat_max = max(dat_set)
    dat = (ss.expit(dat_trans)*(dat_max-dat_min))+dat_min
    
    return dat

def sample(dat, total_sampled_mass, component_mass):
    '''
    Monte-Carlo samples a galactic realization
    '''
    nSystems = mass_weighted_number(dat, total_sampled_mass, component_mass)
    dat_trans = dat_transform_gw(dat)
    
    print nSystems
    if nSystems < 1e7:
        dat_kernel = stats.gaussian_kde(dat_trans)
        gx_sample = dat_kernel.resample(nSystems)
    
        m1 = dat_un_transform_gw(gx_sample[0,:], dat.mass_1)
        m2 = dat_un_transform_gw(gx_sample[1,:], dat.mass_2)
        porb = dat_un_transform_gw(gx_sample[2,:], dat.porb)
        ecc = ss.expit(gx_sample[3,:])

        return np.vstack([m1, m2, porb, ecc])
    else:
        return np.zeros(10)
