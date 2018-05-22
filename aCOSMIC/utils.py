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

'''`utils`
'''
import scipy.integrate
import numpy as np
import scipy.special as ss
import astropy.stats as astrostats

##################################################################################
# DEFINE MIN AND MAX MASS SELECTOR
##################################################################################
def mass_min_max_select(kstar_1, kstar_2):
    '''
    Select a minimum and maximum mass to filter out binaries in the initial
    parameter sample to reduce the number of unneccessary binaries evolved
    in BSE
    '''

    primary_max = 150.0
    secondary_max = 150.0

    primary_min = 0.08
    secondary_min = 0.08
    
    min_mass = [primary_min, secondary_min]
    max_mass = [primary_max, secondary_max]

    if len(kstar_1) == 1:
        # there is a range of final kstar_1s to save
        kstar_1_lo = kstar_1[0]
        kstar_1_hi = kstar_1[0]
    else:
        kstar_1_lo = kstar_1[0]
        kstar_1_hi = kstar_1[1]

    if len(kstar_2) == 1:
        # there is a range of final kstar_1s to save
        kstar_2_lo = kstar_2[0]
        kstar_2_hi = kstar_2[0]
    else:
        kstar_2_lo = kstar_2[0]
        kstar_2_hi = kstar_2[1]

    kstar_lo = [kstar_1_lo, kstar_2_lo]
    kstar_hi = [kstar_1_hi, kstar_2_hi]

    ii = 0
    for k in kstar_lo:
        if k == 14.0:
            min_mass[ii] = 15.0
        elif k == 13.0:
            min_mass[ii] = 6.0
        elif k == 12.0:
            min_mass[ii] = 5.0
        elif k == 11.0:
            min_mass[ii] = 2.0
        elif k == 10.0:
            min_mass[ii] = 0.5
        ii += 1

    ii = 0
    for k in kstar_hi:
        if k == 13.0:
            max_mass[ii] = 50.0
        elif k == 12.0:
            max_mass[ii] = 12.0
        elif k== 11.0:
            max_mass[ii] = 8.0
        elif k == 10.0:
            max_mass[ii] = 5.0
        ii += 1

    return min_mass[0], max_mass[0], min_mass[1], max_mass[1]


def idl_tabulate(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in range(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


def rndm(a, b, g, size):
    '''
    Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b
    '''

    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

def param_transform(dat):
    '''
    Transforms a data set to limits between zero and one
    Leaves some wiggle room on the edges of the data set
    '''

    datMin = min(dat)-0.000001
    datMax = max(dat)+0.000001
    datZeroed = dat-datMin

    datTransformed = datZeroed/((datMax-datMin))
    return datTransformed


def dat_transform(dat, dat_list):
    '''
    Transform a data set to have limits between zero and one using 
    param_transform, then transform to logit space
    '''

    dat_trans = []
    for column in dat_list:
        dat_trans.append(param_transform(dat[column]))
    
    dat_trans = ss.logit(np.vstack([dat_trans]))

    return dat_trans

def dat_un_transform(dat_sample, dat_set, dat_list):
    '''
    Un-transform data that was transformed in dat_transform
    '''
    dat = []
    
    dat_exp = ss.expit(dat_sample)
    for ii,column in zip(range(len(dat_list)),dat_list):
        dat_untrans = dat_exp[ii, :]*\
                   (max(dat_set[column]) - min(dat_set[column])) +\
                    min(dat_set[column])
        dat.append(dat_untrans)
    dat = np.vstack(dat)
    if not np.any(['ecc' in x for x in dat_list]):
        dat = np.vstack([dat, np.zeros(len(dat[0]))])
    return dat

def knuth_bw_selector(dat_list):
    bw_list = []
    for dat in dat_list:
        try:
            bw = astrostats.knuth_bin_width(dat)
        except:
            bw = astrostats.scott_bin_width(dat)
        bw_list.append(bw)
        
    print 'binwidths are: ',bw_list
    return np.min(bw_list)
                      
        


