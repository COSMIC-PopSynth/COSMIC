# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2019)
#
# This file is part of cosmic.
#
# cosmic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cosmic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cosmic.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy.interpolate import interp1d

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__all__ = ['lisa_characteristic_noise']

def lisa_characteristic_noise():
    '''Computes LISA characteristic strain sensitivity curve according to `Cornish and Robson 2018 <https://arxiv.org/pdf/1803.01944.pdf>`_ without the Galactic foreground

    Parameters
    ----------
    none

    Returns
    -------
    LISA_hc : interpolation of LISA sensitivity curve
    '''
    freq = np.logspace(-9,1,10000)
    # note: freq [Hz], L_arm [m], S_n [Hz^-0.5]
    L_arm = 2.5e9
    f_star = 19.09*1e-3

    P_oms = (1.5e-11)**2*(1. + (2.0e-3/freq)**4)
    P_acc = (3.0e-15)**2*(1. + (0.4e-3/freq)**2)*(1. + (freq/(8.0e-3))**4)

    P_n = (P_oms + 2.*(1. + np.cos(freq/f_star)**2)*P_acc/(2.*np.pi*freq)**4)/L_arm**2
    R = 3./10./(1. + 6./10.*(freq/f_star)**2)
    S_n = (P_n/R*freq)**0.5

    LISA_hc = interp1d(freq, S_n)
    return LISA_hc

