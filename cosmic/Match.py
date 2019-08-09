# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017)
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

"""`Match`
"""

import numpy as np
import pandas as pd
import astropy.stats as astroStats
import warnings
from cosmic.utils import dat_transform, filter_bpp_bcm, bcm_conv_select


__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = ['match', 'perform_convergence']


def match(dataCm):
    """Performs the Match calculation in Eq. 1 of Breivik & Larson (2018)

        Parameters
        ----------
        dataCm : list
            List of two cumulative data sets for a single paramter

        Returns
        -------
        match : list
            List of matches for each cumulative data set

        binwidth : float
            Binwidth of histograms used for match computation
    """

    # DEFINE A LIST TO HOLD THE BINNED DATA:
    histo = [[], []]
    histoBinEdges = [[], []]
    # COMPUTE THE BINWIDTH FOR THE MOST COMPLETE DATA SET:
    # NOTE: THIS WILL BE THE BINWIDTH FOR ALL THE HISTOGRAMS IN THE HISTO LIST
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in double_scalars")
        try:
            bw, binEdges = astroStats.knuth_bin_width(np.array(dataCm[0]), return_bins=True)
        except:
            bw, binEdges = astroStats.scott_bin_width(np.array(dataCm[0]), return_bins=True)
    if bw < 1e-4:
        bw = 1e-4
        binEdges = np.arange(binEdges[0], binEdges[-1], bw)

    # BIN THE DATA:
    for i in range(2):
        histo[i], histoBinEdges[i] = astroStats.histogram(dataCm[i], bins = binEdges, density = True)
    # COMPUTE THE MATCH:
    nominator = []
    denominator1 = []
    denominator2 = []
    nominatorSum = []
    denominator1Sum = []
    denominator2Sum = []

    histo2 = histo[1]
    histo1 = histo[0]

    for j in range(len(histo1)):
        nominator.append(histo1[j]*histo2[j])
        denominator1.append((histo1[j]*histo1[j]))
        denominator2.append((histo2[j]*histo2[j]))
    nominatorSum.append(np.sum(nominator))
    denominator1Sum.append(np.sum(denominator1))
    denominator2Sum.append(np.sum(denominator2))

    nominatorSum = np.array(nominatorSum, dtype=np.float128)
    denominator1Sum = np.array(denominator1Sum, dtype=np.float128)
    denominator2Sum = np.array(denominator2Sum, dtype=np.float128)

    binwidth = binEdges[1]-binEdges[0]
    if binwidth < 1e-7:
        match = 1e-9
    else:
        match = np.log10(1-nominatorSum/np.sqrt(denominator1Sum*denominator2Sum))


    return match[0], binwidth;

def perform_convergence(conv_params, bin_states, conv_filter,\
                        bcm_save_tot, bcm_save_last,\
                        bpp_save_tot, final_kstar_1, final_kstar_2, log_file):
    """Performs the convergence calculations for each convergence parameter
       and binary state

       Parameters
       ----------
       conv_params : dict
           List of user supplied convergence parameters
       bin_states : dict
           List of user supplied binary states
       conv_filter : dict
           List of user supplied convergence filters
       bcm_save_tot : DataFrame
           Cumulative data set of bcm arrays
       bcm_save_last : DataFrame
           Most recent data set of bcm array from most recent BSE run
       bpp_save_tot : DataFrame
           Cumulative data set of bpp arrays
       log_file : file write
           File to log matches or if the convergence params are not appropriate
           e.g. eccentricity for a disrupted system

       Returns
       -------
       match : array
           Matches for each cumulative data set
    """

    match_lists = []
    for bin_state in bin_states:
        bcm_save_tot_conv = bcm_save_tot.loc[bcm_save_tot.bin_state==bin_state]
        bcm_save_last_conv = bcm_save_last.loc[bcm_save_last.bin_state==bin_state]
        bpp_save_tot_conv = bpp_save_tot.loc[bpp_save_tot.bin_num.isin(bcm_save_tot_conv.bin_num)]
        if bin_state == 1:
            conv_filter = {'formation' : True}
        elif bin_state == 2:
            conv_filter = {'disruption' : True}    
        else:
            # bin_state == 0
            if (conv_filter['formation'] & conv_filter['1_SN']) or
               (conv_filter['formation'] & conv_filter['2_SN']) or
               (conv_filter['formation'] & conv_filter['disruption']) or
               (conv_filter['formation'] & conv_filter['final_state']) or
               (conv_filter['formation'] & conv_filter['XRB_form']) or
               (conv_filter['1_SN'] & conv_filter['2_SN']) or
               (conv_filter['1_SN'] & conv_filter['final_state']) or
               (conv_filter['1_SN'] & conv_filter['XRB_form']) or
               (conv_filter['2_SN'] & conv_filter['final_state']) or
               (conv_filter['2_SN'] & conv_filter['XRB_form']) or
               (conv_filter['final_state'] & conv_filter['XRB_form']) :
                raise ValueError('You specified bin_state == 0 and set more than one convergence filter to True, please set only one formation filter to True')

        conv_1, conv_2 = bcm_conv_select(bcm_save_tot_conv, bcm_save_last_conv,\
                                         bpp_save_tot_conv, final_kstar_1, \
                                         final_kstar_2, conv_filter)
        
        # Perform the Match calculations for all interested parameters
        # supplied by user in conv_params
        if len(conv_2) > 3:
            match_all = []
            for conv_param in conv_params:
                if (conv_param == 'ecc') and (np.all(conv_1[conv_param] < 1e-7)) and (np.all(conv_1[conv_param] >=0.0) and (bin_state == 0):
                    log_file.write('{0} is circular for all binstate {1} binaries\n'.format(conv_param, \
                                                                                            bin_state))
                    match_all.append(-9)
                elif (conv_param == 'ecc') and (bin_state >= 1):
                    log_file.write('{0} is the same for all binstate {1} binaries\n'.format(conv_param, \
                                                                                            bin_state))
                    match_all.append(-9)
                elif (conv_param == 'porb') and (bin_state >= 1):
                    log_file.write('{0} is the same for all binstate {1} binaries\n'.format(conv_param, \
                                                                                            bin_state))
                    match_all.append(-9)
                elif (conv_param == 'sep') and (bin_state >= 1):
                    log_file.write('{0} is the same for all binstate {1} binaries\n'.format(conv_param, \
                                                                                            bin_state))
                    match_all.append(-9)

                else:
                    match_compute, bw = match([dat_transform(bcm_conv_1, [conv_param])[0].tolist(),\
                                               dat_transform(bcm_conv_2, [conv_param])[0].tolist()])
                    match_all.append(match_compute)

            log_file.write('matches for bin state {0} are: {1}\n'.format(bin_state, match_all))
            log_file.write('Number of binaries is: {0}\n'.format(len(bcm_save_conv)))
            log_file.write('Binwidth is: {0}\n'.format(bw))
            log_file.write('\n')
            match_lists.extend(match_all)

        else:
            log_file.write('The filtered bcm array for bin state: {0} does not have >3 values in it yet\n'.format(bin_state))
            log_file.write('Consider larger Nstep sizes\n')
            log_file.write('\n')

    if len(match_lists) > 1:
        match_save = np.array(match_lists)
    else:
        match_save = []

    return match_save
