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
import astropy.stats as astroStats

__author__ = 'Katelyn Breivik <katie.breivik@gmail.com>'
__credits__ = 'Scott Coughlin <scott.coughlin@ligo.org>'
__all__ = []


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
    mainHisto, binEdges = astroStats.histogram(np.array(dataCm[0]), bins='knuth')
    
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
        match = 1.0
    else:
        match = np.log10(1-(nominatorSum/np.sqrt(denominator1Sum*denominator2Sum)))

        
    return match[0], binwidth;
