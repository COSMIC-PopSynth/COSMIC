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

def match(dataCm, nRuns):
    """Performs the Match calculation in Eq. 1 of Breivik & Larson (2018)

        Parameters
        ----------
        dataCm : list 
            List of cumulative data for a single paramter

        nRuns : int
            Length of the list

        Returns
        -------
        match : list
            List of matches for each cumulative data set
    
        binwidth : float
            Binwidth of histograms used for match computation
    """

    # DEFINE A LIST TO HOLD THE BINNED DATA:   
    histo = [[] for i in range(nRuns)]
    histoBinEdges = [[] for i in range(int(nRuns))]
    # COMPUTE THE BINWIDTH FOR THE MOST COMPLETE DATA SET:
    # NOTE: THIS WILL BE THE BINWIDTH FOR ALL THE HISTOGRAMS IN THE HISTO LIST
    mainHisto, binEdges = astroStats.histogram(np.array(dataCm[len(dataCm)-1]), bins='scott')
    
    # BIN THE DATA:
    for i in range(nRuns):
        histo[i], histoBinEdges[i] = astroStats.histogram(dataCm[i], bins = binEdges, density = True)
    # COMPUTE THE MATCH:
    nominator = []
    denominator1 = []
    denominator2 = []
    nominatorSum = []
    denominator1Sum = []
    denominator2Sum = []
    match = np.zeros(nRuns-1)

    for i in range(1,nRuns):
        histo2 = histo[i]
        histo1 = histo[i-1]

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


    for i in range(nRuns-1):
        if binEdges[1]-binEdges[0] < 1e-7:
            match[i] = 1.0
        else:
            match[i] = (nominatorSum[i]/np.sqrt(denominator1Sum[i]*denominator2Sum[i]))
        
    return match, binEdges[1]-binEdges[0];
