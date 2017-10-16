#! /usr/bin/env python

# Code: Match.py
# Version: 1
# Version changes: PLOT HISTOGRAMS
#
# Input: list of data, # of runs
# Output: list of data that adds cumulatively for each run
#
#
# Edited on:  28 Oct 2015

import numpy as np
import astropy.stats as astroStats



def match(dataCm, nRuns):

    # DEFINE A LIST TO HOLD THE BINNED DATA:   
    histo = [[] for i in range(nRuns)]
    histoBinEdges = [[] for i in range(int(nRuns))]
    # COMPUTE THE BINWIDTH FOR THE MOST COMPLETE DATA SET:
    # NOTE: THIS WILL BE THE BINWIDTH FOR ALL THE HISTOGRAMS IN THE HISTO LIST
    mainHisto, binEdges = astroStats.histogram(np.array(dataCm[len(dataCm)-1]), bins='knuth')
    
    # BIN THE DATA:
    for i in range(nRuns):
        histo[i], histoBinEdges[i] = astroStats.histogram(dataCm[i], bins = binEdges, normed = True)
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
        match[i] = (nominatorSum[i]/np.sqrt(denominator1Sum[i]*denominator2Sum[i]))
        
    return match;
