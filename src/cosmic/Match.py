#  -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2021)
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
import warnings
from cosmic.utils import dat_transform


__author__ = "Katelyn Breivik <katie.breivik@gmail.com>"
__credits__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = ["match", "perform_convergence"]


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
        warnings.filterwarnings(
            "ignore", message="divide by zero encountered in double_scalars"
        )
        try:
            bw, binEdges = astroStats.knuth_bin_width(
                np.array(dataCm[0]), return_bins=True
            )
        except Exception:
            bw, binEdges = astroStats.scott_bin_width(
                np.array(dataCm[0]), return_bins=True
            )
    if bw < 1e-4:
        bw = 1e-4
        binEdges = np.arange(binEdges[0], binEdges[-1], bw)

    # BIN THE DATA:
    for i in range(2):
        histo[i], histoBinEdges[i] = astroStats.histogram(
            dataCm[i], bins=binEdges, density=True
        )
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
        nominator.append(histo1[j] * histo2[j])
        denominator1.append((histo1[j] * histo1[j]))
        denominator2.append((histo2[j] * histo2[j]))
    nominatorSum.append(np.sum(nominator))
    denominator1Sum.append(np.sum(denominator1))
    denominator2Sum.append(np.sum(denominator2))

    nominatorSum = np.array(nominatorSum)
    denominator1Sum = np.array(denominator1Sum)
    denominator2Sum = np.array(denominator2Sum)

    binwidth = binEdges[1] - binEdges[0]
    if binwidth < 1e-7:
        match = 1e-9
    else:
        match = np.log10(1 - nominatorSum / np.sqrt(denominator1Sum * denominator2Sum))

    return match[0], binwidth


def perform_convergence(conv_params, conv_1, conv_2, log_file):
    """Performs the convergence calculations for each convergence parameter
    and binary state

    Parameters
    ----------
    conv_params : dict
        List of user supplied convergence parameters
    conv_1 : DataFrame
        Cumulative data set of conv arrays
    conv_2 : DataFrame
        Most recent data set of conv array from most recent BSE run
    log_file : file write
        File to log matches or if the convergence params are not appropriate
        e.g. eccentricity for a disrupted system

    Returns
    -------
    match : array
        Matches for each cumulative data set
    """
    match_all = []
    for conv_param in conv_params:
        if (conv_param == "ecc") and (np.all(conv_1[conv_param] < 1e-7)):
            log_file.write(
                "{0} is circular or disrupted for all conv binaries\n".format(
                    conv_param
                )
            )
            match_all.append(-9)
        elif (conv_param == "ecc") and (np.all(conv_1[conv_param] == -1.0)):
            log_file.write(
                "{0} is the same for all disrupted binaries\n".format(conv_param)
            )
            match_all.append(-9)
        elif conv_param == "ecc":
            conv_1_ecc = conv_1.loc[conv_1.ecc > 0]
            if len(conv_1_ecc) < 10:
                log_file.write("not enough eccentric binaries to compute match")
                match_all.append(-9)
            else:
                conv_2_ecc = conv_2.loc[conv_2.ecc > 0]
                if len(conv_2_ecc) == len(conv_1_ecc):
                    conv_2_ecc = conv_2_ecc[: int(len(conv_2_ecc) / 2)]
                match_compute, bw = match(
                    [
                        dat_transform(conv_1_ecc, [conv_param])[0].tolist(),
                        dat_transform(conv_2_ecc, [conv_param])[0].tolist(),
                    ]
                )
                match_all.append(match_compute)

        elif (conv_param == "porb") and (
            (np.all(conv_1[conv_param] == 0.0)) or (np.all(conv_1[conv_param] == -1.0))
        ):
            log_file.write(
                "{0} is the same for all converging binaries\n".format(conv_param)
            )
            match_all.append(-9)
        elif (conv_param == "sep") and (
            (np.all(conv_1[conv_param] == 0.0)) or (np.all(conv_1[conv_param] == -1.0))
        ):
            log_file.write(
                "{0} is the same for all converging binaries\n".format(conv_param)
            )
            match_all.append(-9)
        else:
            if len(conv_2) == len(conv_1):
                conv_2 = conv_2[: int(len(conv_2) / 2)]
            match_compute, bw = match(
                [
                    dat_transform(conv_1, [conv_param])[0].tolist(),
                    dat_transform(conv_2, [conv_param])[0].tolist(),
                ]
            )
            match_all.append(match_compute)

    log_file.write("matches for converging population are are: {0}\n".format(match_all))
    log_file.write(
        "Number of binaries in converging population is: {0}\n".format(len(conv_1))
    )
    log_file.write("Binwidth is: {0}\n".format(bw))

    return match_all
