# -*- coding: utf-8 -*-
# Copyright (C) Scott Coughlin (2017 - 2021)
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

import numpy
from . import zdata

__author__ = "Scott Coughlin <scott.coughlin@ligo.org>"
__all__ = ["zcnsts"]

c = numpy.array([3.040581e-01, 8.049509e-02, 8.967485e-02, 8.780198e-02, 2.219170e-02])


def zcnsts(z):
    """Calculate constants based on the metallicty
    *      zpars:  1; M below which hook doesn't appear on MS, Mhook.
    *              2; M above which He ignition occurs non-degenerately, Mhef.
    *              3; M above which He ignition occurs on the HG, Mfgb.
    *              4; M below which C/O ignition doesn't occur, Mup.
    *              5; M above which C ignites in the centre, Mec.
    *              6; value of log D for M<= zpars[2]
    *              7; value of x for Rgb propto M^(-x)
    *              8; value of x for tMS = numpy.maximum(tHOOK,x*tBGB)
    *              9; constant for McHeIf when computing Mc,BGB, mchefl.
    *             10; constant for McHeIf when computing Mc,HeI, mchefl.
    *             11; hydrogen abundance.
    *             12; helium abundance.
    *             13; constant x in rnumpy.minimum = rgb*x**y used by LM CHeB.
    *             14; z**0.4 to be used for WD L formula.
    """

    # initialize arrays
    msp = numpy.zeros(shape=(200, len(z)))
    zpars = numpy.zeros(shape=(20, len(z)))

    lzs = numpy.log10(z / 0.020)
    dlzs = 1.0 / (z * numpy.log(10.0))
    lz = numpy.log10(z)

    zpars[0] = 1.01850 + lzs * (0.16015 + lzs * 0.0892)
    zpars[1] = 1.9950 + lzs * (0.25 + lzs * 0.087)
    zpars[2] = 16.50 * z ** 0.06 / (1.0 + (1.0e-04 / z) ** 1.27)
    zpars[3] = numpy.maximum(6.110440 + 1.02167 * lzs, 5.0)
    zpars[4] = zpars[3] + 1.80
    zpars[5] = 5.370 + lzs * 0.135
    zpars[6] = c[0] + lzs * (c[1] + lzs * (c[2] + lzs * (c[3] + lzs * c[4])))
    zpars[7] = numpy.maximum(
        0.950,
        numpy.maximum(
            0.95 - (10.0 / 3.0) * (z - 0.01),
            numpy.minimum(0.990, 0.98 - (100.0 / 7.0) * (z - 0.001)),
        ),
    )

    # Lzams

    msp[0] = zdata.xz[0] + lzs * (
        zdata.xz[1] + lzs * (zdata.xz[2] + lzs * (zdata.xz[3] + lzs * zdata.xz[4]))
    )
    msp[1] = zdata.xz[5] + lzs * (
        zdata.xz[6] + lzs * (zdata.xz[7] + lzs * (zdata.xz[8] + lzs * zdata.xz[9]))
    )
    msp[2] = zdata.xz[10] + lzs * (
        zdata.xz[11] + lzs * (zdata.xz[12] + lzs * (zdata.xz[13] + lzs * zdata.xz[14]))
    )
    msp[3] = zdata.xz[15] + lzs * (
        zdata.xz[16] + lzs * (zdata.xz[17] + lzs * (zdata.xz[18] + lzs * zdata.xz[19]))
    )
    msp[4] = zdata.xz[20] + lzs * (
        zdata.xz[21] + lzs * (zdata.xz[22] + lzs * (zdata.xz[23] + lzs * zdata.xz[24]))
    )
    msp[5] = zdata.xz[25] + lzs * (
        zdata.xz[26] + lzs * (zdata.xz[27] + lzs * (zdata.xz[28] + lzs * zdata.xz[29]))
    )
    msp[6] = zdata.xz[30] + lzs * (
        zdata.xz[31] + lzs * (zdata.xz[32] + lzs * (zdata.xz[33] + lzs * zdata.xz[34]))
    )

    # Rzams

    msp[7] = zdata.xz[35] + lzs * (
        zdata.xz[36] + lzs * (zdata.xz[37] + lzs * (zdata.xz[38] + lzs * zdata.xz[39]))
    )
    msp[8] = zdata.xz[40] + lzs * (
        zdata.xz[41] + lzs * (zdata.xz[42] + lzs * (zdata.xz[43] + lzs * zdata.xz[44]))
    )
    msp[9] = zdata.xz[45] + lzs * (
        zdata.xz[46] + lzs * (zdata.xz[47] + lzs * (zdata.xz[48] + lzs * zdata.xz[49]))
    )
    msp[10] = zdata.xz[50] + lzs * (
        zdata.xz[51] + lzs * (zdata.xz[52] + lzs * (zdata.xz[53] + lzs * zdata.xz[54]))
    )
    msp[11] = zdata.xz[55] + lzs * (
        zdata.xz[56] + lzs * (zdata.xz[57] + lzs * (zdata.xz[58] + lzs * zdata.xz[59]))
    )
    msp[12] = zdata.xz[60]
    msp[13] = zdata.xz[61] + lzs * (
        zdata.xz[62] + lzs * (zdata.xz[63] + lzs * (zdata.xz[64] + lzs * zdata.xz[65]))
    )
    msp[14] = zdata.xz[66] + lzs * (
        zdata.xz[67] + lzs * (zdata.xz[68] + lzs * (zdata.xz[69] + lzs * zdata.xz[70]))
    )
    msp[15] = zdata.xz[71] + lzs * (
        zdata.xz[72] + lzs * (zdata.xz[73] + lzs * (zdata.xz[74] + lzs * zdata.xz[75]))
    )

    # Tbgb

    msp[16] = zdata.xt[0] + lzs * (
        zdata.xt[1] + lzs * (zdata.xt[2] + lzs * zdata.xt[3])
    )
    msp[17] = zdata.xt[4] + lzs * (
        zdata.xt[5] + lzs * (zdata.xt[6] + lzs * zdata.xt[7])
    )
    msp[18] = zdata.xt[8] + lzs * (
        zdata.xt[9] + lzs * (zdata.xt[10] + lzs * zdata.xt[11])
    )
    msp[19] = zdata.xt[12] + lzs * (
        zdata.xt[13] + lzs * (zdata.xt[14] + lzs * zdata.xt[15])
    )
    msp[20] = zdata.xt[16]

    # dTbgb/dz
    msp[116] = dlzs * (
        zdata.xt[1] + lzs * (2.0 * zdata.xt[2] + 3.0 * lzs * zdata.xt[3])
    )
    msp[117] = dlzs * (
        zdata.xt[5] + lzs * (2.0 * zdata.xt[6] + 3.0 * lzs * zdata.xt[7])
    )
    msp[118] = dlzs * (
        zdata.xt[9] + lzs * (2.0 * zdata.xt[10] + 3.0 * lzs * zdata.xt[11])
    )
    msp[119] = dlzs * (
        zdata.xt[13] + lzs * (2.0 * zdata.xt[14] + 3.0 * lzs * zdata.xt[15])
    )

    # Thook
    msp[21] = zdata.xt[17] + lzs * (
        zdata.xt[18] + lzs * (zdata.xt[19] + lzs * zdata.xt[20])
    )
    msp[22] = zdata.xt[21]
    msp[23] = zdata.xt[22] + lzs * (
        zdata.xt[23] + lzs * (zdata.xt[24] + lzs * zdata.xt[25])
    )
    msp[24] = zdata.xt[26] + lzs * (
        zdata.xt[27] + lzs * (zdata.xt[28] + lzs * zdata.xt[29])
    )
    msp[25] = zdata.xt[30]

    # Ltms
    msp[26] = zdata.xl[0] + lzs * (
        zdata.xl[1] + lzs * (zdata.xl[2] + lzs * (zdata.xl[3] + lzs * zdata.xl[4]))
    )
    msp[27] = zdata.xl[5] + lzs * (
        zdata.xl[6] + lzs * (zdata.xl[7] + lzs * (zdata.xl[8] + lzs * zdata.xl[9]))
    )
    msp[28] = zdata.xl[10] + lzs * (
        zdata.xl[11] + lzs * (zdata.xl[12] + lzs * zdata.xl[13])
    )
    msp[29] = zdata.xl[14] + lzs * (
        zdata.xl[15] + lzs * (zdata.xl[16] + lzs * (zdata.xl[17] + lzs * zdata.xl[18]))
    )
    msp[26] = msp[26] * msp[29]
    msp[27] = msp[27] * msp[29]
    msp[30] = zdata.xl[19] + lzs * (
        zdata.xl[20] + lzs * (zdata.xl[21] + lzs * zdata.xl[22])
    )
    msp[31] = zdata.xl[23] + lzs * (
        zdata.xl[24] + lzs * (zdata.xl[25] + lzs * zdata.xl[26])
    )

    # Lalpha
    m2 = 2.0
    msp[32] = zdata.xl[27] + lzs * (
        zdata.xl[28] + lzs * (zdata.xl[29] + lzs * zdata.xl[30])
    )
    msp[33] = zdata.xl[31] + lzs * (
        zdata.xl[32] + lzs * (zdata.xl[33] + lzs * zdata.xl[34])
    )
    msp[34] = zdata.xl[35] + lzs * (
        zdata.xl[36] + lzs * (zdata.xl[37] + lzs * zdata.xl[38])
    )
    msp[35] = zdata.xl[39] + lzs * (
        zdata.xl[40] + lzs * (zdata.xl[41] + lzs * zdata.xl[42])
    )
    msp[36] = numpy.maximum(0.90, 1.1064 + lzs * (0.415 + 0.18 * lzs))
    msp[37] = numpy.maximum(1.0, 1.19 + lzs * (0.377 + 0.176 * lzs))

    if (z > 0.010).any():
        msp[36, z > 0.010] = numpy.minimum(msp[36, z > 0.010], 1.0)
        msp[37, z > 0.010] = numpy.minimum(msp[37, z > 0.010], 1.10)

    msp[38] = numpy.maximum(0.1450, 0.0977 - lzs * (0.231 + 0.0753 * lzs))
    msp[39] = numpy.minimum(0.240 + lzs * (0.18 + 0.595 * lzs), 0.306 + 0.053 * lzs)
    msp[40] = numpy.minimum(0.330 + lzs * (0.132 + 0.218 * lzs), 0.36250 + 0.062 * lzs)
    msp[41] = (msp[32] + msp[33] * m2 ** msp[35]) / (m2 ** 0.40 + msp[34] * m2 ** 1.9)

    # Lbeta
    msp[42] = zdata.xl[43] + lzs * (
        zdata.xl[44] + lzs * (zdata.xl[45] + lzs * (zdata.xl[46] + lzs * zdata.xl[47]))
    )
    msp[43] = zdata.xl[48] + lzs * (
        zdata.xl[49] + lzs * (zdata.xl[50] + lzs * (zdata.xl[51] + lzs * zdata.xl[52]))
    )
    msp[44] = zdata.xl[53] + lzs * (zdata.xl[54] + lzs * zdata.xl[55])
    msp[45] = numpy.minimum(1.40, 1.5135 + 0.3769 * lzs)
    msp[45] = numpy.maximum(0.63550 - 0.4192 * lzs, numpy.maximum(1.25, msp[45]))

    # Lhook
    msp[46] = zdata.xl[56] + lzs * (
        zdata.xl[57] + lzs * (zdata.xl[58] + lzs * zdata.xl[59])
    )
    msp[47] = zdata.xl[60] + lzs * (
        zdata.xl[61] + lzs * (zdata.xl[62] + lzs * zdata.xl[63])
    )
    msp[48] = zdata.xl[64] + lzs * (
        zdata.xl[65] + lzs * (zdata.xl[66] + lzs * zdata.xl[67])
    )
    msp[49] = zdata.xl[68] + lzs * (
        zdata.xl[69] + lzs * (zdata.xl[70] + lzs * zdata.xl[71])
    )
    msp[50] = numpy.minimum(1.40, 1.5135 + 0.3769 * lzs)
    msp[50] = numpy.maximum(0.63550 - 0.4192 * lzs, numpy.maximum(1.25, msp[50]))

    # Rtms
    msp[51] = zdata.xr[0] + lzs * (
        zdata.xr[1] + lzs * (zdata.xr[2] + lzs * (zdata.xr[3] + lzs * zdata.xr[4]))
    )
    msp[52] = zdata.xr[5] + lzs * (
        zdata.xr[6] + lzs * (zdata.xr[7] + lzs * (zdata.xr[8] + lzs * zdata.xr[9]))
    )
    msp[53] = zdata.xr[10] + lzs * (
        zdata.xr[11] + lzs * (zdata.xr[12] + lzs * (zdata.xr[13] + lzs * zdata.xr[14]))
    )
    msp[54] = zdata.xr[15] + lzs * (
        zdata.xr[16] + lzs * (zdata.xr[17] + lzs * zdata.xr[18])
    )
    msp[55] = zdata.xr[19] + lzs * (
        zdata.xr[20] + lzs * (zdata.xr[21] + lzs * zdata.xr[22])
    )
    msp[51] = msp[51] * msp[53]
    msp[52] = msp[52] * msp[53]
    msp[56] = zdata.xr[23]
    msp[57] = zdata.xr[24] + lzs * (
        zdata.xr[25] + lzs * (zdata.xr[26] + lzs * zdata.xr[27])
    )
    msp[58] = zdata.xr[28] + lzs * (
        zdata.xr[29] + lzs * (zdata.xr[30] + lzs * zdata.xr[31])
    )
    msp[59] = zdata.xr[32] + lzs * (
        zdata.xr[33] + lzs * (zdata.xr[34] + lzs * zdata.xr[35])
    )
    msp[60] = zdata.xr[36] + lzs * (
        zdata.xr[37] + lzs * (zdata.xr[38] + lzs * zdata.xr[39])
    )
    #
    msp[61] = numpy.maximum(
        0.0970 - 0.1072 * (lz + 3.0),
        numpy.maximum(0.097, numpy.minimum(0.1461, 0.14610 + 0.1237 * (lz + 2.0))),
    )
    msp[61] = 10.0 ** msp[61]
    m2 = msp[61] + 0.10
    msp[62] = (msp[51] + msp[52] * msp[61] ** msp[54]) / (msp[53] + msp[61] ** msp[55])
    msp[63] = (
        msp[56] * m2 ** 3 + msp[57] * m2 ** msp[60] + msp[58] * m2 ** (msp[60] + 1.50)
    ) / (msp[59] + m2 ** 5)

    # Ralpha
    msp[64] = zdata.xr[40] + lzs * (
        zdata.xr[41] + lzs * (zdata.xr[42] + lzs * zdata.xr[43])
    )
    msp[65] = zdata.xr[44] + lzs * (
        zdata.xr[45] + lzs * (zdata.xr[46] + lzs * zdata.xr[47])
    )
    msp[66] = zdata.xr[48] + lzs * (
        zdata.xr[49] + lzs * (zdata.xr[50] + lzs * zdata.xr[51])
    )
    msp[67] = zdata.xr[52] + lzs * (
        zdata.xr[53] + lzs * (zdata.xr[54] + lzs * zdata.xr[55])
    )
    msp[68] = zdata.xr[56] + lzs * (
        zdata.xr[57] + lzs * (zdata.xr[58] + lzs * (zdata.xr[59] + lzs * zdata.xr[60]))
    )
    msp[69] = numpy.maximum(0.90, numpy.minimum(1.0, 1.116 + 0.166 * lzs))
    msp[70] = numpy.maximum(
        1.4770 + 0.296 * lzs, numpy.minimum(1.6, -0.308 - 1.046 * lzs)
    )
    msp[70] = numpy.maximum(0.80, numpy.minimum(0.8 - 2.0 * lzs, msp[70]))
    msp[71] = zdata.xr[61] + lzs * (zdata.xr[62] + lzs * zdata.xr[63])
    msp[72] = numpy.maximum(0.0650, 0.0843 - lzs * (0.0475 + 0.0352 * lzs))
    msp[73] = 0.07360 + lzs * (0.0749 + 0.04426 * lzs)
    if (z < 0.0040).any():
        msp[73, z < 0.0040] = numpy.minimum(0.055, msp[73, z < 0.0040])

    msp[74] = numpy.maximum(0.0910, numpy.minimum(0.121, 0.136 + 0.0352 * lzs))
    msp[75] = (msp[64] * msp[70] ** msp[66]) / (msp[65] + msp[70] ** msp[67])

    if (msp[69] > msp[70]).any():
        msp[69, msp[69] > msp[70]] = msp[70, msp[69] > msp[70]]
        msp[74, msp[69] > msp[70]] = msp[75, msp[69] > msp[70]]

    # Rbeta
    msp[76] = zdata.xr[64] + lzs * (
        zdata.xr[65] + lzs * (zdata.xr[66] + lzs * zdata.xr[67])
    )
    msp[77] = zdata.xr[68] + lzs * (
        zdata.xr[69] + lzs * (zdata.xr[70] + lzs * zdata.xr[71])
    )
    msp[78] = zdata.xr[72] + lzs * (
        zdata.xr[73] + lzs * (zdata.xr[74] + lzs * zdata.xr[75])
    )
    msp[79] = zdata.xr[76] + lzs * (
        zdata.xr[77] + lzs * (zdata.xr[78] + lzs * zdata.xr[79])
    )
    msp[80] = zdata.xr[80] + lzs * (zdata.xr[81] + lzs * lzs * zdata.xr[82])
    if (z > 0.010).any():
        msp[80, z > 0.010] = numpy.maximum(msp[80, z > 0.010], 0.95)

    msp[81] = numpy.maximum(
        1.40, numpy.minimum(1.6, 1.6 + lzs * (0.764 + 0.3322 * lzs))
    )

    # Rgamma
    msp[82] = numpy.maximum(
        zdata.xr[83] + lzs * (zdata.xr[84] + lzs * (zdata.xr[85] + lzs * zdata.xr[86])),
        zdata.xr[95] + lzs * (zdata.xr[96] + lzs * zdata.xr[97]),
    )
    msp[83] = numpy.minimum(
        0.0,
        zdata.xr[87] + lzs * (zdata.xr[88] + lzs * (zdata.xr[89] + lzs * zdata.xr[90])),
    )
    msp[83] = numpy.maximum(
        msp[83], zdata.xr[98] + lzs * (zdata.xr[99] + lzs * zdata.xr[100])
    )
    msp[84] = zdata.xr[91] + lzs * (
        zdata.xr[92] + lzs * (zdata.xr[93] + lzs * zdata.xr[94])
    )
    msp[84] = numpy.maximum(0.0, numpy.minimum(msp[84], 7.454 + 9.046 * lzs))
    msp[85] = numpy.minimum(
        zdata.xr[101] + lzs * zdata.xr[102], numpy.maximum(2.0, -13.3 - 18.6 * lzs)
    )
    msp[86] = numpy.minimum(1.50, numpy.maximum(0.4, 2.493 + 1.1475 * lzs))
    msp[87] = numpy.maximum(1.0, numpy.minimum(1.27, 0.8109 - 0.6282 * lzs))
    msp[87] = numpy.maximum(msp[87], 0.63550 - 0.4192 * lzs)
    msp[88] = numpy.maximum(5.855420e-02, -0.27110 - lzs * (0.5756 + 0.0838 * lzs))

    # Rhook
    msp[89] = zdata.xr[103] + lzs * (
        zdata.xr[104] + lzs * (zdata.xr[105] + lzs * zdata.xr[106])
    )
    msp[90] = zdata.xr[107] + lzs * (
        zdata.xr[108] + lzs * (zdata.xr[109] + lzs * zdata.xr[110])
    )
    msp[91] = zdata.xr[111] + lzs * (
        zdata.xr[112] + lzs * (zdata.xr[113] + lzs * zdata.xr[114])
    )
    msp[92] = zdata.xr[115] + lzs * (
        zdata.xr[116] + lzs * (zdata.xr[117] + lzs * zdata.xr[118])
    )
    msp[93] = numpy.minimum(
        1.250, numpy.maximum(1.10, 1.9848 + lzs * (1.1386 + 0.3564 * lzs))
    )
    msp[94] = 0.0630 + lzs * (0.0481 + 0.00984 * lzs)
    msp[95] = numpy.minimum(1.30, numpy.maximum(0.45, 1.2 + 2.45 * lzs))

    # Lneta
    msp[96, z > 0.00090] = 10.0
    msp[96, z < 0.00090] = 20.0

    # converting the below from fortran to python still needs work
    """
    # Lbgb
    gbp[0] = zdata.xg[0]+lzs*(zdata.xg[1]+lzs*(zdata.xg[2]+lzs*zdata.xg[3]))
    gbp[1] = zdata.xg[4]+lzs*(zdata.xg[5]+lzs*(zdata.xg[6]+lzs*zdata.xg[7]))
    gbp[2] = zdata.xg[8]+lzs*(zdata.xg[9]+lzs*(zdata.xg[10]+lzs*zdata.xg[11]))
    gbp[3] = zdata.xg[12]+lzs*(zdata.xg[13]+lzs*(zdata.xg[14]+lzs*zdata.xg[15]))
    gbp[4] = zdata.xg[16]+lzs*(zdata.xg[17]+lzs*zdata.xg[18])
    gbp[5] = zdata.xg[19]+lzs*(zdata.xg[20]+lzs*zdata.xg[21])
    gbp[2] = gbp[2]**gbp[5]
    gbp[6] = zdata.xg[22]
    gbp[7] = zdata.xg[23]

    # Lbagb
    # set gbp[15] = 1.0 until it is reset later with an initial
    # call to Lbagbf using mass = zpars[1] and mhefl = 0.0
    gbp[8] = zdata.xg[24] + lzs*(zdata.xg[25] + lzs*zdata.xg[26])
    gbp[9] = zdata.xg[27] + lzs*(zdata.xg[28] + lzs*zdata.xg[29])
    gbp[10] = 15.0
    gbp[11] = zdata.xg[30]+lzs*(zdata.xg[31]+lzs*(zdata.xg[32]+lzs*zdata.xg[33]))
    gbp[12] = zdata.xg[34]+lzs*(zdata.xg[35]+lzs*(zdata.xg[36]+lzs*zdata.xg[37]))
    gbp[13] = zdata.xg[38]+lzs*(zdata.xg[39]+lzs*(zdata.xg[40]+lzs*zdata.xg[41]))
    gbp[14] = zdata.xg[42]+lzs*zdata.xg[43]
    gbp[11] = gbp[11]**gbp[14]
    gbp[13] = gbp[13]**gbp[14]
    gbp[15] = 1.0

    # Rgb
    gbp[16] = -4.67390-0.9394*lz
    gbp[16] = 10.0**gbp[16]
    gbp[16] = numpy.maximum(gbp[16],-0.041670+55.67*z)
    gbp[16] = numpy.minimum(gbp[16],0.47710-9329.21*z**2.94)
    gbp[17] = numpy.minimum(0.540,0.397+lzs*(0.28826+0.5293*lzs))
    gbp[18] = numpy.maximum(-0.14510,-2.2794-lz*(1.5175+0.254*lz))
    gbp[18] = 10.0**gbp[18]
    if (z > 0.0040):
        gbp[18] = numpy.maximum(gbp[18],0.73070+14265.1*z**3.395)

    gbp[19] = zdata.xg[44]+lzs*(zdata.xg[45]+lzs*(zdata.xg[46]+lzs*(zdata.xg[47]+lzs*(zdata.xg[48]+lzs*zdata.xg[49]))))
    gbp[20] = zdata.xg[50]+lzs*(zdata.xg[51]+lzs*(zdata.xg[52]+lzs*(zdata.xg[53]+lzs*zdata.xg[54])))
    gbp[21] = zdata.xg[55]+lzs*(zdata.xg[56]+lzs*(zdata.xg[57]+lzs*(zdata.xg[58]+lzs*(zdata.xg[59]+lzs*zdata.xg[60]))))
    gbp[22] = zdata.xg[61]+lzs*(zdata.xg[62]+lzs*(zdata.xg[63]+lzs*(zdata.xg[64]+lzs*zdata.xg[65])))

    # Ragb
    gbp[23] = numpy.minimum(0.991640-743.123*z**2.83,
                  1.04220+lzs*(0.13156+0.045*lzs))
    gbp[24] = zdata.xg[66]+lzs*(zdata.xg[67]+lzs*(zdata.xg[68]+lzs*(zdata.xg[69]+ lzs*(zdata.xg[70]+lzs*zdata.xg[71]))))
    gbp[25] = zdata.xg[72]+lzs*(zdata.xg[73]+lzs*(zdata.xg[74]+lzs*(zdata.xg[75]+lzs*zdata.xg[76])))
    gbp[26] = zdata.xg[77]+lzs*(zdata.xg[78]+lzs*(zdata.xg[79]+lzs*(zdata.xg[80]+lzs*(zdata.xg[81]+lzs*zdata.xg[82]))))
    gbp[27] = zdata.xg[83]+lzs*(zdata.xg[84]+lzs*(zdata.xg[85]+lzs*(zdata.xg[86]+lzs*zdata.xg[87])))
    gbp[28] = zdata.xg[88]+lzs*(zdata.xg[89]+lzs*(zdata.xg[90]+lzs*(zdata.xg[91]+lzs*(zdata.xg[92]+lzs*zdata.xg[93]))))
    gbp[29] = zdata.xg[94]+lzs*(zdata.xg[95]+lzs*(zdata.xg[96]+lzs*(zdata.xg[97]+lzs*(zdata.xg[98]+lzs*zdata.xg[99]))))
    m1 = zpars[1] - 0.20
    gbp[30] = gbp[28] + gbp[29]*m1
    gbp[31] = numpy.minimum(gbp[24]/zpars[1]**gbp[25],gbp[26]/zpars[1]**gbp[27])

    # Mchei
    gbp[32] = zdata.xg[100]**4
    gbp[33] = zdata.xg[101]*4.0

    # Mcagb
    gbp[34] = zdata.xg[102]+lzs*(zdata.xg[103]+lzs*(zdata.xg[104]+lzs*zdata.xg[105]))
    gbp[35] = zdata.xg[106]+lzs*(zdata.xg[107]+lzs*(zdata.xg[108]+lzs*zdata.xg[109]))
    gbp[36] = zdata.xg[110]+lzs*zdata.xg[111]
    gbp[34] = gbp[34]**4
    gbp[35] = gbp[35]*4.0
    gbp[36] = gbp[36]**4

    # Lhei
    # set gbp[40] = -1.0 until it is reset later with an initial
    # call to Lheif using mass = zpars[1] and mhefl = 0.0
    gbp[37] = zdata.xh[0]+lzs*zdata.xh[1]
    gbp[38] = zdata.xh[2]+lzs*zdata.xh[3]
    gbp[39] = zdata.xh[4]
    gbp[40] = -1.0
    gbp[41] = zdata.xh[5]+lzs*(zdata.xh[6]+lzs*zdata.xh[7])
    gbp[42] = zdata.xh[8]+lzs*(zdata.xh[9]+lzs*zdata.xh[10])
    gbp[43] = zdata.xh[11]+lzs*(zdata.xh[12]+lzs*zdata.xh[13])
    gbp[41] = gbp[41]**2
    gbp[43] = gbp[43]**2
    # Lhe
    gbp[44] = zdata.xh[14]+lzs*(zdata.xh[15]+lzs*zdata.xh[16])
    if (lzs > -1.0):
        gbp[45] = 1.0 - zdata.xh[18]*(lzs+1.0)**zdata.xh[17]
    else:
        gbp[45] = 1.0

    gbp[46] = zdata.xh[19]+lzs*(zdata.xh[20]+lzs*zdata.xh[21])
    gbp[47] = zdata.xh[22]+lzs*(zdata.xh[23]+lzs*zdata.xh[24])
    gbp[44] = gbp[44]**gbp[47]
    gbp[46] = gbp[46]**gbp[47]
    gbp[45] = gbp[45]/zpars[2]**0.10+(gbp[45]*gbp[46]-gbp[44])/zpars[2]**(gbp[47]+0.10)

    # Rnumpy.minimum

    gbp[48] = zdata.xh[25]+lzs*(zdata.xh[26]+lzs*(zdata.xh[27]+lzs*zdata.xh[28]))
    gbp[49] = zdata.xh[29]+lzs*(zdata.xh[30]+lzs*(zdata.xh[31]+lzs*zdata.xh[32]))
    gbp[50] = zdata.xh[33]+lzs*(zdata.xh[34]+lzs*(zdata.xh[35]+lzs*zdata.xh[36]))
    gbp[51] = 5.0+zdata.xh[37]*z**zdata.xh[38]
    gbp[52] = zdata.xh[39]+lzs*(zdata.xh[40]+lzs*(zdata.xh[41]+lzs*zdata.xh[42]))
    gbp[48] = gbp[48]**gbp[52]
    gbp[50] = gbp[50]**(2.0*gbp[52])

    # The
    # set gbp[56] = -1.0 until it is reset later with an initial
    # call to Thef using mass = zpars[1], mc = 0.0  and mhefl = 0.0
    gbp[53] = zdata.xh[43]+lzs*(zdata.xh[44]+lzs*(zdata.xh[45]+lzs*zdata.xh[46]))
    gbp[54] = zdata.xh[47]+lzs*(zdata.xh[48]+lzs*zdata.xh[49])
    gbp[54] = numpy.maximum(gbp[54],1.0)
    gbp[55] = zdata.xh[50]
    gbp[56] = -1.0
    gbp[57] = zdata.xh[51]+lzs*(zdata.xh[52]+lzs*(zdata.xh[53]+lzs*zdata.xh[54]))
    gbp[58] = zdata.xh[55]+lzs*(zdata.xh[56]+lzs*(zdata.xh[57]+lzs*zdata.xh[58]))
    gbp[59] = zdata.xh[59]+lzs*(zdata.xh[60]+lzs*(zdata.xh[61]+lzs*zdata.xh[62]))
    gbp[60] = zdata.xh[63]+lzs*zdata.xh[64]
    gbp[57] = gbp[57]**gbp[60]
    gbp[59] = gbp[59]**5

    # Tbl
    dum1 = zpars[1]/zpars[2]
    gbp[61] = zdata.xh[65]+lzs*zdata.xh[66]
    gbp[61] = -gbp[61]*numpy.log10(dum1)
    gbp[62] = zdata.xh[67]
    if (lzd > 0.0):
        gbp[63] = 1.0-lzd*(zdata.xh[68]+lzd*(zdata.xh[69]+lzd*zdata.xh[70]))
    else:
        gbp[63] = 1.0

    gbp[64] = 1.0-gbp[63]*dum1**gbp[62]
    gbp[65] = 1.0 - lzd*(zdata.xh[76] + lzd*(zdata.xh[77] + lzd*zdata.xh[78]))
    gbp[66] = zdata.xh[71] + lzs*(zdata.xh[72] + lzs*(zdata.xh[73] + lzs*zdata.xh[74]))
    gbp[67] = zdata.xh[75]

    # Lzahb
    gbp[68] = zdata.xh[79] + lzs*(zdata.xh[80] + lzs*zdata.xh[81])
    gbp[69] = zdata.xh[82] + lzs*(zdata.xh[83] + lzs*zdata.xh[84])
    gbp[70] = 15.0
    gbp[71] = zdata.xh[85]
    gbp[72] = zdata.xh[86]

    # Rzahb
    gbp[74] = zdata.xh[87] + lzs*(zdata.xh[88] + lzs*(zdata.xh[89] + lzs*zdata.xh[90]))
    gbp[75] = zdata.xh[91] + lzs*(zdata.xh[92] + lzs*(zdata.xh[93] + lzs*zdata.xh[94]))
    gbp[76] = zdata.xh[95] + lzs*(zdata.xh[96] + lzs*(zdata.xh[97] + lzs*zdata.xh[98]))

    # finish Lbagb
    mhefl = 0.0
    lx = lbagbf(zpars[1], mhefl)
    gbp[15] = lx

    # finish LHeI
    dum1 = 0.0
    lhefl = lheif(zpars[1],mhefl)
    gbp[40] = (gbp[37]*zpars[1]**gbp[38]-lhefl)/(numpy.exp(zpars[1]*gbp[39])*lhefl)

    # finish THe
    thefl = thef(zpars[1],dum1,mhefl)*tbgbf(zpars[1])
    gbp[56] = (thefl-gbp[53])/(gbp[53]*numpy.exp(gbp[55]*zpars[1]))

    # finish Tblf
    rb = ragbf(zpars[2],lheif(zpars[2],zpars[1]),mhefl)
    rr = 1.0 - rnumpy.minimumf(zpars[2])/rb
    rr = numpy.maximum(rr,1.0e-12)
    gbp[65] = gbp[65]/(zpars[2]**gbp[66]*rr**gbp[67])

    # finish Lzahb
    gbp[73] = lhefl*lHef(zpars[1])

    kw = 0
    tm = 0.0
    tn = 0.0
    star(kw,zpars[1],zpars[1],tm,tn,tscls,lums,GB,zpars)
    zpars[8] = mcgbf(lums[2],GB,lums[5])
    zpars[9] = mcgbf(lums[3],GB,lums[5])

    # set the hydrogen and helium abundances
    zpars[10] = 0.760 - 3.0*z
    zpars[11] = 0.240 + 2.0*z

    # set constant for low-mass CHeB stars
    zpars[12] = (rnumpy.minimumf(zpars[1])/
                rgbf(zpars[1],lzahbf(zpars[1],zpars[8],zpars[1])))

    zpars[13] = z**0.40
    #
    """
    return zpars, msp
