import numpy
import zdata

c = numpy.array([3.040581e-01, 8.049509e-02, 8.967485e-02, 8.780198e-02, 2.219170e-02])

def zcnsts(z):
    """
    *
    *       ------------------------------------------------------------
    *
    *      zpars:  1; M below which hook doesn't appear on MS, Mhook.
    *              2; M above which He ignition occurs non-degenerately, Mhef.
    *              3; M above which He ignition occurs on the HG, Mfgb.
    *              4; M below which C/O ignition doesn't occur, Mup.
    *              5; M above which C ignites in the centre, Mec.
    *              6; value of log D for M<= zpars[3]
    *              7; value of x for Rgb propto M^(-x)
    *              8; value of x for tMS = max(tHOOK,x*tBGB)
    *              9; constant for McHeIf when computing Mc,BGB, mchefl.
    *             10; constant for McHeIf when computing Mc,HeI, mchefl.
    *             11; hydrogen abundance.
    *             12; helium abundance.
    *             13; constant x in rmin = rgb*x**y used by LM CHeB.
    *             14; z**0.4 to be used for WD L formula.
    *
    *       ------------------------------------------------------------
    *
    """

    # initialize arrays
    msp = numpy.zeros(200)
    gbp = numpy.zeros(200)
    zpars = numpy.zeros(20)
    tscls = numpy.zeros(20)
    lums = numpy.zeros(10)
    GB = numpy.zeros(10)

    lzs = numpy.log10(z/0.020)
    dlzs = 1.0/(z*numpy.log(10.0))
    lz = numpy.log10(z)
    lzd = lzs + 1.0

    zpars[1] = 1.01850 + lzs*(0.16015 + lzs*0.0892)
    zpars[2] = 1.9950 + lzs*(0.25 + lzs*0.087)
    zpars[3] = 16.50*z**0.06/(1.0 + (1.0e-04/z)**1.27)
    zpars[4] = max(6.110440 + 1.02167*lzs, 5.0)
    zpars[5] = zpars[4] + 1.80
    zpars[6] = 5.370 + lzs*0.135
    zpars[7] = c[1] + lzs*(c[2] + lzs*(c[3] + lzs*(c[4] + lzs*c[5])))
    zpars[8] = max(0.950,max(0.95-(10.0/3.0)*(z-0.01),
               min(0.990,0.98-(100.0/7.0)*(z-0.001))))

    import pdb
    pdb.set_trace()

    # Lzams

    msp[1] = zdata.xz[1]+lzs*(zdata.xz[2]+lzs*(zdata.xz[3]+lzs*(zdata.xz[4]+lzs*zdata.xz[5])))
    msp[2] = zdata.xz[6]+lzs*(zdata.xz[7]+lzs*(zdata.xz[8]+lzs*(zdata.xz[9]+lzs*zdata.xz[10])))
    msp[3] = zdata.xz[11]+lzs*(zdata.xz[12]+lzs*(zdata.xz[13]+lzs*(zdata.xz[14]+lzs*zdata.xz[15])))
    msp[4] = zdata.xz[16]+lzs*(zdata.xz[17]+lzs*(zdata.xz[18]+lzs*(zdata.xz[19]+lzs*zdata.xz[20])))
    msp[5] = zdata.xz[21]+lzs*(zdata.xz[22]+lzs*(zdata.xz[23]+lzs*(zdata.xz[24]+lzs*zdata.xz[25])))
    msp[6] = zdata.xz[26]+lzs*(zdata.xz[27]+lzs*(zdata.xz[28]+lzs*(zdata.xz[29]+lzs*zdata.xz[30])))
    msp[7] = zdata.xz[31]+lzs*(zdata.xz[32]+lzs*(zdata.xz[33]+lzs*(zdata.xz[34]+lzs*zdata.xz[35])))

    # Rzams

    msp[8] = zdata.xz[36]+lzs*(zdata.xz[37]+lzs*(zdata.xz[38]+lzs*(zdata.xz[39]+lzs*zdata.xz[40])))
    msp[9] = zdata.xz[41]+lzs*(zdata.xz[42]+lzs*(zdata.xz[43]+lzs*(zdata.xz[44]+lzs*zdata.xz[45])))
    msp[10] = zdata.xz[46]+lzs*(zdata.xz[47]+lzs*(zdata.xz[48]+lzs*(zdata.xz[49]+lzs*zdata.xz[50])))
    msp[11] = zdata.xz[51]+lzs*(zdata.xz[52]+lzs*(zdata.xz[53]+lzs*(zdata.xz[54]+lzs*zdata.xz[55])))
    msp[12] = zdata.xz[56]+lzs*(zdata.xz[57]+lzs*(zdata.xz[58]+lzs*(zdata.xz[59]+lzs*zdata.xz[60])))
    msp[13] = zdata.xz[61]
    msp[14] = zdata.xz[62]+lzs*(zdata.xz[63]+lzs*(zdata.xz[64]+lzs*(zdata.xz[65]+lzs*zdata.xz[66])))
    msp[15] = zdata.xz[67]+lzs*(zdata.xz[68]+lzs*(zdata.xz[69]+lzs*(zdata.xz[70]+lzs*zdata.xz[71])))
    msp[16] = zdata.xz[72]+lzs*(zdata.xz[73]+lzs*(zdata.xz[74]+lzs*(zdata.xz[75]+lzs*zdata.xz[76])))

    # Tbgb

    msp[17] = zdata.xt[1]+lzs*(zdata.xt[2]+lzs*(zdata.xt[3]+lzs*zdata.xt[4]))
    msp[18] = zdata.xt[5]+lzs*(zdata.xt[6]+lzs*(zdata.xt[7]+lzs*zdata.xt[8]))
    msp[19] = zdata.xt[9]+lzs*(zdata.xt[10]+lzs*(zdata.xt[11]+lzs*zdata.xt[12]))
    msp[20] = zdata.xt[13]+lzs*(zdata.xt[14]+lzs*(zdata.xt[15]+lzs*zdata.xt[16]))
    msp[21] = zdata.xt[17]

    # dTbgb/dz
    msp[117] = dlzs*(zdata.xt[2]+lzs*(2.0*zdata.xt[3]+3.0*lzs*zdata.xt[4]))
    msp[118] = dlzs*(zdata.xt[6]+lzs*(2.0*zdata.xt[7]+3.0*lzs*zdata.xt[8]))
    msp[119] = dlzs*(zdata.xt[10]+lzs*(2.0*zdata.xt[11]+3.0*lzs*zdata.xt[12]))
    msp[120] = dlzs*(zdata.xt[14]+lzs*(2.0*zdata.xt[15]+3.0*lzs*zdata.xt[16]))

    # Thook
    msp[22] = zdata.xt[18]+lzs*(zdata.xt[19]+lzs*(zdata.xt[20]+lzs*zdata.xt[21]))
    msp[23] = zdata.xt[22]
    msp[24] = zdata.xt[23]+lzs*(zdata.xt[24]+lzs*(zdata.xt[25]+lzs*zdata.xt[26]))
    msp[25] = zdata.xt[27]+lzs*(zdata.xt[28]+lzs*(zdata.xt[29]+lzs*zdata.xt[30]))
    msp[26] = zdata.xt[31]

    # Ltms
    msp[27] = zdata.xl[1]+lzs*(zdata.xl[2]+lzs*(zdata.xl[3]+lzs*(zdata.xl[4]+lzs*zdata.xl[5])))
    msp[28] = zdata.xl[6]+lzs*(zdata.xl[7]+lzs*(zdata.xl[8]+lzs*(zdata.xl[9]+lzs*zdata.xl[10])))
    msp[29] = zdata.xl[11]+lzs*(zdata.xl[12]+lzs*(zdata.xl[13]+lzs*zdata.xl[14]))
    msp[30] = zdata.xl[15]+lzs*(zdata.xl[16]+lzs*(zdata.xl[17]+lzs*(zdata.xl[18]+lzs*zdata.xl[19])))
    msp[27] = msp[27]*msp[30]
    msp[28] = msp[28]*msp[30]
    msp[31] = zdata.xl[20]+lzs*(zdata.xl[21]+lzs*(zdata.xl[22]+lzs*zdata.xl[23]))
    msp[32] = zdata.xl[24]+lzs*(zdata.xl[25]+lzs*(zdata.xl[26]+lzs*zdata.xl[27]))

    # Lalpha
    m2 = 2.0
    msp[33] = zdata.xl[28]+lzs*(zdata.xl[29]+lzs*(zdata.xl[30]+lzs*zdata.xl[31]))
    msp[34] = zdata.xl[32]+lzs*(zdata.xl[33]+lzs*(zdata.xl[34]+lzs*zdata.xl[35]))
    msp[35] = zdata.xl[36]+lzs*(zdata.xl[37]+lzs*(zdata.xl[38]+lzs*zdata.xl[39]))
    msp[36] = zdata.xl[40]+lzs*(zdata.xl[41]+lzs*(zdata.xl[42]+lzs*zdata.xl[43]))
    msp[37] = max(0.90,1.1064+lzs*(0.415+0.18*lzs))
    msp[38] = max(1.0,1.19+lzs*(0.377+0.176*lzs))

    if (z > 0.010):
        msp[37] = min(msp[37],1.0)
        msp[38] = min(msp[38],1.10)

    msp[39] = max(0.1450,0.0977-lzs*(0.231+0.0753*lzs))
    msp[40] = min(0.240+lzs*(0.18+0.595*lzs),0.306+0.053*lzs)
    msp[41] = min(0.330+lzs*(0.132+0.218*lzs),
                  0.36250+0.062*lzs)
    msp[42] = (msp[33]+msp[34]*m2**msp[36])/(m2**0.40+msp[35]*m2**1.9)

    # Lbeta
    msp[43] = zdata.xl[44]+lzs*(zdata.xl[45]+lzs*(zdata.xl[46]+lzs*(zdata.xl[47]+lzs*zdata.xl[48])))
    msp[44] = zdata.xl[49]+lzs*(zdata.xl[50]+lzs*(zdata.xl[51]+lzs*(zdata.xl[52]+lzs*zdata.xl[53])))
    msp[45] = zdata.xl[54]+lzs*(zdata.xl[55]+lzs*zdata.xl[56])
    msp[46] = min(1.40,1.5135+0.3769*lzs)
    msp[46] = max(0.63550-0.4192*lzs,max(1.25,msp[46]))

    # Lhook
    msp[47] = zdata.xl[57]+lzs*(zdata.xl[58]+lzs*(zdata.xl[59]+lzs*zdata.xl[60]))
    msp[48] = zdata.xl[61]+lzs*(zdata.xl[62]+lzs*(zdata.xl[63]+lzs*zdata.xl[64]))
    msp[49] = zdata.xl[65]+lzs*(zdata.xl[66]+lzs*(zdata.xl[67]+lzs*zdata.xl[68]))
    msp[50] = zdata.xl[69]+lzs*(zdata.xl[70]+lzs*(zdata.xl[71]+lzs*zdata.xl[72]))
    msp[51] = min(1.40,1.5135+0.3769*lzs)
    msp[51] = max(0.63550-0.4192*lzs,max(1.25,msp[51]))

    # Rtms
    msp[52] = zdata.xr[1]+lzs*(zdata.xr[2]+lzs*(zdata.xr[3]+lzs*(zdata.xr[4]+lzs*zdata.xr[5])))
    msp[53] = zdata.xr[6]+lzs*(zdata.xr[7]+lzs*(zdata.xr[8]+lzs*(zdata.xr[9]+lzs*zdata.xr[10])))
    msp[54] = zdata.xr[11]+lzs*(zdata.xr[12]+lzs*(zdata.xr[13]+lzs*(zdata.xr[14]+lzs*zdata.xr[15])))
    msp[55] = zdata.xr[16]+lzs*(zdata.xr[17]+lzs*(zdata.xr[18]+lzs*zdata.xr[19]))
    msp[56] = zdata.xr[20]+lzs*(zdata.xr[21]+lzs*(zdata.xr[22]+lzs*zdata.xr[23]))
    msp[52] = msp[52]*msp[54]
    msp[53] = msp[53]*msp[54]
    msp[57] = zdata.xr[24]
    msp[58] = zdata.xr[25]+lzs*(zdata.xr[26]+lzs*(zdata.xr[27]+lzs*zdata.xr[28]))
    msp[59] = zdata.xr[29]+lzs*(zdata.xr[30]+lzs*(zdata.xr[31]+lzs*zdata.xr[32]))
    msp[60] = zdata.xr[33]+lzs*(zdata.xr[34]+lzs*(zdata.xr[35]+lzs*zdata.xr[36]))
    msp[61] = zdata.xr[37]+lzs*(zdata.xr[38]+lzs*(zdata.xr[39]+lzs*zdata.xr[40]))
    #
    msp[62] = max(0.0970-0.1072*(lz+3.0),max(0.097,min(0.1461,
                  0.14610+0.1237*(lz+2.0))))
    msp[62] = 10.0**msp[62]
    m2 = msp[62] + 0.10
    msp[63] = (msp[52]+msp[53]*msp[62]**msp[55])/(msp[54]+msp[62]**msp[56])
    msp[64] = (msp[57]*m2**3+msp[58]*m2**msp[61] + msp[59]*m2**(msp[61]+1.50))/(msp[60]+m2**5)

    # Ralpha
    msp[65] = zdata.xr[41]+lzs*(zdata.xr[42]+lzs*(zdata.xr[43]+lzs*zdata.xr[44]))
    msp[66] = zdata.xr[45]+lzs*(zdata.xr[46]+lzs*(zdata.xr[47]+lzs*zdata.xr[48]))
    msp[67] = zdata.xr[49]+lzs*(zdata.xr[50]+lzs*(zdata.xr[51]+lzs*zdata.xr[52]))
    msp[68] = zdata.xr[53]+lzs*(zdata.xr[54]+lzs*(zdata.xr[55]+lzs*zdata.xr[56]))
    msp[69] = zdata.xr[57]+lzs*(zdata.xr[58]+lzs*(zdata.xr[59]+lzs*(zdata.xr[60]+lzs*zdata.xr[61])))
    msp[70] = max(0.90,min(1.0,1.116+0.166*lzs))
    msp[71] = max(1.4770+0.296*lzs,min(1.6,-0.308-1.046*lzs))
    msp[71] = max(0.80,min(0.8-2.0*lzs,msp[71]))
    msp[72] = zdata.xr[62]+lzs*(zdata.xr[63]+lzs*zdata.xr[64])
    msp[73] = max(0.0650,0.0843-lzs*(0.0475+0.0352*lzs))
    msp[74] = 0.07360+lzs*(0.0749+0.04426*lzs)
    if (z < 0.0040):
        msp[74] = min(0.055,msp[74])

    msp[75] = max(0.0910, min(0.121, 0.136+0.0352*lzs))
    msp[76] = (msp[65]*msp[71]**msp[67])/(msp[66] + msp[71]**msp[68])

    if (msp[70] > msp[71]):
        msp[70] = msp[71]
        msp[75] = msp[76]

    # Rbeta
    msp[77] = zdata.xr[65]+lzs*(zdata.xr[66]+lzs*(zdata.xr[67]+lzs*zdata.xr[68]))
    msp[78] = zdata.xr[69]+lzs*(zdata.xr[70]+lzs*(zdata.xr[71]+lzs*zdata.xr[72]))
    msp[79] = zdata.xr[73]+lzs*(zdata.xr[74]+lzs*(zdata.xr[75]+lzs*zdata.xr[76]))
    msp[80] = zdata.xr[77]+lzs*(zdata.xr[78]+lzs*(zdata.xr[79]+lzs*zdata.xr[80]))
    msp[81] = zdata.xr[81]+lzs*(zdata.xr[82]+lzs*lzs*zdata.xr[83])
    if (z > 0.010):
        msp[81] = max(msp[81], 0.95)

    msp[82] = max(1.40, min(1.6, 1.6+lzs*(0.764+0.3322*lzs)))

    # Rgamma
    msp[83] = max(zdata.xr[84]+lzs*(zdata.xr[85]+lzs*(zdata.xr[86]+lzs*zdata.xr[87])),
                  zdata.xr[96]+lzs*(zdata.xr[97]+lzs*zdata.xr[98]))
    msp[84] = min(0.0, zdata.xr[88]+lzs*(zdata.xr[89]+lzs*(zdata.xr[90]+lzs*zdata.xr[91])))
    msp[84] = max(msp[84], zdata.xr[99]+lzs*(zdata.xr[100]+lzs*zdata.xr[101]))
    msp[85] = zdata.xr[92]+lzs*(zdata.xr[93]+lzs*(zdata.xr[94]+lzs*zdata.xr[95]))
    msp[85] = max(0.0, min(msp[85], 7.454+9.046*lzs))
    msp[86] = min(zdata.xr[102]+lzs*zdata.xr[103], max(2.0, -13.3-18.6*lzs))
    msp[87] = min(1.50,max(0.4,2.493+1.1475*lzs))
    msp[88] = max(1.0,min(1.27,0.8109-0.6282*lzs))
    msp[88] = max(msp[88],0.63550-0.4192*lzs)
    msp[89] = max(5.855420e-02,-0.27110-lzs*(0.5756+0.0838*lzs))

    # Rhook
    msp[90] = zdata.xr[104]+lzs*(zdata.xr[105]+lzs*(zdata.xr[106]+lzs*zdata.xr[107]))
    msp[91] = zdata.xr[108]+lzs*(zdata.xr[109]+lzs*(zdata.xr[110]+lzs*zdata.xr[111]))
    msp[92] = zdata.xr[112]+lzs*(zdata.xr[113]+lzs*(zdata.xr[114]+lzs*zdata.xr[115]))
    msp[93] = zdata.xr[116]+lzs*(zdata.xr[117]+lzs*(zdata.xr[118]+lzs*zdata.xr[119]))
    msp[94] = min(1.250,
              max(1.10,1.9848+lzs*(1.1386+0.3564*lzs)))
    msp[95] = 0.0630 + lzs*(0.0481 + 0.00984*lzs)
    msp[96] = min(1.30,max(0.45,1.2+2.45*lzs))

    # Lneta
    if (z > 0.00090):
        msp[97] = 10.0
    else:
        msp[97] = 20.0

    # Lbgb
    gbp[1] = zdata.xg[1]+lzs*(zdata.xg[2]+lzs*(zdata.xg[3]+lzs*zdata.xg[4]))
    gbp[2] = zdata.xg[5]+lzs*(zdata.xg[6]+lzs*(zdata.xg[7]+lzs*zdata.xg[8]))
    gbp[3] = zdata.xg[9]+lzs*(zdata.xg[10]+lzs*(zdata.xg[11]+lzs*zdata.xg[12]))
    gbp[4] = zdata.xg[13]+lzs*(zdata.xg[14]+lzs*(zdata.xg[15]+lzs*zdata.xg[16]))
    gbp[5] = zdata.xg[17]+lzs*(zdata.xg[18]+lzs*zdata.xg[19])
    gbp[6] = zdata.xg[20]+lzs*(zdata.xg[21]+lzs*zdata.xg[22])
    gbp[3] = gbp[3]**gbp[6]
    gbp[7] = zdata.xg[23]
    gbp[8] = zdata.xg[24]

    # Lbagb
    # set gbp[16] = 1.0 until it is reset later with an initial
    # call to Lbagbf using mass = zpars[2] and mhefl = 0.0
    gbp[9] = zdata.xg[25] + lzs*(zdata.xg[26] + lzs*zdata.xg[27])
    gbp[10] = zdata.xg[28] + lzs*(zdata.xg[29] + lzs*zdata.xg[30])
    gbp[11] = 15.0
    gbp[12] = zdata.xg[31]+lzs*(zdata.xg[32]+lzs*(zdata.xg[33]+lzs*zdata.xg[34]))
    gbp[13] = zdata.xg[35]+lzs*(zdata.xg[36]+lzs*(zdata.xg[37]+lzs*zdata.xg[38]))
    gbp[14] = zdata.xg[39]+lzs*(zdata.xg[40]+lzs*(zdata.xg[41]+lzs*zdata.xg[42]))
    gbp[15] = zdata.xg[43]+lzs*zdata.xg[44]
    gbp[12] = gbp[12]**gbp[15]
    gbp[14] = gbp[14]**gbp[15]
    gbp[16] = 1.0

    # Rgb
    gbp[17] = -4.67390-0.9394*lz
    gbp[17] = 10.0**gbp[17]
    gbp[17] = max(gbp[17],-0.041670+55.67*z)
    gbp[17] = min(gbp[17],0.47710-9329.21*z**2.94)
    gbp[18] = min(0.540,0.397+lzs*(0.28826+0.5293*lzs))
    gbp[19] = max(-0.14510,-2.2794-lz*(1.5175+0.254*lz))
    gbp[19] = 10.0**gbp[19]
    if (z > 0.0040):
        gbp[19] = max(gbp[19],0.73070+14265.1*z**3.395)

    gbp[20] = zdata.xg[45]+lzs*(zdata.xg[46]+lzs*(zdata.xg[47]+lzs*(zdata.xg[48]+lzs*(zdata.xg[49]+lzs*zdata.xg[50]))))
    gbp[21] = zdata.xg[51]+lzs*(zdata.xg[52]+lzs*(zdata.xg[53]+lzs*(zdata.xg[54]+lzs*zdata.xg[55])))
    gbp[22] = zdata.xg[56]+lzs*(zdata.xg[57]+lzs*(zdata.xg[58]+lzs*(zdata.xg[59]+lzs*(zdata.xg[60]+lzs*zdata.xg[61]))))
    gbp[23] = zdata.xg[62]+lzs*(zdata.xg[63]+lzs*(zdata.xg[64]+lzs*(zdata.xg[65]+lzs*zdata.xg[66])))

    # Ragb
    gbp[24] = min(0.991640-743.123*z**2.83,
                  1.04220+lzs*(0.13156+0.045*lzs))
    gbp[25] = zdata.xg[67]+lzs*(zdata.xg[68]+lzs*(zdata.xg[69]+lzs*(zdata.xg[70]+ lzs*(zdata.xg[71]+lzs*zdata.xg[72]))))
    gbp[26] = zdata.xg[73]+lzs*(zdata.xg[74]+lzs*(zdata.xg[75]+lzs*(zdata.xg[76]+lzs*zdata.xg[77])))
    gbp[27] = zdata.xg[78]+lzs*(zdata.xg[79]+lzs*(zdata.xg[80]+lzs*(zdata.xg[81]+lzs*(zdata.xg[82]+lzs*zdata.xg[83]))))
    gbp[28] = zdata.xg[84]+lzs*(zdata.xg[85]+lzs*(zdata.xg[86]+lzs*(zdata.xg[87]+lzs*zdata.xg[88])))
    gbp[29] = zdata.xg[89]+lzs*(zdata.xg[90]+lzs*(zdata.xg[91]+lzs*(zdata.xg[92]+lzs*(zdata.xg[93]+lzs*zdata.xg[94]))))
    gbp[30] = zdata.xg[95]+lzs*(zdata.xg[96]+lzs*(zdata.xg[97]+lzs*(zdata.xg[98]+lzs*(zdata.xg[99]+lzs*zdata.xg[100]))))
    m1 = zpars[2] - 0.20
    gbp[31] = gbp[29] + gbp[30]*m1
    gbp[32] = min(gbp[25]/zpars[2]**gbp[26],gbp[27]/zpars[2]**gbp[28])

    # Mchei
    gbp[33] = zdata.xg[101]**4
    gbp[34] = zdata.xg[102]*4.0

    # Mcagb
    gbp[35] = zdata.xg[103]+lzs*(zdata.xg[104]+lzs*(zdata.xg[105]+lzs*zdata.xg[106]))
    gbp[36] = zdata.xg[107]+lzs*(zdata.xg[108]+lzs*(zdata.xg[109]+lzs*zdata.xg[110]))
    gbp[37] = zdata.xg[111]+lzs*zdata.xg[112]
    gbp[35] = gbp[35]**4
    gbp[36] = gbp[36]*4.0
    gbp[37] = gbp[37]**4

    # Lhei
    # set gbp[41] = -1.0 until it is reset later with an initial
    # call to Lheif using mass = zpars[2] and mhefl = 0.0
    gbp[38] = zdata.xh[1]+lzs*zdata.xh[2]
    gbp[39] = zdata.xh[3]+lzs*zdata.xh[4]
    gbp[40] = zdata.xh[5]
    gbp[41] = -1.0
    gbp[42] = zdata.xh[6]+lzs*(zdata.xh[7]+lzs*zdata.xh[8])
    gbp[43] = zdata.xh[9]+lzs*(zdata.xh[10]+lzs*zdata.xh[11])
    gbp[44] = zdata.xh[12]+lzs*(zdata.xh[13]+lzs*zdata.xh[14])
    gbp[42] = gbp[42]**2
    gbp[44] = gbp[44]**2
    # Lhe
    gbp[45] = zdata.xh[15]+lzs*(zdata.xh[16]+lzs*zdata.xh[17])
    if (lzs > -1.0):
        gbp[46] = 1.0 - zdata.xh[19]*(lzs+1.0)**zdata.xh[18]
    else:
        gbp[46] = 1.0

    gbp[47] = zdata.xh[20]+lzs*(zdata.xh[21]+lzs*zdata.xh[22])
    gbp[48] = zdata.xh[23]+lzs*(zdata.xh[24]+lzs*zdata.xh[25])
    gbp[45] = gbp[45]**gbp[48]
    gbp[47] = gbp[47]**gbp[48]
    gbp[46] = gbp[46]/zpars[3]**0.10+(gbp[46]*gbp[47]-gbp[45])/zpars[3]**(gbp[48]+0.10)

    # Rmin

    gbp[49] = zdata.xh[26]+lzs*(zdata.xh[27]+lzs*(zdata.xh[28]+lzs*zdata.xh[29]))
    gbp[50] = zdata.xh[30]+lzs*(zdata.xh[31]+lzs*(zdata.xh[32]+lzs*zdata.xh[33]))
    gbp[51] = zdata.xh[34]+lzs*(zdata.xh[35]+lzs*(zdata.xh[36]+lzs*zdata.xh[37]))
    gbp[52] = 5.0+zdata.xh[38]*z**zdata.xh[39]
    gbp[53] = zdata.xh[40]+lzs*(zdata.xh[41]+lzs*(zdata.xh[42]+lzs*zdata.xh[43]))
    gbp[49] = gbp[49]**gbp[53]
    gbp[51] = gbp[51]**(2.0*gbp[53])

    # The
    # set gbp[57] = -1.0 until it is reset later with an initial
    # call to Thef using mass = zpars[2], mc = 0.0  and mhefl = 0.0
    gbp[54] = zdata.xh[44]+lzs*(zdata.xh[45]+lzs*(zdata.xh[46]+lzs*zdata.xh[47]))
    gbp[55] = zdata.xh[48]+lzs*(zdata.xh[49]+lzs*zdata.xh[50])
    gbp[55] = max(gbp[55],1.0)
    gbp[56] = zdata.xh[51]
    gbp[57] = -1.0
    gbp[58] = zdata.xh[52]+lzs*(zdata.xh[53]+lzs*(zdata.xh[54]+lzs*zdata.xh[55]))
    gbp[59] = zdata.xh[56]+lzs*(zdata.xh[57]+lzs*(zdata.xh[58]+lzs*zdata.xh[59]))
    gbp[60] = zdata.xh[60]+lzs*(zdata.xh[61]+lzs*(zdata.xh[62]+lzs*zdata.xh[63]))
    gbp[61] = zdata.xh[64]+lzs*zdata.xh[65]
    gbp[58] = gbp[58]**gbp[61]
    gbp[60] = gbp[60]**5

    # Tbl
    dum1 = zpars[2]/zpars[3]
    gbp[62] = zdata.xh[66]+lzs*zdata.xh[67]
    gbp[62] = -gbp[62]*numpy.log10(dum1)
    gbp[63] = zdata.xh[68]
    if (lzd > 0.0):
        gbp[64] = 1.0-lzd*(zdata.xh[69]+lzd*(zdata.xh[70]+lzd*zdata.xh[71]))
    else:
        gbp[64] = 1.0

    gbp[65] = 1.0-gbp[64]*dum1**gbp[63]
    gbp[66] = 1.0 - lzd*(zdata.xh[77] + lzd*(zdata.xh[78] + lzd*zdata.xh[79]))
    gbp[67] = zdata.xh[72] + lzs*(zdata.xh[73] + lzs*(zdata.xh[74] + lzs*zdata.xh[75]))
    gbp[68] = zdata.xh[76]

    # Lzahb
    gbp[69] = zdata.xh[80] + lzs*(zdata.xh[81] + lzs*zdata.xh[82])
    gbp[70] = zdata.xh[83] + lzs*(zdata.xh[84] + lzs*zdata.xh[85])
    gbp[71] = 15.0
    gbp[72] = zdata.xh[86]
    gbp[73] = zdata.xh[87]

    # Rzahb
    gbp[75] = zdata.xh[88] + lzs*(zdata.xh[89] + lzs*(zdata.xh[90] + lzs*zdata.xh[91]))
    gbp[76] = zdata.xh[92] + lzs*(zdata.xh[93] + lzs*(zdata.xh[94] + lzs*zdata.xh[95]))
    gbp[77] = zdata.xh[96] + lzs*(zdata.xh[97] + lzs*(zdata.xh[98] + lzs*zdata.xh[99]))

    # finish Lbagb
    mhefl = 0.0
    lx = lbagbf(zpars[2], mhefl)
    gbp[16] = lx

    # finish LHeI
    dum1 = 0.0
    lhefl = lheif(zpars[2],mhefl)
    gbp[41] = (gbp[38]*zpars[2]**gbp[39]-lhefl)/(numpy.exp(zpars[2]*gbp[40])*lhefl)

    # finish THe
    thefl = thef(zpars[2],dum1,mhefl)*tbgbf(zpars[2])
    gbp[57] = (thefl-gbp[54])/(gbp[54]*numpy.exp(gbp[56]*zpars[2]))

    # finish Tblf
    rb = ragbf(zpars[3],lheif(zpars[3],zpars[2]),mhefl)
    rr = 1.0 - rminf(zpars[3])/rb
    rr = max(rr,1.0e-12)
    gbp[66] = gbp[66]/(zpars[3]**gbp[67]*rr**gbp[68])

    # finish Lzahb
    gbp[74] = lhefl*lHef(zpars[2])

    kw = 0
    tm = 0.0
    tn = 0.0
    star(kw,zpars[2],zpars[2],tm,tn,tscls,lums,GB,zpars)
    zpars[9] = mcgbf(lums[3],GB,lums[6])
    zpars[10] = mcgbf(lums[4],GB,lums[6])

    # set the hydrogen and helium abundances
    zpars[11] = 0.760 - 3.0*z
    zpars[12] = 0.240 + 2.0*z

    # set constant for low-mass CHeB stars
    zpars[13] = (rminf(zpars[2])/
                rgbf(zpars[2],lzahbf(zpars[2],zpars[9],zpars[2])))

    zpars[14] = z**0.40
    #
    return
