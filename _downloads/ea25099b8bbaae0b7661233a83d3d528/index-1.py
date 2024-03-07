from cosmic.utils import a_from_p
from cosmic.sample.initialbinarytable import InitialBinaryTable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
final_kstar = np.linspace(0,14,15)
colors = {'green' : '#1b9e77', 'orange' : '#d95f02', 'purple' : '#7570b3'}
initC_logP, m_sin_logP, m_bin_logP, n_sin_logP, n_bin_logP = InitialBinaryTable.sampler('independent',
                                                                                        final_kstar1=final_kstar,
                                                                                        final_kstar2=final_kstar,
                                                                                        binfrac_model=1.0,
                                                                                        primary_model='kroupa01',
                                                                                        ecc_model='thermal',
                                                                                        porb_model='log_uniform',
                                                                                        qmin=-1,
                                                                                        SF_start=13700.0,
                                                                                        SF_duration=0.0,
                                                                                        met=0.02,
                                                                                        size=100000)
initC_Sana, m_sin_Sana, m_bin_Sana, n_sin_Sana, n_bin_Sana = InitialBinaryTable.sampler('independent',
                                                                                        final_kstar1=final_kstar,
                                                                                        final_kstar2=final_kstar,
                                                                                        binfrac_model=1.0,
                                                                                        primary_model='kroupa01',
                                                                                        ecc_model='sana12',
                                                                                        porb_model='sana12',
                                                                                        qmin=-1,
                                                                                        SF_start=13700.0,
                                                                                        SF_duration=0.0,
                                                                                        met=0.02,
                                                                                        size=100000)
initC_Moe, m_sin_Moe, m_bin_Moe, n_sin_Moe, n_bin_Moe = InitialBinaryTable.sampler('independent',
                                                                                   final_kstar1=final_kstar,
                                                                                   final_kstar2=final_kstar,
                                                                                   binfrac_model=1.0,
                                                                                   primary_model='kroupa01',
                                                                                   ecc_model='sana12',
                                                                                   porb_model='moe19',
                                                                                   qmin=-1,
                                                                                   SF_start=13700.0,
                                                                                   SF_duration=0.0,
                                                                                   met=0.02,
                                                                                   size=100000)
initC_logP['sep'] = a_from_p(p=initC_logP.porb, m1=initC_logP.mass_1, m2=initC_logP.mass_2)
initC_Sana['sep'] = a_from_p(p=initC_Sana.porb, m1=initC_Sana.mass_1, m2=initC_Sana.mass_2)
initC_Moe['sep'] = a_from_p(p=initC_Moe.porb, m1=initC_Moe.mass_1, m2=initC_Moe.mass_2)
fig = plt.figure(figsize = (15,6))
ax1 = plt.subplot(231)
ax2 = plt.subplot(232)
ax3 = plt.subplot(233)
ax4 = plt.subplot(234)
ax5 = plt.subplot(235)
ax6 = plt.subplot(236)
ax1.hist(np.log10(initC_logP.mass_1), bins = 20, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax1.hist(np.log10(initC_Sana.mass_1), bins = 20, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax1.hist(np.log10(initC_Moe.mass_1), bins = 20, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax1.set_xlabel(r'Log$_{10}$(M$_1$/M$_{\odot}$)', size=18)
ax1.set_ylabel('normalized counts', size=18)
ax1.legend(prop={'size' : 18})
ax2.hist(np.log10(initC_logP.porb), bins = 20, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax2.hist(np.log10(initC_Sana.porb), bins = 20, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax2.hist(np.log10(initC_Moe.porb), bins = 20, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax2.set_xlabel(r'Log$_{10}$(P$_{\rm{orb}}$/day)', size=18)
ax3.hist(initC_logP.ecc, bins = 10, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax3.hist(initC_Sana.ecc, bins = 10, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax3.hist(initC_Moe.ecc, bins = 10, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax3.set_xlabel('Eccentricity', size=18)
ax4.hist(initC_logP.mass_2/initC_logP.mass_1, bins = 20, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax4.hist(initC_Sana.mass_2/initC_Sana.mass_1, bins = 20, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax4.hist(initC_Moe.mass_2/initC_Moe.mass_1, bins = 20, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax4.set_xlabel(r'q=M$_1$/M$_2$', size=18)
ax4.set_ylabel('normalized counts', size=18)
ax5.hist(np.log10(initC_logP.sep), bins = 20, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax5.hist(np.log10(initC_Sana.sep), bins = 20, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax5.hist(np.log10(initC_Moe.sep), bins = 20, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax5.set_xlabel(r'Log$_{10}$(a/R$_{\odot}$)', size=18)
ax6.hist(np.log10(initC_logP.sep*(1-initC_logP.ecc)), bins = 20, histtype='step', density=True,
         lw=3, color=colors['purple'], label='independent')
ax6.hist(np.log10(initC_Sana.sep*(1-initC_Sana.ecc)), bins = 20, histtype='step', density=True,
         lw=3, color=colors['orange'], label='Sana+2012')
ax6.hist(np.log10(initC_Moe.sep*(1-initC_Moe.ecc)), bins = 20, histtype='step', density=True,
         lw=3, color=colors['green'], label='Moe+2019')
ax6.set_xlabel(r'Log$_{10}$(a(1-e)/R$_{\odot}$)', size=18)
fig.tight_layout()
fig.show()
