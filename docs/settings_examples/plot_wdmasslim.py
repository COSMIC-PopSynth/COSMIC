"""
wd_mass_lim
===========

This example tests ``wd_mass_lim``, which only matters for a specific,
fairly rare event: two white dwarfs merging together, where the combined
mass pushes the surviving white dwarf over the Chandrasekhar mass. 


``wd_mass_lim=1`` always gives the same fixed neutron star mass, regardless
of how much mass the merger actually delivered. ``wd_mass_lim=0`` instead
lets the neutron star's mass scale with the combined mass of the merger, so
a bigger merger produces a heavier neutron star.

An  white dwarf just under the Chandrasekhar mass is paired with a range
of white dwarf companions, close enough together that they spiral in and
merge within a reasonable amount of time. The resulting neutron star mass
is plotted against the companion's mass, comparing both settings side by
side.
""" 

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

BSEDict = {
    'pts1': 0.001, 'pts2': 0.01, 'pts3': 0.02, 'zsun': 0.014, 'windflag': 3,
    'eddlimflag': 0, 'neta': 0.5, 'bwind': 0.0, 'hewind': 0.5, 'beta': 0.125,
    'xi': 0.5, 'acc2': 1.5, 'LBV_flag': 1, 'alpha1': [1.0, 1.0], 'lambdaf': 0.0,
    'ceflag': 1, 'cekickflag': 2, 'cemergeflag': 1, 'cehestarflag': 0, 'qcflag': 5,
    'qcrit_array': [0.0]*16, 'kickflag': 5, 'sigma': 265.0, 'bhflag': 1,
    'bhsigmafrac': 1.0, 'sigmadiv': -20.0, 'ecsn': 2.25, 'ecsn_mlow': 1.6,
    'aic': 1, 'ussn': 1, 'polar_kick_angle': 90.0,
    'natal_kick_array': [[-100.0, -100.0, -100.0, -100.0, 0.0], [-100.0, -100.0, -100.0, -100.0, 0.0]],
    'mm_mu_ns': 400.0, 'mm_mu_bh': 200.0, 'remnantflag': 4, 'fryer_mass_limit': 0,
    'mxns': 3.0, 'fryer_fmix': 1.0, 'fryer_mcrit_nsbh': 5.75, 'rembar_massloss': 0.5,
    'wd_mass_lim': 1, 'maltsev_mode': 0, 'maltsev_fallback': 0.5, 'maltsev_pf_prob': 0.1,
    'pisn': -2, 'ppi_co_shift': 0.0, 'ppi_extra_ml': 0.0, 'bhspinflag': 0,
    'bhspinmag': 0.0, 'grflag': 1, 'eddfac': 10, 'gamma': -2, 'don_lim': -1,
    'acc_lim': [-1, -1], 'smt_periastron_check': 0, 'tflag': 1, 'ST_tide': 1,
    'fprimc_array': [2.0/21.0]*16, 'ifflag': 1, 'wdflag': 1, 'epsnov': 0.001,
    'bdecayfac': 1, 'bconst': 3000, 'ck': 1000, 'rejuv_fac': 1.0, 'rejuvflag': 0,
    'bhms_coll_flag': 0, 'htpmb': 1, 'ST_cr': 1, 'rtmsflag': 0,
}

outdir = os.path.expanduser('~/projects/')
Z = 0.02
masses = np.round(np.logspace(np.log10(0.5), np.log10(50), 55), 4)
n = len(masses)

def vrotf(m, ST):
    m = np.asarray(m, float)
    if ST > 0:
        hi = (10.0*m**-0.0354)/(0.0389 + m**-7.95)
        lo = (13.4*m**-0.12)/(0.0389 + m**-7.95)
        return np.where(m > 6.35, hi, lo)
    return 330.0*m**3.3/(15.0 + m**3.45)

def run(st):
    BSEDict['ST_tide'] = st
    ibt = InitialBinaryTable.InitialBinaries(
        m1=list(masses), m2=[0.5]*n, porb=[1e9]*n, ecc=[0.0]*n,
        tphysf=[1.0]*n, kstar1=[1]*n, kstar2=[1]*n, metallicity=[Z]*n,
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict)
    first = bpp.reset_index().groupby('bin_num').first()
    return pd.DataFrame({'m': first['mass_1'].values, 'osp': first['omega_spin_1'].values,
                         'rad': first['rad_1'].values}, index=first.index).sort_values('m')

bse = run(0)
strk = run(1)
bse['vrot'] = bse['osp'] * bse['rad'] / 45.35
strk['vrot'] = strk['osp'] * strk['rad'] / 45.35
print('ST_tide=0 omega_spin range:', round(bse['osp'].min(), 3), '-', round(bse['osp'].max(), 1), flush=True)
print('ST_tide=1 omega_spin range:', round(strk['osp'].min(), 3), '-', round(strk['osp'].max(), 1), flush=True)

mm = np.logspace(np.log10(0.5), np.log10(50), 300)
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
ax[0].plot(bse['m'], bse['osp'], marker='o', ms=4, lw=1.6, color='C3', label='ST_tide=0 (BSE)')
ax[0].plot(strk['m'], strk['osp'], marker='o', ms=4, lw=1.6, color='C0', label='ST_tide=1 (StarTrack)')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlabel(r'ZAMS mass [$M_\odot$]')
ax[0].set_ylabel(r'Initial $\omega_{\rm spin}$ [COSMIC units]')
ax[0].set_title('ST_tide: initial spin rate')
ax[0].legend(frameon=False, fontsize=11)
ax[0].grid(alpha=0.25, which='both')

ax[1].plot(mm, vrotf(mm, 0), ls='--', color='0.5', lw=1.2, label='vrotf analytic (BSE)')
ax[1].plot(mm, vrotf(mm, 1), ls='--', color='0.2', lw=1.2, label='vrotf analytic (StarTrack)')
ax[1].plot(bse['m'], bse['vrot'], marker='o', ms=4, lw=0, color='C3', label='ST_tide=0 (COSMIC)')
ax[1].plot(strk['m'], strk['vrot'], marker='o', ms=4, lw=0, color='C0', label='ST_tide=1 (COSMIC)')
ax[1].set_xscale('log')
ax[1].set_xlabel(r'ZAMS mass [$M_\odot$]')
ax[1].set_ylabel(r'Initial $v_{\rm rot}$ [km/s]')
ax[1].set_title('ST_tide: rotation velocity check')
ax[1].legend(frameon=False, fontsize=11)
ax[1].grid(alpha=0.25)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'ST_tide_test.png'), dpi=160)
plt.show()
print('saved:', os.path.join(outdir, 'ST_tide_test.png'), flush=True)
