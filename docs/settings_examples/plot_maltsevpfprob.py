""" 
maltsev_pf_prob
================


This example tests ``maltsev_pf_prob``, which sets the chance that a star
landing in the Maltsev+2025 partial-fallback mass range actually forms a
black hole through partial fallback, instead of just collapsing straight to
a neutron star.

``maltsev_pf_prob=0`` means stars in that range always become a
neutron star. ``maltsev_pf_prob=1`` means they always undergo partial
fallback and become a black hole instead. Values in between apply that
outcome to a matching fraction of stars, chosen at random.

A group of progenitors chosen to sit exactly in the partial-fallback range
is evolved repeatedly at low metallicity,sweeping ``maltsev_pf_prob`` from
0 to 1. The left panel shows the measured fraction of black holes matching
the expected value directly. The other two panels show how the population
splits between a neutron star peak near 1.4 solar masses and a black hole
tail at higher mass, and how that split shifts as ``maltsev_pf_prob``
increases.
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
BSEDict['remnantflag'] = 6
BSEDict['maltsev_mode'] = 0
BSEDict['maltsev_fallback'] = 0.5

outdir = os.path.expanduser('~/projects/')
Z = 1e-4

def run(masses, pf):
    BSEDict['maltsev_pf_prob'] = pf
    n = len(masses)
    ibt = InitialBinaryTable.InitialBinaries(
        m1=list(masses), m2=[1.0]*n, porb=[1e10]*n, ecc=[0.0]*n,
        tphysf=[13700.0]*n, kstar1=[1]*n, kstar2=[1]*n, metallicity=[Z]*n,
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict)
    fin = bcm.reset_index().groupby('bin_num').last()
    prog = initC.reset_index().set_index('bin_num')['mass_1']
    return pd.DataFrame({'prog': prog.loc[fin.index].values, 'mrem': fin['mass_1'].values}, index=fin.index)

# finer mass grid (0.1 Msun steps instead of 0.5) so the deterministic
# NS/BH-outcome histograms below are built from ~5x more progenitor masses,
# smoothing out the discretization bumps in the hump region.
grid = np.round(np.arange(15.0, 55.01, 0.1), 3)
# pf_prob=0.0 and pf_prob=1.0 are both deterministic (the random draw either
# always fails or always succeeds), so these two runs carry zero sampling
# noise -- lo = guaranteed-NS outcome, hi = guaranteed-partial-fallback
# outcome, for every progenitor mass.
lo = run(grid, 0.0).sort_values('prog')
hi = run(grid, 1.0).sort_values('prog')
flip = (lo['mrem'].values < 1.6) & (hi['mrem'].values > 1.6)
window = lo['prog'].values[flip]
lo_win = lo[flip]['mrem'].values
hi_win = hi[flip]['mrem'].values
print('window masses:', len(window))

# Real COSMIC random draws across the full pf_prob range, like the original
# version of this panel -- the dashed y=x line is the theoretical reference,
# the markers are actual stochastic Evolve.evolve() output.
pfvals = np.round(np.arange(0.0, 1.01, 0.1), 2)
Nrep = 30
pop = np.repeat(np.sort(window), Nrep)
frac, err = [], []
for pf in pfvals:
    r = run(pop, pf)
    isbh = r['mrem'].values > 1.6
    p = isbh.mean()
    frac.append(p)
    err.append(np.sqrt(p * (1 - p) / len(r)))
frac, err = np.array(frac), np.array(err)

fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(18, 5.2))
ax0.plot([0, 1], [0, 1], ls='--', color='0.5', lw=1.2, label='$y=x$')
ax0.errorbar(pfvals, frac, yerr=err, marker='o', lw=1.6, capsize=3, color='C0',
             label='COSMIC')
ax0.set_xlabel('maltsev_pf_prob')
ax0.set_ylabel('BH fraction among window systems')
ax0.set_title('maltsev_pf_prob: BH fraction')
ax0.legend(frameon=False, fontsize=11)
ax0.grid(alpha=0.25)
ax0.set_xlim(-0.03, 1.03)
ax0.set_ylim(-0.03, 1.03)

# deterministic (noise-free) expected mass histogram: weighted mix of the
# guaranteed-NS and guaranteed-BH outcomes for every window mass, instead of
# rerunning with a random pf_prob draw per star.
bins = np.linspace(0, 12, 121)
hist_ns, edges = np.histogram(lo_win, bins=bins)
hist_bh, _ = np.histogram(hi_win, bins=bins)
centers = 0.5 * (edges[:-1] + edges[1:])
scale = 1200 / len(window)

colors3 = plt.cm.viridis(np.linspace(0.15, 0.85, 3))
for c, pf in zip(colors3, [0.2, 0.5, 0.8]):
    expected = scale * ((1 - pf) * hist_ns + pf * hist_bh)
    ax1.step(centers, expected, where='mid', lw=1.8, color=c, label=f'pf_prob = {pf}')
    ax2.step(centers, expected, where='mid', lw=1.8, color=c, label=f'pf_prob = {pf}')

ax1.axvline(1.4, ls=':', color='0.4', lw=1.2)
ax1.set_xlim(0, 2.5)
ax1.set_xlabel(r'Remnant mass [$M_\odot$]')
ax1.set_ylabel('expected count (analytical, N=1200)')
ax1.set_title('maltsev_pf_prob: NS peak (0-2.5 $M_\\odot$)')
ax1.legend(frameon=False, fontsize=11)
ax1.grid(alpha=0.25)

ax2.set_xlim(2.5, 12)
ax2.set_ylim(0, 200)
ax2.set_xlabel(r'Remnant mass [$M_\odot$]')
ax2.set_ylabel('expected count (analytical, N=1200)')
ax2.set_title('maltsev_pf_prob: BH tail (2.5-12 $M_\\odot$)')
ax2.legend(frameon=False, fontsize=11)
ax2.grid(alpha=0.25)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'maltsev_pf_prob_test.png'), dpi=160)
plt.show()
