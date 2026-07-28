"""
maltsev_fallback 
================

This example tests ``maltsev_fallback``, which controls how much extra mass
a black hole gets to keep when it forms through "partial fallback." 

``maltsev_fallback=0`` means none of that extra material is kept, so the
remnant ends up close to the smallest possible neutron star mass.
``maltsev_fallback=1`` means all of it is kept, so the remnant grows all the
way up to nearly the star's full core mass. Values in between keep a
proportional amount, scaling linearly from one extreme to the other.

Bare single stars are swept across ZAMS mass 15-55 Msun at low metallicity,
with fallback forced to always happen so the result isn't muddied by
randomness. The left panel shows how the final remnant mass depends on the
progenitor's mass for a few different fallback settings, with the affected
mass range shaded. The right panel zooms in on a handful of progenitors and
shows remnant mass increasing steadily as ``maltsev_fallback`` goes from 0
to 1.""" 



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
BSEDict['maltsev_pf_prob'] = 1.0

outdir = os.path.expanduser('~/projects/')
Z = 1e-4
masses = np.round(np.arange(15.0, 55.01, 0.5), 3)
fvals = [0.0, 0.25, 0.5, 0.75, 1.0]
n = len(masses)

curves = {}
for f in fvals:
    BSEDict['maltsev_fallback'] = f
    ibt = InitialBinaryTable.InitialBinaries(
        m1=list(masses), m2=[1.0]*n, porb=[1e10]*n, ecc=[0.0]*n,
        tphysf=[13700.0]*n, kstar1=[1]*n, kstar2=[1]*n, metallicity=[Z]*n,
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=ibt, BSEDict=BSEDict)
    fin = bcm.reset_index().groupby('bin_num').last()
    prog = initC.reset_index().set_index('bin_num')['mass_1']
    s = pd.Series(fin['mass_1'].values, index=prog.loc[fin.index].values).sort_index()
    curves[f] = s

allm = curves[fvals[0]].index.values
M = pd.DataFrame({f: curves[f].reindex(allm) for f in fvals})
spread = (M.max(axis=1) - M.min(axis=1))
window = spread[spread > 0.1].index.values
top = np.sort(window[np.linspace(0, len(window)-1, min(4, len(window))).astype(int)]) if len(window) else allm[:0]
print('partial-fallback progenitors:', len(window), 'of', len(allm))

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(fvals)))
for c, f in zip(cmap, fvals):
    ax[0].plot(curves[f].index, curves[f].values, color=c, lw=1.8, label=fr'$f_{{\rm fb}}$ = {f}')
if len(window):
    ax[0].axvspan(window.min(), window.max(), color='0.85', alpha=0.5, zorder=0,
                  label='partial-fallback window')
ax[0].set_xlabel(r'ZAMS progenitor mass [$M_\odot$]')
ax[0].set_ylabel(r'Remnant mass [$M_\odot$]')
ax[0].set_title('maltsev_fallback: remnant mass vs progenitor')
ax[0].legend(frameon=False, fontsize=11)
ax[0].grid(alpha=0.25)

# COSMIC converts every remnant mass from baryonic to gravitational via
# baryonic_to_gravitational_mass() (assign_remnant.f) -- applied to the
# *whole* fallback-interpolated mass, not just the f=0/f=1 endpoints:
#   mrem = max(6.6667*(sqrt(1+0.3*mt) - 1), mt - rembar_massloss)
# The fallback law M_rem = M_NS + f_fb*(M_c,tot - M_NS) is linear in the
# baryonic mass mt, but this conversion is piecewise (a kink where the two
# branches cross), so the true f_fb -> gravitational-mass relation is not a
# straight line. Recover the baryonic core mass from the f=1 COSMIC point
# and re-run it through the same conversion to get the correct theory curve.
M_NS = 1.4
rembar_massloss = BSEDict['rembar_massloss']

def baryonic_to_grav(mt):
    mrem_quad = 6.6666667 * (np.sqrt(1.0 + 0.3 * mt) - 1.0)
    return np.maximum(mrem_quad, mt - rembar_massloss)

def grav_to_baryonic(mrem):
    # baryonic_to_grav is monotonic increasing in mt -> invert by bisection
    lo = np.full_like(mrem, M_NS, dtype=float)
    hi = np.full_like(mrem, 300.0, dtype=float)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        too_high = baryonic_to_grav(mid) > mrem
        hi = np.where(too_high, mid, hi)
        lo = np.where(too_high, lo, mid)
    return 0.5 * (lo + hi)

ff = np.linspace(0, 1, 200)
for c, mm in zip(plt.cm.plasma(np.linspace(0.1, 0.85, len(top))), top):
    mco_grav = M.loc[mm, 1.0]
    mco_bary = grav_to_baryonic(np.array([mco_grav]))[0]
    mt_bary_theory = M_NS + ff * (mco_bary - M_NS)
    mrem_theory = baryonic_to_grav(mt_bary_theory)
    ax[1].plot(ff, mrem_theory, color=c, lw=1.2, ls='--', zorder=1)
    ax[1].plot(fvals, M.loc[mm].values, color=c, marker='o', lw=0, ms=7, zorder=2,
               label=fr'{mm:.1f} $M_\odot$ ($M_{{\rm c,bary}}$={mco_bary:.1f})')
ax[1].set_xlabel(r'maltsev_fallback  $f_{\rm fb}$')
ax[1].set_ylabel(r'Remnant mass [$M_\odot$]')
ax[1].set_title(r'maltsev_fallback: remnant mass vs $f_{\rm fb}$')
ax[1].legend(frameon=False, fontsize=11, title='progenitor')
ax[1].grid(alpha=0.25)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'maltsev_fallback_test.png'), dpi=160)
plt.show()
