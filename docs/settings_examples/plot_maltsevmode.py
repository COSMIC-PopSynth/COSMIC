"""
maltsev_mode
============

This example tests ``maltsev_mode``, which controls whether the Maltsev+2025
remnant-mass prescription (``remnantflag=6``) is allowed to be used all the
way down to very low metallicities, or whether it gets capped at some
minimum metallicity instead.

``maltsev_mode=0`` applies the prescription at any metallicity, with no
lower limit. ``maltsev_mode=1`` treats any metallicity below 1/50th of solar
as if it were exactly 1/50th of solar. ``maltsev_mode=2`` does the same
thing but at a higher cutoff, 1/10th of solar. All three modes give
identical results everywhere above their respective cutoffs: they only
disagree for very metal-poor stars, and even then only for progenitor
masses that happen to sit near a boundary in the underlying table.

A grid of ZAMS mass (15-55 Msun) and metallicity (from far below solar to a
few times solar) is evolved once per mode. The left panel is a heatmap
showing where in that grid the three modes give the biggest disagreement in
remnant mass. The right panel takes the single most disagreement-prone
progenitor mass and plots its remnant mass against metallicity for all
three modes side by side.""" 

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

Zsun = 0.014

Z_low_boundary  = Zsun / 50
Z_high_boundary = Zsun / 10

Z_ultra_low   = np.logspace(np.log10(0.0001),                np.log10(Z_low_boundary  * 0.99), 8)
Z_low         = np.logspace(np.log10(Z_low_boundary  * 1.01), np.log10(Z_high_boundary * 0.99), 8)
Z_in_range    = np.logspace(np.log10(Z_high_boundary * 1.01), np.log10(Zsun),                   10)
Z_super_solar = np.logspace(np.log10(Zsun * 1.01),            np.log10(Zsun * 3.0),              6)

Z_grid = np.concatenate([Z_ultra_low, Z_low, Z_in_range, Z_super_solar])
Z_grid = np.sort(np.unique(np.round(Z_grid, 8)))

# maltsev_mode only changes the outcome for progenitors whose core mass sits
# near the M1/M2/M3 boundaries in the Maltsev+25 table -- a single fixed mass
# (the earlier version of this notebook used M1=25 only) can easily miss that
# entirely, which is exactly what happened. Sweep mass too, matching the
# range malsev_fallback_final.ipynb already showed is boundary-sensitive.
mass_grid = np.round(np.linspace(15.0, 55.0, 12), 2)

print(f"Testing {len(mass_grid)} masses x {len(Z_grid)} metallicities x 3 modes "
      f"= {len(mass_grid) * len(Z_grid) * 3} runs")
print(f"Z/Zsun range: {Z_grid.min()/Zsun:.4f} -- {Z_grid.max()/Zsun:.2f}\n")

# stellar_engine belongs in its own SSEDict, not BSEDict -- leaving it inside
# BSEDict forced COSMIC to silently override + warn about it on every single
# evolve() call below.
SSEDict = {'stellar_engine': 'sse'}

BSE_BASE = dict(
    pts1              = 0.001,
    pts2              = 0.01,
    pts3              = 0.02,
    zsun              = Zsun,
    windflag          = 3,
    eddlimflag        = 0,
    neta              = 0.5,
    bwind             = 0.0,
    hewind            = 0.5,
    beta              = 0.125,
    xi                = 0.5,
    acc2              = 1.5,
    LBV_flag          = 1,
    alpha1            = [1.0, 1.0],
    lambdaf           = 0.0,
    ceflag            = 1,
    cekickflag        = 2,
    cemergeflag       = 1,
    cehestarflag      = 0,
    qcflag            = 5,
    qcrit_array       = [0.0]*16,
    kickflag          = 5,
    sigma             = 265.0,
    bhflag            = 1,
    bhsigmafrac       = 1.0,
    sigmadiv          = -20.0,
    ecsn              = 2.25,
    ecsn_mlow         = 1.6,
    aic               = 1,
    ussn              = 1,
    polar_kick_angle  = 90.0,
    natal_kick_array  = [[-100.0, -100.0, -100.0, -100.0, 0.0], [-100.0, -100.0, -100.0, -100.0, 0.0]],
    remnantflag       = 6,
    mxns              = 3.0,
    fryer_mass_limit  = 0,
    fryer_fmix        = 1.0,
    fryer_mcrit_nsbh  = 5.75,
    rembar_massloss   = 0.5,
    wd_mass_lim       = 1,
    maltsev_mode      = 0,
    maltsev_fallback  = 0.5,
    # deterministic: every star in the partial-fallback window actually goes
    # through partial fallback (see malsev_fallback_final.ipynb for why 0.1
    # was wrong here -- it makes ~90% of in-window stars randomly resolve to
    # a plain 1.4 Msun NS, drowning the mode comparison in noise)
    maltsev_pf_prob   = 1.0,
    mm_mu_ns          = 400.0,
    mm_mu_bh          = 200.0,
    bhspinflag        = 0,
    bhspinmag         = 0.0,
    grflag            = 1,
    eddfac            = 10,
    gamma             = -2,
    don_lim           = -1,
    acc_lim           = [-1, -1],
    smt_periastron_check = 0,
    tflag             = 1,
    ST_tide           = 1,
    fprimc_array      = [2.0/21.0]*16,
    ifflag            = 1,
    wdflag            = 1,
    epsnov            = 0.001,
    bdecayfac         = 1,
    bconst            = 3000,
    ck                = 1000,
    rejuv_fac         = 1.0,
    rejuvflag         = 0,
    bhms_coll_flag    = 0,
    htpmb             = 1,
    ST_cr             = 1,
    rtmsflag          = 0,
    pisn              = -2,
    piflag            = 0,
    ppi_co_shift      = 0.0,
    ppi_extra_ml      = 0.0,
)

def make_ibt(M1, Z):
    return InitialBinaryTable.InitialBinaries(
        m1          = [M1],
        m2          = [1.0],
        porb        = [1e10],
        ecc         = [0.0],
        tphysf      = [13700.0],
        kstar1      = [1],
        kstar2      = [1],
        metallicity = [Z],
    )

MODES = [0, 1, 2]
results = []

for mode in MODES:
    BSEDict = {**BSE_BASE, 'maltsev_mode': mode}
    print(f"Running maltsev_mode = {mode} ...")

    for M1 in mass_grid:
        for Z in Z_grid:
            ibt = make_ibt(M1, Z)
            try:
                bpp, bcm, initC, kick_info = Evolve.evolve(
                    initialbinarytable=ibt,
                    BSEDict=BSEDict,
                    SSEDict=SSEDict,
                )
                # evol_type is a numeric end-state code, never the string
                # 'final' -- the last row of bpp for this single binary is
                # its final state.
                final = bpp.iloc[-1]
                results.append({
                    'maltsev_mode': mode,
                    'M1':           M1,
                    'Z':            Z,
                    'Z_over_Zsun':  Z / Zsun,
                    'mass1':        final['mass_1'],
                })
            except Exception as e:
                print(f"  WARN: mode={mode}, M1={M1}, Z={Z:.2e} failed: {e}")
                results.append({
                    'maltsev_mode': mode,
                    'M1':           M1,
                    'Z':            Z,
                    'Z_over_Zsun':  Z / Zsun,
                    'mass1':        np.nan,
                })

df = pd.DataFrame(results)
df.to_csv(os.path.expanduser('~/projects/maltsev_mode_results.csv'), index=False)
print(f"\nSaved results: {len(df)} rows")
print(f"Successful runs: {df['mass1'].notna().sum()} / {len(df)}\n")

if df['mass1'].isna().all():
    print("All runs failed — check BSEDict settings before plotting.")
else:
    pivot = df.pivot_table(index=['M1', 'Z_over_Zsun'], columns='maltsev_mode', values='mass1')
    pivot.columns = [f'mode_{c}' for c in pivot.columns]
    mode_cols = [c for c in pivot.columns if c.startswith('mode_')]
    pivot['max_delta'] = pivot[mode_cols].max(axis=1) - pivot[mode_cols].min(axis=1)

    n_sensitive = (pivot['max_delta'] > 0.05).sum()
    print(f"(mass, Z) points where maltsev_mode changes the outcome (>0.05 Msun difference): "
          f"{n_sensitive} / {len(pivot)}")

    in_range = pivot.reset_index()
    in_range = in_range[in_range['Z_over_Zsun'].between(0.1, 1.0)]
    max_diff_inrange = in_range['max_delta'].max()
    print(f"Sanity check — biggest difference within the calibrated range Z/Zsun=[1/10,1]: "
          f"{max_diff_inrange:.6f} Msun")
    print("  ✓ PASS: all modes identical inside the calibrated range"
          if max_diff_inrange < 1e-6 else
          "  ✗ WARN: modes differ inside the calibrated range — unexpected!")

    heat = pivot['max_delta'].unstack('Z_over_Zsun')

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax0 = axes[0]
    Zvals = heat.columns.values
    Mvals = heat.index.values
    im = ax0.pcolormesh(Zvals, Mvals, heat.values, shading='nearest', cmap='viridis')
    ax0.set_xscale('log')
    for zb in (1/50, 1/10, 1.0):
        ax0.axvline(zb, color='white', ls='--', lw=1, alpha=0.8)
    ax0.set_xlabel('Z / Z☉')
    ax0.set_ylabel(r'ZAMS mass [$M_\odot$]')
    ax0.set_title('maltsev_mode: sensitivity map')
    fig.colorbar(im, ax=ax0, label=r'biggest difference [$M_\odot$]')

    best_M1 = pivot['max_delta'].idxmax()[0]
    best_delta = pivot['max_delta'].max()
    print(f"\nMost mode-sensitive progenitor: M1={best_M1:.1f} Msun (biggest difference {best_delta:.3f} Msun)")

    colors = {0: '#2196F3', 1: '#FF9800', 2: '#4CAF50'}
    labels = {
        0: 'mode 0 – no metallicity limit',
        1: 'mode 1 – floor at Z/Zsun=1/50',
        2: 'mode 2 – floor at Z/Zsun=1/10',
    }

    ax1 = axes[1]
    sub = df[df['M1'] == best_M1]
    for mode in MODES:
        s = sub[sub['maltsev_mode'] == mode].sort_values('Z_over_Zsun')
        ax1.plot(s['Z_over_Zsun'], s['mass1'], color=colors[mode], marker='o', ms=4,
                  lw=2, label=labels[mode])
    ax1.set_xscale('log')
    for zb in (1/50, 1/10, 1.0):
        ax1.axvline(zb, color='0.6', ls='--', lw=1)
    ax1.set_xlabel('Z / Z☉')
    ax1.set_ylabel(r'Remnant mass [$M_\odot$]')
    ax1.set_title(f'maltsev_mode (M₁={best_M1:.1f} M$_\\odot$)')
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(os.path.expanduser('~/projects/maltsev_mode_test.png'), dpi=150, bbox_inches='tight')
    print("Figure saved: ~/projects/maltsev_mode_test.png")
    plt.show()
