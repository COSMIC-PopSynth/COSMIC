"""
ppi_co_shift / ppi_extra_ml
============================
This example tests ``ppi_co_shift`` and ``ppi_extra_ml``, two settings that
adjust how much extra mass a very massive star loses right before it
collapses, in the mass range where pulsational pair-instability (PPI) mass
loss applies (``pisn=-4``, the Renzo+2022 prescription).

``ppi_co_shift`` shifts the star's core mass, in solar masses, before it is
looked up in the mass-loss table: a positive value makes the star behave
as if its core were more massive than it actually is when deciding how much
mass to lose. ``ppi_extra_ml`` instead adds a fixed extra amount of mass
loss, in solar masses, directly on top of whatever the table already
predicts.

A grid of very massive single stars is evolved at two metallicities, once
for each setting while the other is held at its default value of zero. Each
panel plots the resulting black hole mass against the star's initial mass,
for several values of the setting being tested.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

n_grid = 30
masses = np.linspace(20.0, 100.0, n_grid)

SSEDict = {'stellar_engine': 'sse'}

BSEDict_base = {
    "pts1": 0.001, "pts2": 0.01, "pts3": 0.02, "zsun": 0.014,
    "windflag": -1,
    "eddlimflag": 0, "neta": 0.5, "bwind": 0.0, "hewind": 0.5, "beta": 0.125,
    "xi": 0.5, "acc2": 1.5, "LBV_flag": 1, "alpha1": [1.0, 1.0],
    "lambdaf": 0.0, "ceflag": 1, "cekickflag": 2, "cemergeflag": 1,
    "cehestarflag": 0, "qcflag": 5,
    "qcrit_array": [0.0]*16,
    "kickflag": 5, "sigma": 265.0, "bhflag": 1, "bhsigmafrac": 1.0,
    "sigmadiv": -20.0, "ecsn": 2.25, "ecsn_mlow": 1.6, "aic": 1, "ussn": 1,
    "polar_kick_angle": 90.0,
    "natal_kick_array": [[-100.0,-100.0,-100.0,-100.0,0.0],[-100.0,-100.0,-100.0,-100.0,0.0]],
    "mm_mu_ns": 400.0, "mm_mu_bh": 200.0,
    "remnantflag": 4,
    "fryer_mass_limit": 0,
    "mxns": 3.0, "fryer_fmix": 1.0,
    "fryer_mcrit_nsbh": 5.75, "rembar_massloss": 0.5, "wd_mass_lim": 1,
    "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
    "pisn": -4,
    "ppi_co_shift": 0.0,
    "ppi_extra_ml": 0.0,
    "bhspinflag": 0, "bhspinmag": 0.0,
    "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
    "acc_lim": [-1, -1], "smt_periastron_check": 0, "tflag": 1, "ST_tide": 1,
    "fprimc_array": [2.0/21.0]*16,
    "ifflag": 1, "wdflag": 1, "epsnov": 0.001, "bdecayfac": 1,
    "bconst": 3000, "ck": 1000, "rejuv_fac": 1.0, "rejuvflag": 0,
    "bhms_coll_flag": 0, "htpmb": 1, "ST_cr": 1, "rtmsflag": 0
}

def get_remnant_masses(BSEDict, metallicity, n_grid):
    binary_grid = InitialBinaryTable.InitialBinaries(
        m1=masses,
        m2=np.ones(n_grid) * 0.1,
        porb=np.ones(n_grid) * 100000.0,
        ecc=np.zeros(n_grid),
        tphysf=np.ones(n_grid) * 13700.0,
        kstar1=np.ones(n_grid),
        kstar2=np.ones(n_grid),
        metallicity=np.ones(n_grid) * metallicity
    )
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict
    )
    init_masses, final_masses = [], []
    for i in range(n_grid):
        rows = bpp.loc[i]
        init_mass = rows.iloc[0]['mass_1']
        remnant_rows = rows[rows['kstar_1'].isin([13, 14])]
        if len(remnant_rows) > 0:
            init_masses.append(init_mass)
            final_masses.append(remnant_rows.iloc[-1]['mass_1'])
    return init_masses, final_masses

shift_values = {
    -10.0: 'red',
    -5.0:  'orange',
     0.0:  'green',
     5.0:  'blue',
     10.0: 'purple',
}

extra_ml_values = {
    0.0:  'green',
    2.0:  'blue',
    5.0:  'orange',
    10.0: 'red',
    20.0: 'purple',
}

metallicities = {
    0.002: '-',   # low metallicity solid line
    0.02:  '--',  # high metallicity dashed line
}

fig, axes = plt.subplots(1, 2, figsize=(19, 6))

remnantflag = 4
label_rf = 'Delayed (remnantflag=4)'

# ── panel 1: ppi_co_shift ────────────────────────────────────────────────
ax = axes[0]
for shift, color in shift_values.items():
    for metallicity, linestyle in metallicities.items():
        BSEDict = BSEDict_base.copy()
        BSEDict['remnantflag'] = remnantflag
        BSEDict['ppi_co_shift'] = shift
        BSEDict['ppi_extra_ml'] = 0.0
        init_masses, final_masses = get_remnant_masses(BSEDict, metallicity, n_grid)
        label = f'shift={shift}, Z={metallicity}'
        if len(init_masses) > 0:
            ax.plot(init_masses, final_masses, color=color,
                   linestyle=linestyle, marker='o', label=label)

ax.set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
ax.set_ylabel('Final Remnant Mass (Msun)', fontsize=11)
ax.set_title('ppi_co_shift (remnantflag=4)', fontsize=13)
ax.legend(fontsize=9, ncol=2, loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
ax.grid(True, alpha=0.3)

# ── panel 2: ppi_extra_ml ────────────────────────────────────────────────
ax = axes[1]
for extra_ml, color in extra_ml_values.items():
    for metallicity, linestyle in metallicities.items():
        BSEDict = BSEDict_base.copy()
        BSEDict['remnantflag'] = remnantflag
        BSEDict['ppi_co_shift'] = 0.0
        BSEDict['ppi_extra_ml'] = extra_ml
        init_masses, final_masses = get_remnant_masses(BSEDict, metallicity, n_grid)
        label = f'extra_ml={extra_ml}, Z={metallicity}'
        if len(init_masses) > 0:
            ax.plot(init_masses, final_masses, color=color,
                   linestyle=linestyle, marker='o', label=label)

ax.set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
ax.set_ylabel('Final Remnant Mass (Msun)', fontsize=11)
ax.set_title('ppi_extra_ml (remnantflag=4)', fontsize=13)
ax.legend(fontsize=9, ncol=2, loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ppi_flags_comprehensive.png', dpi=150)
plt.show()
