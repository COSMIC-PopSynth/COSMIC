""" 
remnantflag
=========== 

This example tests ``remnantflag``, which chooses which published recipe is
used to turn a collapsing star's mass into the mass of the neutron star or
black hole it leaves behind.

The eight settings correspond to eight different published recipes: Hurley+
2000, Belczynski+2002, Belczynski+2008, the Fryer+2012 rapid and delayed
models, Mandel & Mueller 2020, Maltsev+2025, and Fryer+2022.

A grid of single stars is evolved once for each recipe, and each panel
plots the resulting neutron star or black hole mass against the star's
initial mass for that recipe. A shaded band marks the traditional "mass
gap" between 2 and 5 solar masses, so it's easy to see which recipes leave
that gap empty and which ones fill it in.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

n_grid = 30
masses = np.linspace(5.0, 100.0, n_grid)  # extended down to 5 Msun

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
    "pisn": -2, "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0, "bhspinflag": 0,
    "bhspinmag": 0.0, "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
    "acc_lim": [-1, -1], "smt_periastron_check": 0, "tflag": 1, "ST_tide": 1,
    "fprimc_array": [2.0/21.0]*16,
    "ifflag": 1, "wdflag": 1, "epsnov": 0.001, "bdecayfac": 1,
    "bconst": 3000, "ck": 1000, "rejuv_fac": 1.0, "rejuvflag": 0,
    "bhms_coll_flag": 0, "htpmb": 1, "ST_cr": 1, "rtmsflag": 0
}

remnant_flags = {
    0: ('Hurley+2000',           'blue',   0.002),  # use low metallicity
    1: ('Belczynski+2002',       'orange', 0.002),
    2: ('Belczynski+2008',       'green',  0.002),
    3: ('Fryer+2012 rapid',      'red',    0.002),
    4: ('Fryer+2012 delayed',    'purple', 0.002),
    5: ('Mandel & Mueller 2020', 'brown',  0.002),
    6: ('Maltsev+2025',          'pink',   0.02),   # use solar metallicity
    7: ('Fryer+2022',            'gray',   0.002),
}

fig, axes = plt.subplots(4, 2, figsize=(14, 20))
axes = axes.flatten()

for idx, (rflag, (label, color, metallicity)) in enumerate(remnant_flags.items()):
    ax = axes[idx]

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

    BSEDict = BSEDict_base.copy()
    BSEDict['remnantflag'] = rflag

    # remnantflag=5 needs specific settings per docs
    if rflag == 5:
        BSEDict['mxns'] = 2.0
        BSEDict['rembar_massloss'] = 0.1
        BSEDict['kickflag'] = 6

    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict
    )

    init_masses = []
    final_masses = []
    kstar_types = []
    for i in range(n_grid):
        rows = bpp.loc[i]
        init_mass = rows.iloc[0]['mass_1']
        remnant_rows = rows[rows['kstar_1'].isin([13, 14])]
        if len(remnant_rows) > 0:
            final_mass = remnant_rows.iloc[-1]['mass_1']
            kstar = remnant_rows.iloc[-1]['kstar_1']
            init_masses.append(init_mass)
            final_masses.append(final_mass)
            kstar_types.append(kstar)

    ns_init  = [init_masses[i] for i in range(len(kstar_types)) if kstar_types[i] == 13]
    ns_final = [final_masses[i] for i in range(len(kstar_types)) if kstar_types[i] == 13]
    bh_init  = [init_masses[i] for i in range(len(kstar_types)) if kstar_types[i] == 14]
    bh_final = [final_masses[i] for i in range(len(kstar_types)) if kstar_types[i] == 14]

    if len(ns_init) > 0:
        ax.plot(ns_init, ns_final, 'o-', color='steelblue', label='Neutron Star')
    if len(bh_init) > 0:
        ax.plot(bh_init, bh_final, 'o-', color=color, label='Black Hole')

    ax.axhspan(2, 5, alpha=0.1, color='gray', label='Mass gap (2-5 Msun)')
    ax.set_xlabel('Initial Stellar Mass (Msun)', fontsize=10)
    ax.set_ylabel('Final Remnant Mass (Msun)', fontsize=10)
    ax.set_title(f'remnantflag={rflag}: {label} (Z={metallicity})', fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    print(f"remnantflag={rflag} ({label}): {len(ns_init)} NS, {len(bh_init)} BH found")

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig('remnantflag.png', dpi=150)
plt.show()

