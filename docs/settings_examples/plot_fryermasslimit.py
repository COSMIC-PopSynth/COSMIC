"""
fryer_mass_limit
=================
This example tests ``fryer_mass_limit``, which controls which mass the
Fryer+2012 remnant-mass prescriptions (``remnantflag=3`` rapid,
``remnantflag=4`` delayed) use as mass that fallback
accretion is built on top of.

``fryer_mass_limit=0`` uses the star's total mass at core collapse.
``fryer_mass_limit=1`` uses the star's CO core mass at core collapse instead.

A mass grid of single stars is evolved once per ``fryer_mass_limit`` value,
under both the rapid and delayed prescriptions, at two metallicities
(Z=0.002 and Z=0.02):  four panels in total. Each panel plots final BH mass
vs. initial ZAMS mass; overlapping curves mean the setting made no
difference for that area, while a visible split means the choice of
limiting mass changed the final BH mass.
""" 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

n_grid = 20

# two mass ranges - one for NS, one for BH
ns_masses = np.linspace(8.0, 25.0, n_grid)  
bh_masses = np.linspace(20.0, 100.0, n_grid)  

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

def get_remnant_masses(masses, BSEDict, metallicity, n_grid):
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
    bh_init, bh_final, ns_init, ns_final = [], [], [], []
    for i in range(n_grid):
        rows = bpp.loc[i]
        init_mass = rows.iloc[0]['mass_1']
        bh_rows = rows[rows['kstar_1'] == 14]
        ns_rows = rows[rows['kstar_1'] == 13]
        if len(bh_rows) > 0:
            bh_init.append(init_mass)
            bh_final.append(bh_rows.iloc[-1]['mass_1'])
        if len(ns_rows) > 0:
            ns_init.append(init_mass)
            ns_final.append(ns_rows.iloc[-1]['mass_1'])
    return bh_init, bh_final, ns_init, ns_final

# 2 rows (remnantflag 4 and 3), 2 columns (BH low Z, BH high Z) -- NS panels dropped
fig, axes = plt.subplots(2, 2, figsize=(11, 10))

for row, remnantflag in enumerate([4, 3]):
    label_rf = 'Delayed (remnantflag=4)' if remnantflag == 4 else 'Rapid (remnantflag=3)'

    for col, (remnant_type, metallicity) in enumerate([
        ('BH', 0.002),
        ('BH', 0.02),
    ]):
        ax = axes[row, col]
        label_z = f'Z={metallicity}'

        # use different mass range depending on remnant type
        masses = bh_masses if remnant_type == 'BH' else ns_masses

        # blue drawn thick/solid, red thin/dashed on top -- so an exact overlap
        # (fryer_mass_limit having no effect) is visible as a red dashed line
        # running down the middle of a thick blue line, instead of red just
        # hiding blue entirely.
        for fml, color, label, lw, ls, ms in [(0, 'blue', 'fryer_mass_limit=0 (total mass)', 4.5, '-', 9),
                                               (1, 'red',  'fryer_mass_limit=1 (core mass)', 1.8, '--', 4)]:
            BSEDict = BSEDict_base.copy()
            BSEDict['fryer_mass_limit'] = fml
            BSEDict['remnantflag'] = remnantflag

            bh_init, bh_final, ns_init, ns_final = get_remnant_masses(
                masses, BSEDict, metallicity, n_grid
            )

            if remnant_type == 'BH' and len(bh_init) > 0:
                ax.plot(bh_init, bh_final, 'o', linestyle=ls, color=color, lw=lw, markersize=ms, label=label)
            elif remnant_type == 'NS' and len(ns_init) > 0:
                ax.plot(ns_init, ns_final, 'o', linestyle=ls, color=color, lw=lw, markersize=ms, label=label)

        ax.set_xlabel('Initial Stellar Mass (Msun)')
        ax.set_ylabel(f'Final {remnant_type} Mass (Msun)')
        ax.set_title(f'fryer_mass_limit: {label_rf}, {label_z}')
        ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig('fryer_mass_limit_comprehensive.png', dpi=150)
plt.show()
