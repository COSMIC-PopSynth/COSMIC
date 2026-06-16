import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve

n_grid = 30
masses = np.linspace(8.0, 100.0, n_grid)

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
    "mxns": 3.0,
    "fryer_fmix": 1.0,
    "fryer_mcrit_nsbh": 5.75, "rembar_massloss": 0.5,
    "wd_mass_lim": 1,
    "maltsev_mode": 0, "maltsev_fallback": 0.5, "maltsev_pf_prob": 0.1,
    "pisn": -2, "ppi_co_shift": 0.0, "ppi_extra_ml": 0.0, "bhspinflag": 0,
    "bhspinmag": 0.0, "grflag": 1, "eddfac": 10, "gamma": -2, "don_lim": -1,
    "acc_lim": [-1, -1], "smt_periastron_check": 0, "tflag": 1, "ST_tide": 1,
    "fprimc_array": [2.0/21.0]*16,
    "ifflag": 1, "wdflag": 1, "epsnov": 0.001, "bdecayfac": 1,
    "bconst": 3000, "ck": 1000, "rejuv_fac": 1.0, "rejuvflag": 0,
    "bhms_coll_flag": 0, "htpmb": 1, "ST_cr": 1, "rtmsflag": 0
}

binary_grid = InitialBinaryTable.InitialBinaries(
    m1=masses,
    m2=np.ones(n_grid) * 0.1,
    porb=np.ones(n_grid) * 100000.0,
    ecc=np.zeros(n_grid),
    tphysf=np.ones(n_grid) * 13700.0,
    kstar1=np.ones(n_grid),
    kstar2=np.ones(n_grid),
    metallicity=np.ones(n_grid) * 0.002
)

def get_remnant_masses(BSEDict, n_grid):
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=binary_grid, BSEDict=BSEDict, SSEDict=SSEDict
    )
    ns_init, ns_final = [], []
    bh_init, bh_final = [], []
    for i in range(n_grid):
        rows = bpp.loc[i]
        init_mass = rows.iloc[0]['mass_1']
        remnant_rows = rows[rows['kstar_1'].isin([13, 14])]
        if len(remnant_rows) > 0:
            final_mass = remnant_rows.iloc[-1]['mass_1']
            kstar = remnant_rows.iloc[-1]['kstar_1']
            if kstar == 13:
                ns_init.append(init_mass)
                ns_final.append(final_mass)
            elif kstar == 14:
                bh_init.append(init_mass)
                bh_final.append(final_mass)
    return ns_init, ns_final, bh_init, bh_final

# two modes — positive and negative
positive_values = {
    0.1:  'red',
    0.5:  'orange',   # default
    1.0:  'green',
    2.0:  'blue',
}

negative_values = {
    -0.1: 'red',
    -0.3: 'orange',
    -0.5: 'green',
    -0.7: 'blue',
}

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Effect of rembar_massloss on Remnant Masses', fontsize=14)

# ── top row: positive values ──────────────────────────────────────────────────
for rml, color in positive_values.items():
    BSEDict = BSEDict_base.copy()
    BSEDict['rembar_massloss'] = rml
    ns_init, ns_final, bh_init, bh_final = get_remnant_masses(BSEDict, n_grid)
    label = f'rembar={rml} (default)' if rml == 0.5 else f'rembar={rml}'

    if len(ns_init) > 0:
        axes[0, 0].plot(ns_init, ns_final, 'o-', color=color, label=label)
    if len(bh_init) > 0:
        axes[0, 1].plot(bh_init, bh_final, 'o-', color=color, label=label)

    print(f"rembar_massloss={rml}: {len(ns_init)} NS, {len(bh_init)} BH")

axes[0, 0].set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
axes[0, 0].set_ylabel('Final NS Mass (Msun)', fontsize=11)
axes[0, 0].set_title('NS Masses — Positive values\n(higher = more neutrino mass loss)', fontsize=11)
axes[0, 0].legend(fontsize=9)
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
axes[0, 1].set_ylabel('Final BH Mass (Msun)', fontsize=11)
axes[0, 1].set_title('BH Masses — Positive values\n(higher = more neutrino mass loss)', fontsize=11)
axes[0, 1].legend(fontsize=9)
axes[0, 1].grid(True, alpha=0.3)

# ── bottom row: negative values ───────────────────────────────────────────────
for rml, color in negative_values.items():
    BSEDict = BSEDict_base.copy()
    BSEDict['rembar_massloss'] = rml
    ns_init, ns_final, bh_init, bh_final = get_remnant_masses(BSEDict, n_grid)
    label = f'rembar={rml}'

    if len(ns_init) > 0:
        axes[1, 0].plot(ns_init, ns_final, 'o-', color=color, label=label)
    if len(bh_init) > 0:
        axes[1, 1].plot(bh_init, bh_final, 'o-', color=color, label=label)

    print(f"rembar_massloss={rml}: {len(ns_init)} NS, {len(bh_init)} BH")

axes[1, 0].set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
axes[1, 0].set_ylabel('Final NS Mass (Msun)', fontsize=11)
axes[1, 0].set_title('NS Masses — Negative values\n(more negative = less mass kept)', fontsize=11)
axes[1, 0].legend(fontsize=9)
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].set_xlabel('Initial Stellar Mass (Msun)', fontsize=11)
axes[1, 1].set_ylabel('Final BH Mass (Msun)', fontsize=11)
axes[1, 1].set_title('BH Masses — Negative values\n(more negative = less mass kept)', fontsize=11)
axes[1, 1].legend(fontsize=9)
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('rembar_massloss.png', dpi=150)
plt.show()
