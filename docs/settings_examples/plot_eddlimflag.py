"""
``eddlimflag``
==============

This example tests ``eddlimflag``, which changes the metallicity dependence of
wind mass loss for stars near the Eddington limit.  The same grid of very wide
binaries is evolved at each metallicity for both flag values:

* ``eddlimflag = 0``: no Eddington-dependent metallicity correction, the COSMIC
  default.
* ``eddlimflag = 1``: Giacobbo et al. (2018) correction.

The companions are low mass and the orbits are very wide, so the primary evolves
effectively as a single massive star.  The figure compares the primary mass at
the last non-compact stellar stage before compact-object formation; this keeps
the endpoint definition consistent across the grid.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
import generate_default_bsedict

try:
    plt.style.use("../_static/gallery.mplstyle")
except OSError:
    pass

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


EDDLIMFLAG_VALUES = [0, 1]
EDDLIMFLAG_LABELS = {
    0: "0: default",
    1: "1: Eddington correction",
}

METALLICITIES = [0.0002, 0.002, 0.014]
METALLICITY_LABELS = {
    0.0002: "Z = 0.0002",
    0.002: "Z = 0.002",
    0.014: "Z = 0.014",
}

MASS_GRID = np.concatenate(
    [
        np.arange(20.0, 151.0, 5.0),
        np.array([170.0, 200.0]),
    ]
)
FINAL_TPHYS = 13700.0
WIDE_PERIOD = 1.0e6
COMPANION_MASS = 0.1


def make_wide_binary_grid(metallicity):
    """Build a wide-binary grid that isolates primary wind evolution."""
    n_systems = len(MASS_GRID)
    return InitialBinaryTable.InitialBinaries(
        m1=MASS_GRID,
        m2=np.ones(n_systems) * COMPANION_MASS,
        porb=np.ones(n_systems) * WIDE_PERIOD,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * FINAL_TPHYS,
        kstar1=np.ones(n_systems),
        kstar2=np.ones(n_systems),
        metallicity=np.ones(n_systems) * metallicity,
    )


def evolve_grid(metallicity, eddlimflag):
    """Evolve one mass-metallicity grid while changing only eddlimflag."""
    settings = copy.deepcopy(BSEDict)
    settings["eddlimflag"] = eddlimflag
    settings["random_seed"] = 1

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_wide_binary_grid(metallicity),
            BSEDict=settings,
            dtp=10.0,
            randomseed=1,
        )

    return (
        bpp.sort_values(["bin_num", "tphys"]).reset_index(drop=True),
        bcm.sort_values(["bin_num", "tphys"]).reset_index(drop=True),
    )


def endpoint_row(group):
    """Use the primary mass just before compact-object formation."""
    ordered = group.sort_values(["tphys", "evol_type"])
    primary_sn = ordered[ordered["evol_type"].eq(15)]
    if len(primary_sn):
        return primary_sn.iloc[0]

    stellar_rows = ordered[ordered["kstar_1"].lt(13)]
    if len(stellar_rows):
        return stellar_rows.iloc[-1]

    return ordered.iloc[-1]


def summarize_final_state(bpp, metallicity, eddlimflag):
    """Collect consistently defined primary endpoint quantities for one grid."""
    endpoint = pd.DataFrame(
        [endpoint_row(group) for _, group in bpp.groupby("bin_num", sort=True)]
    ).set_index("bin_num")

    summary = pd.DataFrame(index=endpoint.index)
    summary["initial_mass_1"] = MASS_GRID[summary.index.astype(int)]
    summary["metallicity"] = metallicity
    summary["eddlimflag"] = eddlimflag
    summary["endpoint_mass_1"] = endpoint["mass_1"].to_numpy()
    summary["endpoint_kstar_1"] = endpoint["kstar_1"].astype(int).to_numpy()
    summary["endpoint_evol_type"] = endpoint["evol_type"].astype(int).to_numpy()
    summary["endpoint_tphys"] = endpoint["tphys"].to_numpy()
    summary["total_mass_lost"] = (
        summary["initial_mass_1"] - summary["endpoint_mass_1"]
    )

    # Verify that the setup remains an effectively single-star wind test.
    summary["RLOF_or_CE"] = (
        bpp.groupby("bin_num")["evol_type"]
        .apply(lambda values: bool(values.isin([3, 7]).any()))
        .values
    )
    return summary.reset_index(drop=True)


summaries = []
for metallicity in METALLICITIES:
    for eddlimflag in EDDLIMFLAG_VALUES:
        bpp, bcm = evolve_grid(metallicity, eddlimflag)
        summaries.append(summarize_final_state(bpp, metallicity, eddlimflag))

summary = pd.concat(summaries, ignore_index=True)

if summary["RLOF_or_CE"].any():
    raise RuntimeError("The selected wide-binary grid experienced RLOF or CE.")

difference_table = summary.pivot_table(
    index=["initial_mass_1", "metallicity"],
    columns="eddlimflag",
    values="endpoint_mass_1",
)
difference_table["abs_difference"] = (
    difference_table[1] - difference_table[0]
).abs()

if difference_table["abs_difference"].max() < 1.0:
    raise RuntimeError(
        "eddlimflag curves overlap too strongly; extend the high-mass grid."
    )

fig, ax = plt.subplots(figsize=(10.0, 5.0), constrained_layout=True)
colors = plt.get_cmap("tab10").colors
metallicity_colors = dict(zip(METALLICITIES, colors[: len(METALLICITIES)]))
line_styles = {0: "-", 1: "--"}
line_widths = {0: 2.5, 1: 2.4}

for metallicity in METALLICITIES:
    for eddlimflag in EDDLIMFLAG_VALUES:
        subset = summary[
            summary["metallicity"].eq(metallicity)
            & summary["eddlimflag"].eq(eddlimflag)
        ].sort_values("initial_mass_1")
        ax.plot(
            subset["initial_mass_1"],
            subset["endpoint_mass_1"],
            color=metallicity_colors[metallicity],
            linestyle=line_styles[eddlimflag],
            lw=line_widths[eddlimflag],
            marker=None,
        )

ax.set_xlabel(r"Initial primary mass [$M_\odot$]", fontsize=12)
ax.set_ylabel(r"Pre-SN primary mass [$M_\odot$]", fontsize=12)
ax.grid(alpha=0.3)
ax.tick_params(labelsize=10)
ax.set_xlim(MASS_GRID.min(), MASS_GRID.max())

metallicity_handles = [
    Line2D([0], [0], color=metallicity_colors[metallicity], lw=2.5)
    for metallicity in METALLICITIES
]
flag_handles = [
    Line2D([0], [0], color="0.2", linestyle=line_styles[flag], lw=line_widths[flag])
    for flag in EDDLIMFLAG_VALUES
]

metallicity_legend = ax.legend(
    metallicity_handles,
    [METALLICITY_LABELS[metallicity] for metallicity in METALLICITIES],
    loc="upper left",
    fontsize=8.5,
    frameon=True,
    ncol=1,
    handlelength=2.7,
)
ax.add_artist(metallicity_legend)
ax.legend(
    flag_handles,
    [EDDLIMFLAG_LABELS[flag] for flag in EDDLIMFLAG_VALUES],
    loc="upper center",
    bbox_to_anchor=(0.56, 1.02),
    fontsize=8.5,
    frameon=True,
    ncol=2,
    handlelength=2.7,
)

plt.show()
