"""
``neta``
========

This example shows the effect of the ``neta`` setting on low- and
intermediate-mass single-star evolution. ``neta`` is the Reimers wind
coefficient, so it mainly affects cool evolved stars on the red giant branch
and asymptotic giant branch.

COSMIC requires positive ``neta`` values, so ``neta = 0.01`` is used as a
nearly-off control. The left panel focuses on the narrow mass range where the
final mass is most sensitive to ``neta``. The right panel shows the same final
mass comparison over the full initial-mass grid.
"""

import sys

sys.path.append("..")
import generate_default_bsedict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)
plt.style.use("../_static/gallery.mplstyle")
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.title_fontsize": 8,
})


neta_values = [0.01, 0.25, 0.5, 1.0, 2.0]
masses = np.round(np.arange(0.8, 8.05, 0.1), 2)


def make_single_star_grid():
    """Build the low/intermediate-mass single-star grid."""
    n_systems = len(masses)
    return InitialBinaryTable.InitialBinaries(
        m1=masses,
        m2=np.zeros(n_systems),
        porb=np.ones(n_systems) * -1.0,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * 13700.0,
        kstar1=np.ones(n_systems),
        kstar2=np.zeros(n_systems),
        metallicity=np.ones(n_systems) * 0.014,
    )


def neta_label(neta):
    """Format labels so the default and nearly-off cases are explicit."""
    if neta == 0.01:
        return "0.01 (near off)"
    if neta == 0.5:
        return "0.5 (default)"
    return f"{neta:g}"


def summarize_final_state(bcm, neta):
    """Collect final masses for one neta value."""
    final_bcm = bcm.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()

    summary = pd.DataFrame(index=final_bcm.index)
    summary["neta"] = neta
    summary["initial_mass_1"] = masses[summary.index.astype(int)]
    summary["final_mass_1"] = final_bcm["mass_1"]
    summary["mass_lost_1"] = summary["initial_mass_1"] - summary["final_mass_1"]
    summary["final_kstar_1"] = final_bcm["kstar_1"]
    return summary.reset_index(drop=True)


def evolve_neta_grid(grid):
    """Run the single-star grid for each neta value."""
    summaries = []

    for neta in neta_values:
        settings = BSEDict.copy()
        settings["neta"] = neta
        settings["random_seed"] = 1

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=grid,
            BSEDict=settings,
            dtp=25,
        )
        summaries.append(summarize_final_state(bcm, neta))

    return pd.concat(summaries, ignore_index=True)


def plot_final_mass(ax, data, title):
    """Plot final mass as a function of initial mass for all neta values."""
    for neta in neta_values:
        subset = data[data["neta"].eq(neta)].sort_values("initial_mass_1")
        ax.plot(
            subset["initial_mass_1"],
            subset["final_mass_1"],
            marker="o",
            markersize=4,
            label=neta_label(neta),
        )

    ax.set_xlabel("Initial mass [$M_\\odot$]")
    ax.set_ylabel("Final mass [$M_\\odot$]")
    ax.set_title(title)
    ax.legend(title="neta", fontsize=8, title_fontsize=9, markerscale=1.2)
    ax.grid(True, alpha=0.3)


grid = make_single_star_grid()
summaries = evolve_neta_grid(grid)
zoom = summaries[
    (summaries["initial_mass_1"] >= 0.9)
    & (summaries["initial_mass_1"] <= 1.5)
]

fig, axes = plt.subplots(1, 2, figsize=(10, 4.0))

plot_final_mass(
    axes[0],
    zoom,
    "Sensitive transition region",
)
plot_final_mass(
    axes[1],
    summaries,
    "Full mass grid",
)

plt.tight_layout()
plt.savefig("neta_flagtest.png", dpi=150)
plt.show()
