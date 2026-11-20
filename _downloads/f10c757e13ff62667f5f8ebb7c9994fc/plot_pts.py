"""
``pts1``, ``pts2``, and ``pts3``
================================

This example shows how the ``pts1``, ``pts2``, and ``pts3`` timestep modifiers
affect single-star evolution in COSMIC. These are numerical timestep-resolution
parameters, not physical prescriptions. Smaller values use finer timesteps.

The broad grid compares final stellar masses from 1 to 50 solar masses for
three timestep choices. The zoomed grid repeats the comparison from 22 to 25
solar masses and includes one extra-fine timestep choice to highlight a more
timestep-sensitive mass range.
"""

import sys

sys.path.append("..")
import generate_default_bsedict

import matplotlib.pyplot as plt
import numpy as np

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)
plt.style.use("../_static/gallery.mplstyle")


N_GRID = 50
TPHYSF = 13700.0
METALLICITY = 0.014

masses = np.linspace(1.0, 50.0, N_GRID)
zoom_masses = np.linspace(22.0, 25.0, 31)

pts_values = [
    {"pts1": 0.0005, "pts2": 0.005, "pts3": 0.01},
    {"pts1": 0.001, "pts2": 0.01, "pts3": 0.02},
    {"pts1": 0.002, "pts2": 0.02, "pts3": 0.04},
]

zoom_pts_values = [
    {"pts1": 0.00025, "pts2": 0.0025, "pts3": 0.005},
    {"pts1": 0.0005, "pts2": 0.005, "pts3": 0.01},
    {"pts1": 0.001, "pts2": 0.01, "pts3": 0.02},
    {"pts1": 0.002, "pts2": 0.02, "pts3": 0.04},
]


def make_single_star_grid(mass_grid):
    """Build a COSMIC single-star grid using the standard single-star setup."""
    n_systems = len(mass_grid)
    return InitialBinaryTable.InitialBinaries(
        m1=mass_grid,
        m2=np.zeros(n_systems),
        porb=np.ones(n_systems) * -1.0,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * TPHYSF,
        kstar1=np.ones(n_systems),
        kstar2=np.zeros(n_systems),
        metallicity=np.ones(n_systems) * METALLICITY,
    )


def pts_label(pts):
    """Format a legend label from the three timestep modifiers."""
    return f"pts=({pts['pts1']}, {pts['pts2']}, {pts['pts3']})"


def extract_final_mass_1(bcm, n_systems):
    """Return final primary masses, ordered by COSMIC binary number."""
    final = (
        bcm.sort_values(["bin_num", "tphys"])
        .groupby("bin_num")["mass_1"]
        .last()
        .reindex(np.arange(n_systems))
    )

    if final.isna().any():
        missing = final[final.isna()].index.to_list()
        raise ValueError(f"Missing final bcm rows for bin_num values: {missing}")

    return final.to_numpy()


def final_masses_for_pts(mass_grid, pts_grid):
    """Evolve one mass grid for each pts choice and collect final masses."""
    binary_grid = make_single_star_grid(mass_grid)
    final_masses = {}

    for pts in pts_grid:
        settings = BSEDict.copy()
        settings.update(pts)

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=binary_grid,
            BSEDict=settings,
            dtp=0,
        )

        final_masses[pts_label(pts)] = extract_final_mass_1(bcm, len(mass_grid))

    return final_masses


def plot_final_mass_panel(ax, mass_grid, final_masses, title):
    """Plot final mass as a function of initial mass for a set of pts choices."""
    for label, final_mass in final_masses.items():
        ax.plot(mass_grid, final_mass, marker="o", markersize=4, label=label)

    ax.set_xlabel("Initial Mass [$M_\\odot$]", fontsize=12)
    ax.set_ylabel("Final Mass [$M_\\odot$]", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)


final_masses = final_masses_for_pts(masses, pts_values)
zoom_final_masses = final_masses_for_pts(zoom_masses, zoom_pts_values)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

plot_final_mass_panel(
    axes[0],
    masses,
    final_masses,
    "Effect of pts Timesteps on Final Mass",
)
plot_final_mass_panel(
    axes[1],
    zoom_masses,
    zoom_final_masses,
    "Zoomed Final Mass Near Sensitive Region",
)

plt.tight_layout()
plt.savefig("pts_flagtest.png", dpi=150)
plt.show()
