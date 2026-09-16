"""
``acc2``
========

This example tests ``acc2``, the Bondi-Hoyle wind-accretion factor.  The same
detached, wind-fed binary is evolved for every ``acc2`` value while keeping the
wind velocity fixed with ``beta = 0.125``.  The plot shows the companion mass
gained from wind accretion, so all curves start at zero and separate only when
the wind-fed phase begins.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
import generate_default_bsedict

try:
    plt.style.use("../_static/gallery.mplstyle")
except OSError:
    pass

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


ACC2_VALUES = [0.5, 1.0, 1.5, 3.0, 5.0]
ACC2_LABELS = {
    0.5: "0.5",
    1.0: "1.0",
    1.5: "1.5 (default)",
    3.0: "3.0",
    5.0: "5.0",
}

INITIAL_BINARY = {
    "m1": 25.0,
    "m2": 10.0,
    "porb": 10000.0,
    "ecc": 0.0,
    "metallicity": 0.014,
    "tphysf": 7.75,
}

FIXED_BETA = 0.125
PLOT_DTP = 0.01


def make_binary():
    """Build the detached wind-fed binary used for every acc2 value."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([INITIAL_BINARY["m1"]]),
        m2=np.array([INITIAL_BINARY["m2"]]),
        porb=np.array([INITIAL_BINARY["porb"]]),
        ecc=np.array([INITIAL_BINARY["ecc"]]),
        tphysf=np.array([INITIAL_BINARY["tphysf"]]),
        kstar1=np.ones(1),
        kstar2=np.ones(1),
        metallicity=np.array([INITIAL_BINARY["metallicity"]]),
    )


def evolve_acc2(acc2):
    """Evolve the binary while changing only acc2."""
    settings = copy.deepcopy(BSEDict)
    settings["beta"] = FIXED_BETA
    settings["acc2"] = acc2

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_binary(),
            BSEDict=settings,
            dtp=PLOT_DTP,
            randomseed=1,
        )

    return (
        bpp.sort_values("tphys").reset_index(drop=True),
        bcm.sort_values("tphys").reset_index(drop=True),
    )


def assert_detached(acc2, bpp, bcm):
    """Require max(RRLO_1) < 1 and max(RRLO_2) < 1 for every run."""
    max_rrlo_1 = max(bpp["RRLO_1"].max(), bcm["RRLO_1"].max())
    max_rrlo_2 = max(bpp["RRLO_2"].max(), bcm["RRLO_2"].max())
    if max_rrlo_1 >= 1.0 or max_rrlo_2 >= 1.0 or bpp["evol_type"].eq(3).any():
        raise RuntimeError(
            f"acc2={acc2:g} reached Roche-lobe overflow; choose a wider binary."
        )


results = {}
for acc2 in ACC2_VALUES:
    bpp, bcm = evolve_acc2(acc2)
    assert_detached(acc2, bpp, bcm)
    results[acc2] = bcm

m2_gained = [bcm["mass_2"].iloc[-1] - bcm["mass_2"].iloc[0] for bcm in results.values()]
if max(m2_gained) - min(m2_gained) < 0.1:
    raise RuntimeError("Selected binary does not show a visible acc2 wind-accretion effect.")

primary_mass_lost = [
    bcm["mass_1"].iloc[0] - bcm["mass_1"].iloc[-1] for bcm in results.values()
]
if max(primary_mass_lost) - min(primary_mass_lost) > 0.05:
    raise RuntimeError("Primary wind mass loss is not sufficiently similar across acc2 runs.")


fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
colors = plt.get_cmap("plasma")(np.linspace(0.12, 0.88, len(ACC2_VALUES)))
line_widths = [1.7, 1.9, 2.2, 2.5, 2.8]

for zorder, (acc2, color, lw) in enumerate(zip(ACC2_VALUES, colors, line_widths), start=3):
    bcm = results[acc2]
    delta_m2 = bcm["mass_2"] - bcm["mass_2"].iloc[0]
    ax.plot(
        bcm["tphys"],
        delta_m2,
        color=color,
        lw=lw,
        marker=None,
        label=ACC2_LABELS[acc2],
        zorder=zorder,
    )

ax.set_xlim(6.5, INITIAL_BINARY["tphysf"])
ax.grid(alpha=0.25)
ax.tick_params(labelsize=10)
ax.set_xlabel("Time [Myr]", fontsize=12)
ax.set_ylabel(r"Companion mass gained [$M_\odot$]", fontsize=12)
ax.set_title("Effect of acc2 on detached wind accretion", fontsize=13, pad=10)
ax.legend(
    loc="upper left",
    ncol=2,
    fontsize=8,
    frameon=True,
    handlelength=1.8,
)

plt.show()
