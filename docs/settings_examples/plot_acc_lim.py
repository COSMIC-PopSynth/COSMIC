"""
acc_lim
=======

This example tests ``acc_lim``, which limits how much mass each star can accrete
during Roche-lobe overflow.  The same initial binary is evolved for each
``acc_lim`` setting.  The selected binary undergoes stable Roche-lobe overflow
from the primary, so the accretor mass growth is a direct diagnostic of the
accretion prescription; the separation panel shows the accompanying orbital
response.

The positive settings are fraction-based accretion efficiencies, while negative
settings select COSMIC's built-in accretion-limit prescriptions.  In particular,
``[1.0, 1.0]`` and the default ``[-1, -1]`` are not generally equivalent.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

from cosmic.evolve import Evolve
from cosmic.sample import InitialBinaryTable

sys.path.append("..")
sys.path.append("/Users/lukewilner/cosmic-code/COSMIC-testing/docs")
import generate_default_bsedict

try:
    plt.style.use("../_static/gallery.mplstyle")
except OSError:
    pass

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


ACC_LIM_CASES = [
    (-1.0, -1.0),
    (-2.0, -2.0),
    (-3.0, -3.0),
    (-4.0, -4.0),
    (0.0, 0.0),
    (0.5, 0.5),
    (1.0, 1.0),
]

ACC_LIM_LABELS = {
    (-1.0, -1.0): "[-1, -1] default",
    (-2.0, -2.0): "[-2, -2]",
    (-3.0, -3.0): "[-3, -3]",
    (-4.0, -4.0): "[-4, -4]",
    (0.0, 0.0): "[0, 0]",
    (0.5, 0.5): "[0.5, 0.5]",
    (1.0, 1.0): "[1, 1]",
}

LINE_WIDTHS = {
    (-1.0, -1.0): 2.8,
    (-2.0, -2.0): 2.6,
    (-3.0, -3.0): 2.1,
    (-4.0, -4.0): 1.9,
    (0.0, 0.0): 2.4,
    (0.5, 0.5): 2.2,
    (1.0, 1.0): 2.0,
}

SEARCH_INITIALS = {
    "m1": 3.0,
    "m2": 1.2,
    "periods": np.array([2.55, 2.70, 2.852150339128391, 3.00, 3.20]),
}

FINAL_TPHYS = 500.0
METALLICITY = 0.014
PLOT_DTP = 0.05
RLOF_START_EVOL_TYPE = 3
RLOF_END_EVOL_TYPE = 4
CE_EVOL_TYPE = 7
TERMINAL_EVOL_TYPE = 6


def make_binary(m1, m2, porb, tphysf=FINAL_TPHYS):
    """Build one initial binary."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([m1]),
        m2=np.array([m2]),
        porb=np.array([porb]),
        ecc=np.zeros(1),
        tphysf=np.ones(1) * tphysf,
        kstar1=np.ones(1),
        kstar2=np.ones(1),
        metallicity=np.ones(1) * METALLICITY,
    )


def evolve_acc_lim(m1, m2, porb, acc_lim, dtp=PLOT_DTP):
    """Evolve one binary while changing only acc_lim."""
    settings = copy.deepcopy(BSEDict)
    settings["acc_lim"] = list(acc_lim)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=make_binary(m1, m2, porb),
            BSEDict=settings,
            dtp=dtp,
            randomseed=1,
        )

    return (
        bpp.sort_values("tphys").reset_index(drop=True),
        bcm.sort_values("tphys").reset_index(drop=True),
    )


def terminal_time(bpp):
    """Return the first terminal-event time, or NaN if none is recorded."""
    terminal_rows = bpp[bpp["evol_type"].eq(TERMINAL_EVOL_TYPE)]
    if terminal_rows.empty:
        return np.nan
    return terminal_rows["tphys"].iloc[0]


def first_rlof_interval(bpp):
    """Return the first RLOF start and end times."""
    starts = bpp[bpp["evol_type"].eq(RLOF_START_EVOL_TYPE)]
    if starts.empty:
        return np.nan, np.nan

    start_time = starts["tphys"].iloc[0]
    ends = bpp[
        bpp["evol_type"].eq(RLOF_END_EVOL_TYPE) & bpp["tphys"].ge(start_time)
    ]
    if ends.empty:
        return start_time, bpp["tphys"].max()
    return start_time, ends["tphys"].iloc[0]


def nearest_row(df, time):
    """Return the row closest to a requested time."""
    return df.iloc[(df["tphys"] - time).abs().argsort()[:1]].iloc[0]


def interaction_flags(bpp, bcm):
    """Return basic interaction/outcome flags from COSMIC output."""
    final = bcm.iloc[-1]
    return {
        "rlof": bool(bpp["evol_type"].eq(RLOF_START_EVOL_TYPE).any()),
        "ce": bool(bpp["evol_type"].eq(CE_EVOL_TYPE).any()),
        "merger": bool(
            bpp["evol_type"].eq(TERMINAL_EVOL_TYPE).any()
            or final["bin_state"] != 0
            or final["sep"] <= 0
        ),
    }


def rlof_diagnostics(bpp):
    """Measure donor loss and accretor growth during the first RLOF episode."""
    rlof_start, rlof_end = first_rlof_interval(bpp)
    if not np.isfinite(rlof_start):
        return {
            "rlof_start": np.nan,
            "rlof_end": np.nan,
            "donor": np.nan,
            "accretor": np.nan,
            "donor_loss": 0.0,
            "accretor_gain": 0.0,
            "systemic_mass_loss": 0.0,
            "effective_efficiency": np.nan,
        }

    start_row = nearest_row(bpp, rlof_start)
    end_row = nearest_row(bpp, rlof_end)

    donor = 1 if start_row["RRLO_1"] >= start_row["RRLO_2"] else 2
    accretor = 2 if donor == 1 else 1

    donor_start = start_row[f"mass_{donor}"]
    donor_end = end_row[f"mass_{donor}"]
    accretor_start = start_row[f"mass_{accretor}"]
    accretor_end = end_row[f"mass_{accretor}"]

    donor_loss = max(0.0, donor_start - donor_end)
    accretor_gain = max(0.0, accretor_end - accretor_start)
    systemic_loss = max(0.0, donor_loss - accretor_gain)
    efficiency = accretor_gain / donor_loss if donor_loss > 0 else np.nan

    return {
        "rlof_start": rlof_start,
        "rlof_end": rlof_end,
        "donor": donor,
        "accretor": accretor,
        "donor_loss": donor_loss,
        "accretor_gain": accretor_gain,
        "systemic_mass_loss": systemic_loss,
        "effective_efficiency": efficiency,
    }


def quick_metrics(m1, m2, porb):
    """Run the acc_lim cases on a trial binary and return selection metrics."""
    rows = []
    for acc_lim in ACC_LIM_CASES:
        bpp, bcm = evolve_acc_lim(m1, m2, porb, acc_lim, dtp=0.1)
        flags = interaction_flags(bpp, bcm)
        diagnostics = rlof_diagnostics(bpp)
        rows.append({**flags, **diagnostics})
    return rows


def select_binary():
    """Find a stable-RLOF binary that shows measurable acc_lim differences."""
    m1 = SEARCH_INITIALS["m1"]
    m2 = SEARCH_INITIALS["m2"]

    for porb in SEARCH_INITIALS["periods"]:
        metrics = quick_metrics(m1, m2, porb)

        if not all(row["rlof"] for row in metrics):
            continue
        if any(row["ce"] or row["merger"] for row in metrics):
            continue

        rlof_durations = [
            row["rlof_end"] - row["rlof_start"]
            for row in metrics
            if np.isfinite(row["rlof_start"])
        ]
        accretor_gains = [row["accretor_gain"] for row in metrics]

        if min(rlof_durations) < 0.5:
            continue
        if max(accretor_gains) - min(accretor_gains) < 0.5:
            continue
        if max(row["donor_loss"] for row in metrics) < 0.5:
            continue

        return {"m1": m1, "m2": m2, "porb": porb}

    raise RuntimeError("No clean stable-RLOF acc_lim-sensitive binary found.")


def choose_time_window(results):
    """Zoom around the RLOF episode while preserving before/after context."""
    starts = [row["diagnostics"]["rlof_start"] for row in results.values()]
    ends = [row["diagnostics"]["rlof_end"] for row in results.values()]
    finite_starts = [time for time in starts if np.isfinite(time)]
    finite_ends = [time for time in ends if np.isfinite(time)]

    if not finite_starts or not finite_ends:
        return 0.0, max(row["bcm"]["tphys"].max() for row in results.values())

    x_min = max(0.0, min(finite_starts) - 1.0)
    x_max = max(finite_ends) + 1.4
    return x_min, x_max


selected = select_binary()

results = {}
for acc_lim in ACC_LIM_CASES:
    bpp, bcm = evolve_acc_lim(
        selected["m1"],
        selected["m2"],
        selected["porb"],
        acc_lim,
        dtp=PLOT_DTP,
    )
    diagnostics = rlof_diagnostics(bpp)
    results[acc_lim] = {"bpp": bpp, "bcm": bcm, "diagnostics": diagnostics}


fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.0), sharex=True)
colors = plt.get_cmap("tab10").colors[: len(ACC_LIM_CASES)]

for index, (acc_lim, color) in enumerate(zip(ACC_LIM_CASES, colors)):
    bcm = results[acc_lim]["bcm"]
    accretor = int(results[acc_lim]["diagnostics"]["accretor"])
    label = ACC_LIM_LABELS[acc_lim]

    # Plot narrower or duplicated tracks later so overlapped cases remain visible.
    zorder = 3 + index
    linewidth = LINE_WIDTHS[acc_lim]

    axes[0].plot(
        bcm["tphys"],
        bcm[f"mass_{accretor}"],
        color=color,
        lw=linewidth,
        label=label,
        marker=None,
        zorder=zorder,
    )
    axes[1].plot(
        bcm["tphys"],
        bcm["sep"],
        color=color,
        lw=linewidth,
        marker=None,
        zorder=zorder,
    )

default_start = results[(-1.0, -1.0)]["diagnostics"]["rlof_start"]
default_end = results[(-1.0, -1.0)]["diagnostics"]["rlof_end"]
x_min, x_max = choose_time_window(results)

for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=10)
    ax.axvline(default_start, color="0.35", lw=1.0)
    ax.axvline(default_end, color="0.35", lw=1.0)

axes[1].set_yscale("log")
axes[1].yaxis.set_major_locator(FixedLocator([10, 20, 30, 40, 50]))
axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
axes[1].yaxis.set_minor_formatter(NullFormatter())
axes[0].set_ylabel(r"Accretor mass [$M_\odot$]", fontsize=12)
axes[1].set_ylabel(r"Orbital separation [$R_\odot$]", fontsize=12)
axes[1].set_xlabel("Time [Myr]", fontsize=12)

axes[0].set_title("Effect of acc_lim on stable Roche-lobe overflow", fontsize=13, pad=10)
axes[0].legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.18),
    fontsize=8,
    frameon=True,
    ncol=4,
    handlelength=1.7,
)

sep_y_bottom = axes[1].get_ylim()[0]
event_label_y = sep_y_bottom * 1.25
axes[1].text(
    default_start,
    event_label_y,
    "default RLOF starts",
    ha="left",
    va="bottom",
    fontsize=8,
    color="0.25",
)
axes[1].text(
    default_end,
    event_label_y,
    "default RLOF ends",
    ha="left",
    va="bottom",
    fontsize=8,
    color="0.25",
)

fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.9])
plt.show()
