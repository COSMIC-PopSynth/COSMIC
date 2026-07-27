"""
don_lim
=======

This example tests ``don_lim``, which selects the donor mass-loss-rate
prescription during Roche-lobe overflow.  The same initial binary is evolved
twice while changing only ``don_lim``:

* ``don_lim = -1``: Hurley et al. (2002), the COSMIC default.
* ``don_lim = -2``: Claeys et al. (2014).

The selected binary undergoes stable Roche-lobe overflow from the primary
without immediately entering a common envelope or merging.  Donor mass versus
time is the direct diagnostic; orbital separation shows the resulting binary
response.
"""

import copy
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np

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


DON_LIM_VALUES = [-1, -2]
DON_LIM_LABELS = {
    -1: "-1: Hurley+2002 default",
    -2: "-2: Claeys+2014",
}
LINE_WIDTHS = {-1: 2.8, -2: 2.2}

SEARCH_INITIALS = {
    "m1": 3.0,
    "m2": 0.9,
    "periods": np.array([2.55, 2.70, 2.7850600718448604, 2.95, 3.15]),
}

FINAL_TPHYS = 500.0
METALLICITY = 0.014
PLOT_DTP = 0.03
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


def evolve_don_lim(m1, m2, porb, don_lim, dtp=PLOT_DTP):
    """Evolve one binary while changing only don_lim."""
    settings = copy.deepcopy(BSEDict)
    settings["don_lim"] = don_lim

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
    """Return whether RLOF, CE, or merger-like termination occurred."""
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


def rlof_diagnostics(bpp, bcm):
    """Measure donor mass loss and accretor growth during the first RLOF episode."""
    rlof_start, rlof_end = first_rlof_interval(bpp)
    if not np.isfinite(rlof_start):
        return {
            "rlof_start": np.nan,
            "rlof_end": np.nan,
            "duration": np.nan,
            "donor": np.nan,
            "accretor": np.nan,
            "donor_mass_before": np.nan,
            "donor_mass_after": np.nan,
            "donor_loss": 0.0,
            "average_donor_mdot": np.nan,
            "peak_donor_mdot": np.nan,
            "accretor_gain": 0.0,
        }

    start_row = nearest_row(bpp, rlof_start)
    end_row = nearest_row(bpp, rlof_end)
    donor = 1 if start_row["RRLO_1"] >= start_row["RRLO_2"] else 2
    accretor = 2 if donor == 1 else 1

    donor_before = start_row[f"mass_{donor}"]
    donor_after = end_row[f"mass_{donor}"]
    accretor_before = start_row[f"mass_{accretor}"]
    accretor_after = end_row[f"mass_{accretor}"]
    duration = rlof_end - rlof_start
    donor_loss = max(0.0, donor_before - donor_after)
    accretor_gain = max(0.0, accretor_after - accretor_before)

    rlof_bcm = bcm[bcm["tphys"].between(rlof_start, rlof_end)].copy()
    peak_mdot = np.nan
    if len(rlof_bcm) > 2:
        masses = rlof_bcm[f"mass_{donor}"].to_numpy()
        times = rlof_bcm["tphys"].to_numpy()
        dt = np.diff(times)
        dm = -np.diff(masses)
        valid = dt > 0
        if valid.any():
            peak_mdot = np.nanmax(dm[valid] / dt[valid])

    return {
        "rlof_start": rlof_start,
        "rlof_end": rlof_end,
        "duration": duration,
        "donor": donor,
        "accretor": accretor,
        "donor_mass_before": donor_before,
        "donor_mass_after": donor_after,
        "donor_loss": donor_loss,
        "average_donor_mdot": donor_loss / duration if duration > 0 else np.nan,
        "peak_donor_mdot": peak_mdot,
        "accretor_gain": accretor_gain,
    }


def trial_metrics(m1, m2, porb):
    """Run both prescriptions for one trial binary."""
    metrics = []
    for don_lim in DON_LIM_VALUES:
        bpp, bcm = evolve_don_lim(m1, m2, porb, don_lim, dtp=0.08)
        flags = interaction_flags(bpp, bcm)
        diagnostics = rlof_diagnostics(bpp, bcm)
        metrics.append({**flags, **diagnostics, "final_sep": bcm.iloc[-1]["sep"]})
    return metrics


def select_binary():
    """Find a stable-RLOF binary with a visible don_lim response."""
    m1 = SEARCH_INITIALS["m1"]
    m2 = SEARCH_INITIALS["m2"]

    for porb in SEARCH_INITIALS["periods"]:
        metrics = trial_metrics(m1, m2, porb)

        if not all(row["rlof"] for row in metrics):
            continue
        if any(row["ce"] or row["merger"] for row in metrics):
            continue
        if min(row["duration"] for row in metrics) < 0.5:
            continue
        if max(row["donor_loss"] for row in metrics) < 0.5:
            continue

        duration_difference = abs(metrics[0]["duration"] - metrics[1]["duration"])
        final_sep_difference = abs(metrics[0]["final_sep"] - metrics[1]["final_sep"])
        accretor_gain_difference = abs(
            metrics[0]["accretor_gain"] - metrics[1]["accretor_gain"]
        )
        if duration_difference + 0.05 * final_sep_difference + accretor_gain_difference < 0.05:
            continue

        return {"m1": m1, "m2": m2, "porb": porb}

    raise RuntimeError("No clean stable-RLOF don_lim-sensitive binary found.")


def choose_time_window(results):
    """Zoom around the full RLOF episode with a little context."""
    starts = [row["diagnostics"]["rlof_start"] for row in results.values()]
    ends = [row["diagnostics"]["rlof_end"] for row in results.values()]
    finite_starts = [time for time in starts if np.isfinite(time)]
    finite_ends = [time for time in ends if np.isfinite(time)]

    if not finite_starts or not finite_ends:
        return 0.0, max(row["bcm"]["tphys"].max() for row in results.values())

    return max(0.0, min(finite_starts) - 0.8), max(finite_ends) + 1.0


selected = select_binary()

results = {}
for don_lim in DON_LIM_VALUES:
    bpp, bcm = evolve_don_lim(
        selected["m1"],
        selected["m2"],
        selected["porb"],
        don_lim,
        dtp=PLOT_DTP,
    )
    diagnostics = rlof_diagnostics(bpp, bcm)
    results[don_lim] = {"bpp": bpp, "bcm": bcm, "diagnostics": diagnostics}


fig, axes = plt.subplots(2, 1, figsize=(8.4, 5.8), sharex=True)
colors = plt.get_cmap("tab10").colors[: len(DON_LIM_VALUES)]

for index, (don_lim, color) in enumerate(zip(DON_LIM_VALUES, colors)):
    bcm = results[don_lim]["bcm"]
    diagnostics = results[don_lim]["diagnostics"]
    donor = int(diagnostics["donor"])

    axes[0].plot(
        bcm["tphys"],
        bcm[f"mass_{donor}"],
        color=color,
        lw=LINE_WIDTHS[don_lim],
        marker=None,
        label=DON_LIM_LABELS[don_lim],
        zorder=4 + index,
    )
    axes[1].plot(
        bcm["tphys"],
        bcm["sep"],
        color=color,
        lw=LINE_WIDTHS[don_lim],
        marker=None,
        zorder=4 + index,
    )

    for ax in axes:
        ax.axvline(diagnostics["rlof_start"], color=color, lw=0.9, alpha=0.35)
        ax.axvline(diagnostics["rlof_end"], color=color, lw=0.9, alpha=0.35)

x_min, x_max = choose_time_window(results)
for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=10)

axes[0].set_ylabel(r"Donor mass [$M_\odot$]", fontsize=12)
axes[1].set_ylabel(r"Orbital separation [$R_\odot$]", fontsize=12)
axes[1].set_xlabel("Time [Myr]", fontsize=12)
axes[0].set_title("don_lim: donor mass loss during stable RLOF", fontsize=13, pad=12)

axes[0].legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.18),
    fontsize=9,
    frameon=True,
    ncol=2,
    handlelength=2.0,
)

axes[1].text(
    0.02,
    0.08,
    "thin colored lines mark each run's RLOF start/end",
    transform=axes[1].transAxes,
    fontsize=8,
    color="0.25",
    ha="left",
    va="bottom",
    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
)

fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.9])
plt.show()
