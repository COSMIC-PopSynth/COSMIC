"""
``ceflag`` and ``cemergeflag``
==============================

This example compares two common-envelope-related COSMIC flags using one
representative binary per flag. Each panel evolves the same initial binary
multiple times while changing only the flag shown in that panel, then plots
orbital period as a function of time.

Single-binary examples are useful for these flags because the important
behavior often happens at a boundary between common-envelope survival and
merger.

Each panel varies one COSMIC flag while holding the initial binary fixed.
Tracks show orbital-period evolution; endpoint markers indicate whether the
binary ultimately survives or terminates through merger/disruption. The time
axis is limited to the interaction window where the flag differences are most
visible.
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

BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)


EXAMPLES = [
    {
        "title": "ceflag: CE orbital energy",
        "flag_key": "ceflag",
        "flags": {
            0: "0: core",
            1: "1: total (default)",
        },
        "initial": {
            "m1": 32.0,
            "m2": 8.0,
            "porb": 200.0,
            "kstar1": 1,
            "kstar2": 1,
            "tphysf": 300.0,
        },
        "annotation": "same binary, different CE energy prescription",
        "plot_end": 7.0,
    },
    {
        "title": "cemergeflag: CE merger treatment",
        "flag_key": "cemergeflag",
        "flags": {
            0: "0: alternate",
            1: "1: default",
        },
        "initial": {
            "m1": 8.0,
            "m2": 8.0,
            "porb": 1500.0,
            "kstar1": 1,
            "kstar2": 1,
            "tphysf": 500.0,
        },
        "annotation": "0 permits CE evolution; 1 forces merger for this system",
        "plot_end": 45.0,
    },
]


def make_binary(initial):
    """Build one InitialBinaryTable row from an example dictionary."""
    return InitialBinaryTable.InitialBinaries(
        m1=np.array([initial["m1"]]),
        m2=np.array([initial["m2"]]),
        porb=np.array([initial["porb"]]),
        ecc=np.zeros(1),
        tphysf=np.ones(1) * initial["tphysf"],
        kstar1=np.ones(1) * initial["kstar1"],
        kstar2=np.ones(1) * initial["kstar2"],
        metallicity=np.ones(1) * 0.014,
    )


def evolve_example(example, flag_value):
    """Evolve one example binary for one flag value."""
    settings = copy.deepcopy(BSEDict)
    settings[example["flag_key"]] = flag_value

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="At least one of your initial binaries is starting in Roche Lobe Overflow:*",
            category=UserWarning,
        )
        bpp, bcm, initC, _ = Evolve.evolve(
            initialbinarytable=make_binary(example["initial"]),
            BSEDict=settings,
            dtp=0.2,
            randomseed=1,
        )

    return bpp, bcm


def orbital_period_track(bpp, bcm):
    """Combine sampled evolution rows with event rows into one real track."""
    sampled_rows = bcm[["tphys", "porb", "bin_state"]].copy()
    sampled_rows["source_order"] = 0

    event_rows = bpp[["tphys", "porb"]].copy()
    event_rows["bin_state"] = np.nan
    event_rows["source_order"] = np.arange(len(event_rows)) + 1

    track = pd.concat([sampled_rows, event_rows], ignore_index=True)
    track = track[track["porb"].gt(0)].drop_duplicates(
        subset=["tphys", "porb", "bin_state"], keep="first"
    )
    track = track.sort_values(["tphys", "source_order"]).reset_index(drop=True)

    time = track["tphys"].to_numpy()
    period = track["porb"].to_numpy()
    if not np.all(np.diff(time) >= 0):
        raise RuntimeError("Orbital-period track is not time sorted.")
    if not np.all(np.isfinite(period)) or not np.all(period > 0):
        raise RuntimeError("Orbital-period track contains invalid period values.")

    return track


def final_status(bcm):
    """Return a short final-status string."""
    final_row = bcm.sort_values("tphys").iloc[-1]
    if final_row["bin_state"] == 0 and final_row["porb"] > 0:
        return "survives"
    if final_row["bin_state"] == 1:
        return "merges"
    return "disrupts"


def plotted_track(track, plot_end):
    """Return only real COSMIC rows inside the panel time window."""
    if track.empty:
        return track

    plot_track = track[track["tphys"].le(plot_end)].copy()
    if plot_track.empty:
        plot_track = track.iloc[[0]].copy()

    return plot_track


def mark_final_outcome(
    ax,
    track,
    status,
    color,
    marker_size=70,
    zorder=5,
):
    """Mark the plotted endpoint with a consistent simple symbol."""
    if track.empty:
        return

    final_point = track.iloc[-1]
    if status == "survives":
        ax.scatter(
            final_point["tphys"],
            final_point["porb"],
            marker="o",
            s=marker_size * 0.75,
            facecolor="white",
            edgecolor=color,
            linewidths=1.8,
            zorder=zorder,
        )
    elif status == "merges":
        ax.scatter(
            final_point["tphys"],
            final_point["porb"],
            marker="x",
            s=marker_size,
            linewidths=2.2,
            color=color,
            zorder=zorder,
        )
    else:
        ax.scatter(
            final_point["tphys"],
            final_point["porb"],
            marker="^",
            s=marker_size,
            facecolor="white",
            edgecolor=color,
            linewidths=2.0,
            zorder=zorder,
        )


def first_event_time(bpp, event_type):
    """Return the first time for an event type, or NaN if absent."""
    rows = bpp[bpp["evol_type"].eq(event_type)]
    if rows.empty:
        return np.nan
    return rows["tphys"].iloc[0]


def first_rlo_time(bpp):
    """Return the first recorded RLO time, or NaN if absent."""
    rrlo = bpp[["RRLO_1", "RRLO_2"]].max(axis=1)
    rows = bpp[rrlo.ge(1.0)]
    if rows.empty:
        return np.nan
    return rows["tphys"].iloc[0]


def primary_event(bpp):
    """Return the single event label/time most relevant to the panel."""
    ce_time = first_event_time(bpp, 7)
    if np.isfinite(ce_time):
        return "CE/RLO", ce_time

    rlo_time = first_rlo_time(bpp)
    if np.isfinite(rlo_time):
        return "CE/RLO", rlo_time

    return "", np.nan


def annotate_event(ax, label, event_time):
    """Draw one readable event marker."""
    if not label or not np.isfinite(event_time):
        return

    ax.axvline(event_time, color="0.35", ls=":", lw=1.4, zorder=4)


def line_width(flag_value):
    """Use a consistent width for the two CE-flag tracks."""
    return 2.0


def endpoint_marker_size(flag_value):
    """Use a consistent endpoint marker size."""
    return 70


fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
axes = axes.ravel()

for ax, example in zip(axes, EXAMPLES):
    records = []
    reference_record = None

    for flag_value, label in example["flags"].items():
        bpp, bcm = evolve_example(example, flag_value)
        track = orbital_period_track(bpp, bcm)
        status = final_status(bcm)
        record = {
            "flag_value": flag_value,
            "label": label,
            "bpp": bpp,
            "bcm": bcm,
            "track": track,
            "status": status,
        }
        records.append(record)

        if flag_value == max(example["flags"]):
            reference_record = record

    event_label, event_time = primary_event(reference_record["bpp"])
    all_plotted_tracks = []
    plot_end = example["plot_end"]
    for plot_index, record in enumerate(records):
        track = plotted_track(record["track"], plot_end)
        all_plotted_tracks.append(track)
        linewidth = line_width(record["flag_value"])
        zorder = 3 + record["flag_value"]
        line, = ax.plot(
            track["tphys"],
            track["porb"],
            lw=linewidth,
            ls="-",
            label=record["label"],
            zorder=zorder,
        )
        record["color"] = line.get_color()
        mark_final_outcome(
            ax,
            track,
            record["status"],
            line.get_color(),
            marker_size=endpoint_marker_size(record["flag_value"]),
            zorder=5,
        )

    annotate_event(ax, event_label, event_time)
    combined = pd.concat(all_plotted_tracks, ignore_index=True)
    ax.set_yscale("log")
    if not combined.empty:
        max_time = combined["tphys"].max()
        ax.set_xlim(0, max_time * 1.08)
        min_porb = combined["porb"].min()
        max_porb = combined["porb"].max()
        ax.set_ylim(min_porb * 0.65, max_porb * 1.6)
    ax.set_title(example["title"], fontsize=11)
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Orbital period [days]")
    ax.grid(True, which="both", alpha=0.25)
    flag_legend = ax.legend(
        fontsize=7,
        loc="best",
        title=example["flag_key"],
        title_fontsize=8,
    )
    ax.text(
        0.02,
        0.03,
        example["annotation"],
        transform=ax.transAxes,
        fontsize=7,
        color="0.25",
    )
fig.legend(
    handles=[
        Line2D([0], [0], marker="o", color="0.2", markerfacecolor="white",
               lw=0, markersize=6, label="survives"),
        Line2D([0], [0], marker="x", color="0.2",
               lw=0, markersize=7, label="merges"),
        Line2D([0], [0], marker="^", color="0.2", markerfacecolor="white",
               lw=0, markersize=7, label="disrupts"),
        Line2D([0], [0], color="0.35", ls=":", lw=1.4, label="CE/RLO"),
    ],
    loc="lower center",
    ncol=5,
    bbox_to_anchor=(0.5, -0.02),
    fontsize=7,
)
plt.show()
