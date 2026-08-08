"""
``LBV_flag``
============

This example shows the effect of the ``LBV_flag`` on the evolution of massive stars.

The ``LBV_flag`` controls which prescription to use for LBV-like mass loss for
stars that are above the Humphreys-Davidson limit. The first figures show the
evolution of a few single stars on the HR diagram and highlight the LBV regime
where we'd expect to see differences with this flag. The final figure compares
the final primary mass over a wider initial-mass grid, showing how the LBV wind
prescription changes the amount of mass retained by massive single stars.
"""

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
## You'll want to edit this part locally to use your own BSEDict and style sheet!

import sys
sys.path.append("..")
import generate_default_bsedict
BSEDict = generate_default_bsedict.get_default_BSE_settings(to_python=True)

import matplotlib.pyplot as plt
plt.style.use("../_static/gallery.mplstyle")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as consts

from cosmic.sample import InitialBinaryTable
from cosmic.evolve import Evolve
from cosmic.output import kstar_translator

def LBV_limit(T_eff):
    """ Compute the luminosity of the LBV limit given a certain temperature """
    return np.maximum(np.sqrt(1e10 * u.Rsun**2 * u.Lsun * (4 * np.pi * consts.sigma_sb * (T_eff * u.K)**4)).to(u.Lsun), np.repeat(6e5 * u.Lsun, len(T_eff)))


# create a small grid of single stars
hr_masses = [30, 33, 36, 40, 45, 50]
grid = InitialBinaryTable.InitialBinaries(
    m1=hr_masses,
    m2=np.zeros(len(hr_masses)),
    porb=np.ones(len(hr_masses))*-1.0,
    ecc=np.ones(len(hr_masses))*0.0,
    tphysf=np.ones(len(hr_masses))*13700.0,
    kstar1=np.ones(len(hr_masses)),
    kstar2=np.ones(len(hr_masses)),
    metallicity=np.ones(len(hr_masses))*0.014*0.01
)

# loop over different flag choices
for LBV_flag, label in zip([0, 1, 2], ["No LBV mass loss", "LBV mass loss following Hurley+2000",
                                       "LBV mass loss following Belczynski+2008"]):
    # evolve with updated BSEDict
    settings = BSEDict.copy()
    settings["LBV_flag"] = LBV_flag
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=grid, BSEDict=settings,
        dtp=0
    )

    # plot a HR diagram of the evolution of these stars, colour coded by kstar
    fig, ax = plt.subplots(figsize=(10, 16))

    for i, bin_num in enumerate(np.unique(bcm["bin_num"])):
        log_teff = np.log10(bcm[bcm["bin_num"] == bin_num]["teff_1"].values)
        log_lum = np.log10(bcm[bcm["bin_num"] == bin_num]["lum_1"].values)
        kstar = bcm[bcm["bin_num"] == bin_num]["kstar_1"].values

        mask = kstar < 10

        ax.scatter(log_teff[mask], log_lum[mask], c=[kstar_translator[k]["colour"] for k in kstar[mask]],
                   s=5)
        ax.plot(log_teff[mask], log_lum[mask], color='grey', alpha=0.5, lw=0.5, zorder=-1)

        # annotate the first point with the initial mass
        ax.annotate(f'{bcm[bcm["bin_num"] == bin_num]["mass_1"].values[0]:.0f} M$_\\odot$', (log_teff[0], log_lum[0]),
                    fontsize=14, ha="right", color='black', va='top')
        
    # annotate the various kstar types
    for k in bcm["kstar_1"].unique():
        if k < 10:
            ax.scatter([], [], color=kstar_translator[k]["colour"], label=kstar_translator[k]["short"], s=50)
    ax.legend()

    # make some fairly reasonable temperature range
    T_eff_range = np.logspace(3.5, 5.4, 1000)
    
    # plot a line at the HD limit and fill the area
    ax.plot(np.log10(T_eff_range), np.log10(LBV_limit(T_eff_range).value), color="grey", linestyle="dotted")
    ax.fill_between(np.log10(T_eff_range), np.log10(LBV_limit(T_eff_range).value), 9, color="black", lw=2, alpha=0.03)
    
    # annotate the regime with a custom location
    ax.annotate("LBV Regime", xy=(0.77, 0.93), xycoords="axes fraction", ha="center", va="center",
                fontsize=20, color="grey", zorder=10)

    ax.set(
        xlabel=r"$\log_{10}(T_{\rm eff} \, [\rm K])$",
        ylabel=r"$\log_{10}(L \, [\rm L_{\odot}])$",
        title=f"LBV_flag = {LBV_flag}\n{label}",
        ylim=(5, 6.2)
    )

    ax.invert_xaxis()
    
    plt.show()


# compare the final masses over a wider initial-mass grid
LBV_flag_labels = {
    0: "0: no LBV winds",
    1: "1: Hurley+2000",
    2: "2: Belczynski+2008",
}
LBV_flag_values = list(LBV_flag_labels)
mass_grid = np.arange(25.0, 81.0, 1.0)


def make_final_mass_grid():
    """Build the massive single-star grid used for the final-mass comparison."""
    n_systems = len(mass_grid)
    return InitialBinaryTable.InitialBinaries(
        m1=mass_grid,
        m2=np.zeros(n_systems),
        porb=np.ones(n_systems) * -1.0,
        ecc=np.zeros(n_systems),
        tphysf=np.ones(n_systems) * 13700.0,
        kstar1=np.ones(n_systems),
        kstar2=np.zeros(n_systems),
        metallicity=np.ones(n_systems) * 0.014,
    )


def summarize_final_masses(bcm, LBV_flag):
    """Collect final primary masses for one LBV_flag value."""
    final_bcm = bcm.sort_values(["bin_num", "tphys"]).groupby("bin_num").last()

    summary = pd.DataFrame(index=final_bcm.index)
    summary["LBV_flag"] = LBV_flag
    summary["initial_mass_1"] = mass_grid[summary.index.astype(int)]
    summary["final_mass_1"] = final_bcm["mass_1"]
    return summary.reset_index(drop=True)


def evolve_final_mass_grid(final_mass_grid):
    """Run the same single-star grid for each LBV_flag value."""
    summaries = []

    for LBV_flag in LBV_flag_values:
        settings = BSEDict.copy()
        settings["LBV_flag"] = LBV_flag
        settings["random_seed"] = 1

        bpp, bcm, initC, kick_info = Evolve.evolve(
            initialbinarytable=final_mass_grid,
            BSEDict=settings,
            dtp=0,
        )
        summaries.append(summarize_final_masses(bcm, LBV_flag))

    return pd.concat(summaries, ignore_index=True)


final_mass_grid = make_final_mass_grid()
final_mass_summaries = evolve_final_mass_grid(final_mass_grid)

fig, ax = plt.subplots(figsize=(8.2, 6.0))

for LBV_flag in LBV_flag_values:
    subset = final_mass_summaries[final_mass_summaries["LBV_flag"].eq(LBV_flag)]
    ax.plot(
        subset["initial_mass_1"],
        subset["final_mass_1"],
        label=LBV_flag_labels[LBV_flag],
    )

ax.set_xlabel("Initial primary mass [$M_\\odot$]")
ax.set_ylabel("Final primary mass [$M_\\odot$]", fontsize=14)
ax.set_title("Effect of LBV_flag on final primary mass", fontsize=15)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
