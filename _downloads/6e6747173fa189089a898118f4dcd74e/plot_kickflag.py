"""
kickflag
========

This example demonstrates how changing the `kickflag` parameter can affect the natal kick velocity distributions of compact objects.
We separate this into neutron stars and black holes.
"""

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve
from cosmic.output import COSMICOutput


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

flag_labels = {
    1: r"kickflag=1: Hobbs+2005 (Maxwellian, $\sigma=265\ \rm km\ s^{-1}$)",
    2: r"kickflag=2: Giacobbo & Mapelli (2020) Eq. 1",
    3: r"kickflag=3: Giacobbo & Mapelli (2020) Eq. 2",
    4: r"kickflag=4: Bray & Eldridge (2016)",
    5: r"kickflag=5: Disberg & Mandel (2025)",
    6: r"kickflag=6: Mandel & Müller (2020)",
    7: r"kickflag=7: Janka (2017) (Eq. 2 in Chattaraj+2026)",
    8: r"kickflag=8: Richards+2023"
}

BINS = np.linspace(0, 800, 40)

# sample a population of neutron stars and black holes

ibt = InitialBinaryTable.sampler(
    'independent', [13, 14], [13, 14],
    binfrac_model=1.0, primary_model='kroupa01',
    ecc_model='sana12', porb_model='sana12',
    qmin=-1, SF_start=13700.0, SF_duration=0.0,
    met=0.002, size=2000
)[0]

ns_kick_results = {}
bh_kick_results = {}

# go through each kickflag
for flag in flag_labels:
    start_time = time.perf_counter()
    print(f"Evolving population with kickflag = {flag}...")

    # update BSEDict with new kickflag
    BSEDict['kickflag'] = flag

    # evolve a population and store output in results
    bpp, bcm, initC, kick_info = Evolve.evolve(
        initialbinarytable=ibt, 
        BSEDict=BSEDict,
        progress=True,
        nproc=4
    )
    results = COSMICOutput(bpp=bpp, bcm=bcm, initC=initC, kick_info=kick_info)

    # assign each row in the bpp a unique row_num
    results.bpp["row_num"] = np.arange(len(results.bpp))

    # find the rows that produce the first and second SNe
    first_pre_SN_rows = results.bpp[results.bpp["evol_type"] == 15]
    second_pre_SN_rows = results.bpp[results.bpp["evol_type"] == 16]

    # convert them to the post SN rows by increasing row_num by one to get the remnant type
    first_SN_type = results.bpp[results.bpp["row_num"].isin(first_pre_SN_rows["row_num"] + 1)]["kstar_1"].values
    second_SN_type = results.bpp[results.bpp["row_num"].isin(second_pre_SN_rows["row_num"] + 1)]["kstar_2"].values

    # split kicks by the first and second
    first_kick = results.kick_info[results.kick_info["star"] == 1]["natal_kick"].values
    second_kick = results.kick_info[results.kick_info["star"] == 2]["natal_kick"].values

    # use information on remnant type to separate into NSs and BHs
    ns_kick_results[flag] = np.concatenate([first_kick[first_SN_type == 13], second_kick[second_SN_type == 13]])
    bh_kick_results[flag] = np.concatenate([first_kick[first_SN_type == 14], second_kick[second_SN_type == 14]])
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time    
    print(f"   [Time taken: {elapsed_time:.4f} seconds]")

# plot the results for NSs and BHs separately
for kick_results, title in zip([ns_kick_results, bh_kick_results],
                               ["Neutron star kicks", "Black hole kicks"]):

    # separate into four flags at a time to clean up the plot
    for flags_to_plot, subtitle in zip([[1, 2, 3, 4], [5, 6, 7, 8]],
                                       ["kickflags 1-4", "kickflags 5-8"]):

        fig, ax = plt.subplots()

        for flag in flags_to_plot:
            ax.hist(kick_results[flag], bins=BINS, label=flag_labels[flag], density=True,
                    histtype="step", lw=2, color=f"C{flag - 1}")
        
        ax.set(
            xlabel=r"Natal kick, $v_{\rm kick}$ [$\rm km\ s^{-1}$]",
            ylabel=r"Probability density",
            title=f"{title} ({subtitle})"
        )

        ax.grid(True, linestyle=":", alpha=0.5, color="gray")

        ax.legend(
            loc="upper right", 
            fontsize=8.5, 
            facecolor='white', 
            edgecolor='none',
        )

        plt.tight_layout()

        plt.show()
