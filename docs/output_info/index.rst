.. _output_info:

###############################################################
Describing the output of COSMIC/BSE: Columns names/Values/Units
###############################################################

Evolutionary states of stars/binaries
-------------------------------------

Since COSMIC uses BSE as it's core binary evolution algorithm, the output of COSMIC follows most of the same conventions as BSE. The kstar values and evolution stages are nearly identical to their BSE counterparts. 

The kstar value specifies the evolutionary state of the star:

.. _kstar-table:

.. table:: Evolutionary State of the Star

    =====   ==================
    kstar   evolutionary state
    =====   ==================
    0       Main Sequence (MS), < 0.7 :math:`{\mathrm{M}_\odot}`
    1       MS, > 0.7 :math:`{\mathrm{M}_\odot}`
    2       Hertzsprung Gap
    3       First Giant Branch
    4       Core Helium Burning
    5       Early Asymptotic Giant Branch (AGB)
    6       Thermally Pulsing AGB
    7       Naked Helium Star MS
    8       Naked Helium Star Hertzsprung Gap
    9       Naked Helium Star Giant Branch
    10      Helium White Dwarf
    11      Carbon/Oxygen White Dwarf
    12      Oxygen/Neon White Dwarf
    13      Neutron Star
    14      Black Hole
    15      Massless Remnant
    =====   ==================

The evolutionary changes of the binary are logged in the evol_type column, which is filled with integer values. The key for each integer is listed below:

.. _evolve-type-table:

.. table:: Evolve Type

    =========   =====================
    evol_type   evolutionary change
    =========   =====================
    1           initial state
    2           kstar change
    3           begin Roche lobe overflow
    4           end Roche lobe overlow
    5           contact
    6           coalescence
    7           begin common envelope
    8           end common envelope
    9           no remnant leftover
    10          max evolution time
    11          binary disruption
    12          begin symbiotic phase
    13          end symbiotic phase
    14          blue straggler
    15          supernova of primary
    16          supernova of secondary
    =========   =====================


BPP DataFrame
-------------

This `pandas.DataFrame` tracks a selection of binary parameters at key evolutionary changes. 
Entries are added with changes in the :ref:`evolve-type-table`. 

=============  =====================================================
``tphys``      Evolution time [:math:`{\rm{Myr}}`]
``mass_1``     Primary mass [:math:`{\mathrm{M}_\odot}`]
``mass_2``     Secondary mass [:math:`{\mathrm{M}_\odot}`]
``kstar_1``    Evolutionary state of primary (see :ref:`kstar-table`)
``kstar_2``    Evolutionary state of secondary (see :ref:`kstar-table`)
``sep``        Semimajor axis [:math:`{\mathrm{R}_\odot}`]
``porb``       Orbital period [:matH:`{\rm{days}}`]
``ecc``        Eccentricity
``RROL_1``     Primary radius divided by Roche lobe radius
``RROL_2``     Secondary radius divided by Roche lobe radius
``evol_type``  Key moment in evolution (see :ref:`evolve-type-table`)
``Vsys_1``     Change in systemic velocity due to first SN [:math:`{\rm{km/s}}`]
``Vsys_2``     Change in systemic velocity due to second SN [:math:`{\rm{km/s}}`]
``SNkick``     Magnitude of supernova natal kick [:math:`{\rm{km/s}}`]
``SNtheta``    Angular change in orbital plane due to supernova [:math:`{\rm{degrees}}`]
``aj_1``       Effective age of the primary [:math:`{\rm{Myr}}`]
``aj_2``       Effective age of the secondary [:math:`{\rm{Myr}}`]
``tms_1``      Primary main sequence lifetime [:math:`{\rm{Myr}}`]
``tms_2``      Secondary main sequence lifetime [:math:`{\rm{Myr}}`]
``massc_1``    Primary core mass [:math:`{\mathrm{M}_\odot}`]
``massc_2``    Secondary core mass [:math:`{\mathrm{M}_\odot}`]
``rad_1``      Primary radius [:math:`{\mathrm{R}_\odot}`]
``rad_2``      Secondary radius [:math:`{\mathrm{R}_\odot}`]
``bin_num``    Unique binary index that is consistent across initial conditions, BCM and BPP DataFrames
=============  =====================================================



BCM DataFrame
-------------
This `pandas.DataFrame` provides several binary parameters at user-specified timesteps in the evolution.
By default, COSMIC saves only the first and last timestep in the BCM DataFrame. 
All values with a `_1` label refer to the primary; the BCM DataFrame also includes the same column for the secondary with `_1` replaced by `_2`

=================  =====================================================
``tphys``          Evolution time [:math:`\rm{Myr}`]
``kstar_1``        Evolutionary state of primary (see :ref:`kstar-table`)
``mass0_1``        Previous evolutionary stage primary mass [:math:`{\mathrm{M}_\odot}`]
``mass_1``         Primary mass [:math:`{\mathrm{M}_\odot}`]
``lumin_1``        Primary luminosity [:math:`{\mathrm{L}_\odot}`]
``rad_1``          Primary radius [:math:`{\mathrm{R}_\odot}`]
``teff_1``         Primary effective temperature [:math:`{\rm{K}}`]
``massc_1``        Primary core mass [:math:`{\mathrm{M}_\odot}`]
``radc_1``         Primary core radius [:math:`{\mathrm{R}_\odot}`]
``menv_1``         Primary envelope mass [:math:`{\mathrm{M}_\odot}`]
``renv_1``         Primary envelope radius [:math:`{\mathrm{R}_\odot}`]
``epoch_1``        Primary epoch [:math:`\rm{Myr}`]
``ospin_1``        Primary spin [:math:`\rm{rad/yr}`] 
``deltam_1``       Primary mass transfer rate [:math:`{\mathrm{M}_\odot/\rm{yr}}`]
``RROL_1``         Primary radius divided by Roche lobe radius
``porb``           Orbital period [:math:`\rm{days}`]
``sep``            Semimajor axis [:math:`\mathrm{R}_{\odot}`]
``ecc``            Eccentricity
``B_0_1``          Initial neutron star magnetic field [:math:`{\rm{G}}`]
``SNkick_1``       Magnitude of first natal kick [:math:`{\rm{km/s}}`]
``Vsys_final``     Final systemic velocity magnitude [:math:`{\rm{km/s}}`]
``SNtheta_final``  Final systemic velocity angle [:math:`{\rm{degrees}}`]
``SN_1``           Supernova type:

                    1: Fe Core-collapse SN

                    2: Electron capture SN

                    3: Ultra-stripped supernovae (these happen whenever a He-star undergoes a CE with a compact companion)

                    4: Accretion induced collapse SN

                    5: Merger induced collapse

                    6: Pulsational-pair instability

                    7: Pair instability SN
``bin_state``      State of the binary: 0 [binary], 1 [merged], 2 [disrupted]
``merger_type``    String of the kstar's in the merger, '-001' if not merged
``bin_num``        Unique binary index that is consistent across initial conditions, BCM and BPP DataFrames 
=================  =====================================================
