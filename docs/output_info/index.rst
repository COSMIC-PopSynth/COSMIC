.. _output_info:

###############################################################
Describing the output of COSMIC/BSE: Columns names/Values/Units
###############################################################

What do all the stupid numbers mean?
------------------------------------

Have you ever asked yourself "where are all the black holes"? 
Well, it's as easy as finding a system with kstar=14. 
If you want to learn about all types of stars in the Cosmos, you can
use the easily-accessible keys in the tables below. Shibby. 

The kstar specifies the evolutionary state of the star:

.. _kstar-table:

.. table:: Evolutionary State of the Star

    =====   ==================
    kstar   evolutionary state
    =====   ==================
    0       Main Sequence (MS), < 0.7[:math:`{\mathrm{M}_\odot}`]
    1       MS, > 0.7[:math:`{\mathrm{M}_\odot}`]
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

In addition, you may be intersted in the properties of the binary at
key moments in its evolution. These are also tracked by integer values
whose real meaning relates as follows:

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

This `pandas.DataFrame` tracks certain astrophysical information about your
binary at key moments in its evolution. 
Entries are added with changes in the :ref:`evolve-type-table`. 

=============  =====================================================
``tphys``      Evolution time [Myr]
``mass_1``     Primary mass [:math:`{\mathrm{M}_\odot}`]
``mass_2``     Secondary mass [:math:`{\mathrm{M}_\odot}`]
``kstar_1``    Evolutionary state of primary (see :ref:`kstar-table`)
``kstar_2``    Evolutionary state of secondary (see :ref:`kstar-table`)
``sep``        Semimajor axis [:math:`{\mathrm{R}_\odot}`]
``porb``       Orbital period [days]
``ecc``        Eccentricity
``RROL_1``     Stellar radius divided by Roche lobe radius for primary
``RROL_2``     Stellar radius divided by Roche lobe radius for secondary
``evol_type``  Key moment in evolution (see :ref:`evolve-type-table`)
``Vsys_1``     Change in systemic velocity due to first SN [km/s]
``Vsys_2``     Change in systemic velocity due to second SN [km/s]
``SNkick``     Magnitude of supernova natal kick [km/s]
``SNtheta``    Angular change in orbital plane due to supernova [degrees]
``aj_1``       Effective age of the primary [Myr]
``aj_2``       Effective age of the secondary [Myr]
``tms_1``      Main sequence lifetime [Myr]
``tms_2``      Main sequence lifetime [Myr]
``massc_1``    Core mass of primary [:math:`{\mathrm{M}_\odot}`]
``massc_2``    Core mass of secondary [:math:`{\mathrm{M}_\odot}`]
``rad_1``      Stellar radius of primary [:math:`{\mathrm{R}_\odot}`]
``rad_2``      Stellar radius of secondary [:math:`{\mathrm{R}_\odot}`]
``bin_num``    Unique binary index that is consistent across initial conditions, BCM and BPP DataFrames
=============  =====================================================



BCM DataFrame
-------------
This `pandas.DataFrame` provides detailed parameters of the binary at user-specified timesteps in the evolution.
By default, only the first and last timestep are output in the BCM DataFrame. 

=================  =====================================================
``tphys``          Evolution time [Myr]
``kstar_1``        Evolutionary state of primary (see :ref:`kstar-table`)
``mass0_1``        ??
``mass_1``         ??
``lumin_1``        ??
``rad_1``          ??
``teff_1``         ??
``massc_1``        ??
``radc_1``         ??
``menv_1``         ??
``renv_1``         ??
``epoch_1``        ??
``ospin_1``        ??
``deltam_1``       ??
``RROL_1``         ??
``kstar_2``        Evolutionary state of secondary (see :ref:`kstar-table`)
``mass0_2``        ??
``mass_2``         ??
``lumin_2``        ??
``rad_2``          ??
``teff_2``         ??
``massc_2``        ??
``radc_2``         ??
``menv_2``         ??
``renv_2``         ??
``epoch_2``        ??
``ospin_2``        ??
``deltam_2``       ??
``RROL_2``         ??
``porb``           ??
``sep``            ??
``ecc``            ??
``B_0_1``          ??
``B_0_2``          ??
``SNkick_1``       ??
``SNkick_2``       ??
``Vsys_final``     ??
``SNtheta_final``  ??
``SN_1``           ??
``SN_2``           ??
``bin_state``      ??
``merger_type``    ??
``bin_num``        An index that allows the user to track infomation about this
                   system across the initial conditions, BCM and BPP DataFrames
=================  =====================================================
