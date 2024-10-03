.. _output_info:

############################
Understanding COSMIC outputs
############################

Since COSMIC uses BSE as it's core binary evolution algorithm, the output of COSMIC follows most of the same conventions as BSE. The kstar values and evolution stages are nearly identical to their BSE counterparts.

``kstar`` - Evolutionary states of stars
----------------------------------------

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

``evol_type`` - Binary evolution events
----------------------------------------

The evolutionary changes of the binary are logged in the evol_type column, which is filled with integer values. The key for each integer is listed below:

.. _evolve-type-table:

.. table:: Evolve Type

    =========   =====================
    evol_type   evolutionary change
    =========   =====================
    1           initial state
    2           kstar change
    3           begin Roche lobe overflow
    4           end Roche lobe overflow
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
    100         RLOF interpolation timeout error
    =========   =====================


``bpp`` - Key evolutionary change table
---------------------------------------

This :class:`pandas.DataFrame` tracks a selection of binary parameters at key evolutionary changes.
Entries are added with changes in the :ref:`evolve-type-table`.
All values with a ``_1`` label refer to the primary (i.e., the initially more massive component); the bpp DataFrame also includes the same column for the secondary with ``_1`` replaced by ``_2``

================  =====================================================
``tphys``         Evolution time [:math:`{\rm{Myr}}`]
``mass_1``        Primary mass [:math:`{\mathrm{M}_\odot}`]
``kstar_1``       Evolutionary state of primary (see :ref:`kstar-table`)
``sep``           Semimajor axis [:math:`{\mathrm{R}_\odot}`]
``porb``          Orbital period [:math:`{\rm{days}}`]
``ecc``           Eccentricity
``RRLO_1``        Primary radius in units of Roche lobe radii
``evol_type``     Key moment in evolution (see :ref:`evolve-type-table`)
``aj_1``          Effective age of the primary [:math:`{\rm{Myr}}`]
``tms_1``         Primary main sequence lifetime [:math:`{\rm{Myr}}`]
``massc_1``       Primary core mass [:math:`{\mathrm{M}_\odot}`]
``rad_1``         Primary radius [:math:`{\mathrm{R}_\odot}`]
``mass0_1``       Previous epoch mass [:math:`{\mathrm{M}_\odot}`]
``lum_1``         Luminosity of the primary [:math:`{\mathrm{L}_\odot}`]
``teff_1``        Effective temperature of the primary [:math:`{\mathrm{K}}`]
``radc_1``        Primary core radius  [:math:`{\mathrm{R}_\odot}`]
``menv_1``        Mass of the envelope of the primary [:math:`{\mathrm{M}_\odot}`]
``renv_1``        Radius of the envelope of the primary [:math:`{\mathrm{R}_\odot}`]
``omega_spin_1``  Angular velocity of the primary [:math:`{\rm{yr}}^{-1}`]
``B_1``           Neutron star magnetic field [:math:`{\rm{G}}`]
``bacc_1``        (only for pulsars) :math:`\delta{\mathrm{M}_\odot}` during accretion, see Equation 7 in COSMIC paper
``tacc_1``        Accretion duration (used for magnetic field decay) [:math:`{\rm{Myr}}`]
``epoch_1``       Time spent in current evolutionary epoch [:math:`{\rm{Myr}}`]
``bhspin_1``      Black hole spin magnitude [dimensionless]
``bin_num``       Unique binary index that is consistent across initial conditions, bcm, bpp, and kick_info DataFrames
================  =====================================================



``bcm`` - User-specified timestep table
---------------------------------------

This :class:`pandas.DataFrame` provides several binary parameters at user-specified timesteps in the evolution.
By default, COSMIC saves only the first and last timestep in the bcm DataFrame.
All values with a ``_1`` label refer to the primary; the bcm DataFrame also includes the same column for the secondary with ``_1`` replaced by ``_2``

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
``omega_spin_1``   Primary spin [:math:`\rm{rad/yr}`]
``deltam_1``       Primary mass transfer rate [:math:`{\mathrm{M}_\odot/\rm{yr}}`]
``RRLO_1``         Primary radius in units of Roche lobe radii
``porb``           Orbital period [:math:`\rm{days}`]
``sep``            Semimajor axis [:math:`\mathrm{R}_{\odot}`]
``ecc``            Eccentricity
``B_1``            Neutron star magnetic field [:math:`{\rm{G}}`]
``SN_1``           Supernova type:

                    1: Iron core-collapse supernova

                    2: Electron capture supernova

                    3: Ultra-stripped supernova (these happen whenever a He-star undergoes a common envelope with a compact companion)

                    4: Accretion induced collapse supernova

                    5: Merger induced collapse

                    6: Pulsational pair instability

                    7: Pair instability supernova
``bin_state``      State of the binary: 0 [binary], 1 [merged], 2 [disrupted]
``merger_type``    String of the kstar's in the merger. For example, two neutron stars that merged will be '1313'. Set to '-001' if binary has not merged. 
``bin_num``        Unique binary index that is consistent across initial conditions, bcm. bpp, and kick_info DataFrames
=================  =====================================================

``kick_info`` - Table of natal kick information
-----------------------------------------------

kick_info is a (2,17) array that tracks information about supernova
kicks. This allows us to track the total change to the systemic
velocity and the total change in the orbital plane tilt after both
supernovae, as well as reproduce systems.
The first row contains information about the first supernova that
occurs, the second row the second supernova.
Note that some values the second row will take into account the
effect of the first SN (e.g., kick_info[2,10] is the total systemic
velocity after both supernovae).

==========================================  ========================================================================================================================================
``kick_info[i,1]: star``                    whether the exploding star is the primary (1) or secondary (2)
``kick_info[i,2]: disrupted``               whether the system was disrupted from the SN (0=no, 1=yes)
``kick_info[i,3]: natal_kick``              magnitude of the natal kick [:math:`{\rm{km/s}}`]
``kick_info[i,4]: phi``                     polar angle of explosion (in the frame of the exploding star) [:math:`{\rm{degrees}}`]
``kick_info[i,5]: theta``                   azimuthal angle of explosion (in the frame of the exploding star) [:math:`{\rm{degrees}}`]
``kick_info[i,6]: mean anomaly``            mean anomaly at time of explosion [:math:`{\rm{degrees}}`]
``kick_info[i,7]: delta_vsysx_1``           change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (x-component)
``kick_info[i,8]: delta_vsysy_1``           change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (y-component)
``kick_info[i,9]: delta_vsysz_1``           change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (z-component)
``kick_info[i,10]: vsys_1_total``           magnitude of systemic velocity of the binary if bound, or magnitude of total velocity of star=1 if disrupted, accounting for both SNe
``kick_info[i,11]: delta_vsysx_2``          change in 3D velocity of the star=2 if system is disrupted (x-component)
``kick_info[i,12]: delta_vsysy_2``          change in 3D velocity of the star=2 if system is disrupted (y-component)
``kick_info[i,13]: delta_vsysz_2``          change in 3D velocity of the star=2 if system is disrupted (z-component)
``kick_info[i,14]: vsys_2_total``           magnitude of velocity of star=2 if disrupted, accounting for both SNe [:math:`{\rm{km/s}}`]
``kick_info[i,15]: delta_theta_total``      angular change in orbital plane due to supernovae, relative to the pre-SN1 orbital plane [:math:`{\rm{degrees}}`]
``kick_info[i,16]: omega``                  azimuthal angle of the orbital plane w.r.t. spins [:math:`{\rm{degrees}}`]
``kick_info[i,17]: randomseed``             random seed at the start of call to kick.f

==========================================  ========================================================================================================================================
