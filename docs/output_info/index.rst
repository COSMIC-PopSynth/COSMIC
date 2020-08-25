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


bpp
---

This `pandas.DataFrame` tracks a selection of binary parameters at key evolutionary changes.
Entries are added with changes in the :ref:`evolve-type-table`.
All values with a `_1` label refer to the primary; the bpp DataFrame also includes the same column for the secondary with `_1` replaced by `_2`

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
``bacc_1``        (only for pulsars) :math:`\delta{\mathrm{M}_\odot}` during accretion see Equation 7 in COSMIC paper
``tacc_1``        Accretion duration (used for magnetic field decay) [:math:`{\rm{Myr}}`]
``epoch_1``
``bhspin_1``      Black hole spin magnitude [unitless]
``bin_num``       Unique binary index that is consistent across initial conditions, bcm and bpp DataFrames
================  =====================================================



bcm
---
This `pandas.DataFrame` provides several binary parameters at user-specified timesteps in the evolution.
By default, COSMIC saves only the first and last timestep in the bcm DataFrame.
All values with a `_1` label refer to the primary; the bcm DataFrame also includes the same column for the secondary with `_1` replaced by `_2`

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

                    1: Fe Core-collapse SN

                    2: Electron capture SN

                    3: Ultra-stripped supernovae (these happen whenever a He-star undergoes a CE with a compact companion)

                    4: Accretion induced collapse SN

                    5: Merger induced collapse

                    6: Pulsational-pair instability

                    7: Pair instability SN
``bin_state``      State of the binary: 0 [binary], 1 [merged], 2 [disrupted]
``merger_type``    String of the kstar's in the merger, '-001' if not merged
``bin_num``        Unique binary index that is consistent across initial conditions, bcm and bpp DataFrames
=================  =====================================================

kick_info
---------
kick_info is a (2,17) array that tracks information about the supernova
kicks. This allows us to track the total change to the systemic
velocity and the total change in the orbital plane tilt after both
supernovae, as well as reproduce systems.
The first row contains information about the first supernova that
occurs, the second row the second supernova.
Note that some values the second row will take into account the
effect of the first SN (e.g., kick_info[2,10] is the total systemic
velocity after both supernovae).

  * kick_info[i,1]: snstar of exploding star
  * kick_info[i,2]: disrupted (0=no, 1=yes)
  * kick_info[i,3]: magnitude of the natal kick
  * kick_info[i,4-5]: phi and theta (in the frame of the exploding star)
  * kick_info[i,6]: mean anomaly
  * kick_info[i,7-9]: change in 3D systemic velocity of the binary, or the change in 3D velocity of snstar=1 if the system is disrupted
  * kick_info[i,10]: magnitude of systemic velocity of the binary if bound or magnitude of total velocity of snstar=1 if disrupted, accounting for both SNe
  * kick_info[i,11-13]: change in 3D velocity of the snstar=2 if system is disrupted
  * kick_info[i,14]: magnitude of velocity of snstar=2 if disrupted,  accounting for both SNe
  * kick_info[i,15]: (total) tilt of the orbital plane after each SN w.r.t. the original angular momentum axis after each SN
  * kick_info[i,16]: azimuthal angle of the orbital plane w.r.t. spins
  * kick_info[i,17]: random seed at the start of call to kick.f

=====================  =====================================================
``star``               snstar of exploding star
``disrupted``          disrupted (0=no, 1=yes)
``natal_kick``         magnitude of the natal kick [:math:`{\rm{km/s}}`]
``phi``                of explosion (in the frame of the exploding star) [:math:`{\rm{degrees}}`]
``theta``              of explosion (in the frame of the exploding star) [:math:`{\rm{degrees}}`]
``mean anomaly``       mean anomaly [:math:`{\rm{degrees}}`]
``delta_vsysx_1``      change in 3D systemic velocity of the binary, or the change in 3D velocity of snstar=1 if the system is disrupted (x)
``delta_vsysy_1``      change in 3D systemic velocity of the binary, or the change in 3D velocity of snstar=1 if the system is disrupted (y)
``delta_vsysz_1``      change in 3D systemic velocity of the binary, or the change in 3D velocity of snstar=1 if the system is disrupted (z)
``vsys_1_total``       magnitude of systemic velocity of the binary if bound or magnitude of total velocity of snstar=1 if disrupted, accounting for both SNe
``delta_vsysx_2``      change in 3D velocity of the snstar=2 if system is disrupted (x)
``delta_vsysy_2``      change in 3D velocity of the snstar=2 if system is disrupted (y)
``delta_vsysz_2``      change in 3D velocity of the snstar=2 if system is disrupted (z)
``vsys_2_total``       magnitude of velocity of snstar=2 if disrupted,  accounting for both SNe [:math:`{\rm{km/s}}`]
``delta_theta_total``  Angular change in orbital plane due to supernova [:math:`{\rm{degrees}}`]
``omega``              azimuthal angle of the orbital plane w.r.t. spins [:math:`{\rm{degrees}}`]
``randomseed``         random seed at the start of call to kick.f

=====================  =====================================================
