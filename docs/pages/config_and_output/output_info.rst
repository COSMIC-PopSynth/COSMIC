.. _output_info:

############################
Understanding COSMIC outputs
############################

Since COSMIC uses BSE as it's core binary evolution algorithm, the output of COSMIC follows most of the same conventions as BSE. The kstar values and evolution stages are nearly identical to their BSE counterparts.

You can learn more about how to analyze the output of COSMIC in the :ref:`analysis_interface` section.

``kstar`` - Evolutionary states of stars
----------------------------------------

The kstar value specifies the evolutionary state of the star:

.. _kstar-table:

.. table::
    :class: center-first-col

    ========================   ===========
    Stellar type (``kstar``)   Description
    ========================   ===========
    0                          Main Sequence (MS), < 0.7 :math:`{\mathrm{M}_\odot}`
    1                          MS, > 0.7 :math:`{\mathrm{M}_\odot}`
    2                          Hertzsprung Gap
    3                          First Giant Branch
    4                          Core Helium Burning
    5                          Early Asymptotic Giant Branch (AGB)
    6                          Thermally Pulsing AGB
    7                          Naked Helium Star MS
    8                          Naked Helium Star Hertzsprung Gap
    9                          Naked Helium Star Giant Branch
    10                         Helium White Dwarf
    11                         Carbon/Oxygen White Dwarf
    12                         Oxygen/Neon White Dwarf
    13                         Neutron Star
    14                         Black Hole
    15                         Massless Remnant
    ========================   ===========

``evol_type`` - Binary evolution events
----------------------------------------

The evolutionary changes of the binary are logged in the evol_type column, which is filled with integer values. The key for each integer is listed below:

.. _evolve-type-table:

.. table::
    :class: center-first-col

    ==================================  ================================
    Evolutionary event (``evol_type``)  Description
    ==================================  ================================
    1                                   Initial state
    2                                   Stellar type (``kstar``) change
    3                                   Begin Roche lobe overflow
    4                                   End Roche lobe overflow
    5                                   Contact
    6                                   Coalescence
    7                                   Begin common envelope
    8                                   End common envelope
    9                                   No remnant leftover
    10                                  Maximum evolution time
    11                                  Binary disruption
    12                                  Begin symbiotic phase
    13                                  End symbiotic phase
    14                                  Blue straggler
    15                                  Supernova of primary
    16                                  Supernova of secondary
    100                                 RLOF interpolation timeout error
    ==================================  ================================


``bpp`` - Key evolutionary change table
---------------------------------------

This :class:`pandas.DataFrame` tracks a selection of binary parameters at key evolutionary changes.
Entries are added with changes in the :ref:`evolve-type-table`.
All values with a ``_1`` label refer to the primary (i.e., the initially more massive component); the bpp DataFrame also includes the same column for the secondary with ``_1`` replaced by ``_2``

====================  =====================================================
Column name           Description
====================  =====================================================
``tphys``             Evolution time [:math:`{\rm{Myr}}`]
``mass_1``            Primary mass [:math:`{\mathrm{M}_\odot}`]
``kstar_1``           Evolutionary state of primary (see :ref:`kstar-table`)
``sep``               Semimajor axis [:math:`{\mathrm{R}_\odot}`]
``porb``              Orbital period [:math:`{\rm{days}}`]
``ecc``               Eccentricity
``RRLO_1``            Primary radius in units of Roche lobe radii
``evol_type``         Key moment in evolution (see :ref:`evolve-type-table`)
``aj_1``              Effective age of the primary [:math:`{\rm{Myr}}`]
``tms_1``             Primary main sequence lifetime [:math:`{\rm{Myr}}`]
``massc_he_layer_1``  Primary helium core mass [:math:`{\mathrm{M}_\odot}`]
``massc_co_layer_1``  Primary carbon-oxygen core mass [:math:`{\mathrm{M}_\odot}`]
``rad_1``             Primary radius [:math:`{\mathrm{R}_\odot}`]
``mass0_1``           Previous epoch mass [:math:`{\mathrm{M}_\odot}`]
``lum_1``             Luminosity of the primary [:math:`{\mathrm{L}_\odot}`]
``teff_1``            Effective temperature of the primary [:math:`{\mathrm{K}}`]
``radc_1``            Primary core radius  [:math:`{\mathrm{R}_\odot}`]
``menv_1``            Mass of the envelope of the primary [:math:`{\mathrm{M}_\odot}`]
``renv_1``            Radius of the envelope of the primary [:math:`{\mathrm{R}_\odot}`]
``omega_spin_1``      Angular velocity of the primary [:math:`{\rm{yr}}^{-1}`]
``B_1``               Neutron star magnetic field [:math:`{\rm{G}}`]
``bacc_1``            (only for pulsars) :math:`\delta{\mathrm{M}_\odot}` during accretion, see Equation 7 in COSMIC paper
``tacc_1``            Accretion duration (used for magnetic field decay) [:math:`{\rm{Myr}}`]
``epoch_1``           Time spent in current evolutionary epoch [:math:`{\rm{Myr}}`]
``bhspin_1``          Black hole spin magnitude [dimensionless]
``bin_num``           Unique binary index that is consistent across initial conditions, bcm, bpp, and kick_info DataFrames
====================  =====================================================



``bcm`` - User-specified timestep table
---------------------------------------

This :class:`pandas.DataFrame` provides several binary parameters at user-specified timesteps in the evolution.
By default, COSMIC saves only the first and last timestep in the bcm DataFrame.
All values with a ``_1`` label refer to the primary; the bcm DataFrame also includes the same column for the secondary with ``_1`` replaced by ``_2``

====================  =====================================================
Column name           Description
====================  =====================================================
``tphys``              Evolution time [:math:`\rm{Myr}`]
``kstar_1``            Evolutionary state of primary (see :ref:`kstar-table`)
``mass0_1``            Previous evolutionary stage primary mass [:math:`{\mathrm{M}_\odot}`]
``mass_1``             Primary mass [:math:`{\mathrm{M}_\odot}`]
``lumin_1``            Primary luminosity [:math:`{\mathrm{L}_\odot}`]
``rad_1``              Primary radius [:math:`{\mathrm{R}_\odot}`]
``teff_1``             Primary effective temperature [:math:`{\rm{K}}`]
``massc_he_layer_1``  Primary helium core mass [:math:`{\mathrm{M}_\odot}`]
``massc_co_layer_1``  Primary carbon-oxygen core mass [:math:`{\mathrm{M}_\odot}`]
``radc_1``             Primary core radius [:math:`{\mathrm{R}_\odot}`]
``menv_1``             Primary envelope mass [:math:`{\mathrm{M}_\odot}`]
``renv_1``             Primary envelope radius [:math:`{\mathrm{R}_\odot}`]
``epoch_1``            Primary epoch [:math:`\rm{Myr}`]
``omega_spin_1``       Primary spin [:math:`\rm{rad/yr}`]
``deltam_1``           Primary mass transfer rate [:math:`{\mathrm{M}_\odot/\rm{yr}}`]
``RRLO_1``             Primary radius in units of Roche lobe radii
``porb``               Orbital period [:math:`\rm{days}`]
``sep``                Semimajor axis [:math:`\mathrm{R}_{\odot}`]
``ecc``                Eccentricity
``B_1``                Neutron star magnetic field [:math:`{\rm{G}}`]
``SN_1``               Supernova type:
    
                        1: Iron core-collapse supernova
    
                        2: Electron capture supernova
    
                        3: Ultra-stripped supernova (these happen whenever a He-star undergoes a common envelope with a compact companion)
    
                        4: Accretion induced collapse supernova
    
                        5: Merger induced collapse
    
                        6: Pulsational pair instability
    
                        7: Pair instability supernova
``bin_state``          State of the binary: -1 [single], 0 [binary], 1 [merged], 2 [disrupted]
``merger_type``        String of the kstar's in the merger. For example, two neutron stars that merged will be '1313'. Set to '-001' if binary has not merged. 
``bin_num``            Unique binary index that is consistent across initial conditions, bcm. bpp, and kick_info DataFrames
====================  =====================================================

``kick_info`` - Table of natal kick information
-----------------------------------------------

kick_info is a (2,19) array that tracks information about supernova
kicks. This allows us to track the total change to the systemic
velocity and the total change in the orbital plane tilt after both
supernovae, as well as reproduce systems.
Each row contains information about the supernova of one star: the
first row is for the primary (star 1) and the second row for the
secondary (star 2), whichever of them explodes first (as for
``natal_kick_array``).
Some values in the row of the *second* supernova take into account
the effect of the first (e.g., ``vsys_1_total`` is then the total
systemic velocity after both supernovae), and its velocities are given
in the frame of the first supernova.

==========================================  ==================  ========================================================================================================================================
Index                                       Name                Description
==========================================  ==================  ========================================================================================================================================
``kick_info[i,1]``                          ``star``            whether the exploding star is the primary (1) or secondary (2)
``kick_info[i,2]``                          ``disrupted``       whether the system was disrupted from the SN (0=no, 1=yes)
``kick_info[i,3]``                          ``natal_kick``      magnitude of the natal kick [:math:`{\rm{km/s}}`]
``kick_info[i,4]``                          ``phi``             angle of the natal kick out of the pre-SN orbital plane, from -90 to 90 (see below) [:math:`{\rm{degrees}}`]
``kick_info[i,5]``                          ``theta``           azimuthal angle of the natal kick in the pre-SN orbital plane, measured from the x-axis (see below) [:math:`{\rm{degrees}}`]
``kick_info[i,6]``                          ``mean anomaly``    mean anomaly at time of explosion [:math:`{\rm{degrees}}`]
``kick_info[i,7]``                          ``delta_vsysx_1``   change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (x-component)
``kick_info[i,8]``                          ``delta_vsysy_1``   change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (y-component)
``kick_info[i,9]``                          ``delta_vsysz_1``   change in 3D systemic velocity of the binary, or the change in 3D velocity of star=1 if the system is disrupted (z-component)
``kick_info[i,10]``                         ``vsys_1_total``    magnitude of systemic velocity of the binary if bound, or magnitude of total velocity of star=1 if disrupted, accounting for both SNe
``kick_info[i,11]``                         ``delta_vsysx_2``   change in 3D velocity of the star=2 if system is disrupted (x-component)
``kick_info[i,12]``                         ``delta_vsysy_2``   change in 3D velocity of the star=2 if system is disrupted (y-component)
``kick_info[i,13]``                         ``delta_vsysz_2``   change in 3D velocity of the star=2 if system is disrupted (z-component)
``kick_info[i,14]``                         ``vsys_2_total``    magnitude of velocity of star=2 if disrupted, accounting for both SNe [:math:`{\rm{km/s}}`]
``kick_info[i,15]``                         ``theta_euler``     tilt of the orbital plane due to the SN: the angle between the orbital angular momentum before and after it, from 0 to 180 (see below) [:math:`{\rm{degrees}}`]
``kick_info[i,16]``                         ``phi_euler``       angle in the pre-SN orbital plane from the x-axis to the line of nodes, about which the orbit tilts (see below) [:math:`{\rm{degrees}}`]
``kick_info[i,17]``                         ``psi_euler``       angle in the post-SN orbital plane from the line of nodes to the new periapsis (see below) [:math:`{\rm{degrees}}`]
``kick_info[i,18]``                         ``randomseed``      random seed at the start of call to kick.f
``kick_info[i,19]``                         ``tphys``           time of the supernova in Myr

==========================================  ==================  ========================================================================================================================================

What do these kick angles mean?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Understanding the reference frame and angles of supernova kicks can be rather confusing (priv. comm. Wagg). Keep reading to gain a better understanding of the geometry of supernova kicks and how they are stored in COSMIC's ``kick_info`` array.

.. dropdown:: More detailed explanation of euler angles and geometry

    The result of a supernova kick of course depends on the kick magnitude and mass loss, but it's also important to know the orbital phase (:math:`M` (``mean_anomaly``)) and kick direction (:math:`\theta` (``theta``), :math:`\phi` (``phi``)). A kick in the direction of motion at periapsis is more likely to disrupt the binary than a kick of the same magnitude at apoapsis, and a kick out of the orbital plane can tilt the orbit without disrupting it. At the same time, it's important to account for the fact that the second supernova's kick is applied in a frame that has already been tilted by the first supernova. The Euler angles describe how the orbit is reoriented by each supernova, and are used to rotate the velocity changes from the second supernova into the frame of the first. ``theta_euler`` is the tilt of the orbit, the angle between the orbital angular momentum before (:math:`\hat{z}`) and after (:math:`\hat{h}`) the supernova. ``phi_euler`` is the angle, in the pre-SN orbital plane, from the :math:`x`-axis to the line of nodes :math:`\hat{z} \times \hat{h}` (where the old and new orbital planes intersect, and the axis about which the orbit tilts). ``psi_euler`` (:math:`\psi_E`) is the angle, in the post-SN orbital plane, from the line of nodes to the new periapsis of the exploding star's orbit relative to its companion.

    It's also worth mentioning that beyond the natal kick, the Blaauw kick resulting from symmetric mass loss can also significantly affect the binary's evolution. The exploding star of mass :math:`m_1` loses the ejecta :math:`\Delta m` and leaves a remnant of mass :math:`m_1 - \Delta m`. Without a natal kick, the binary is disrupted if it loses more than :math:`r/2a` of its total mass :math:`M_i = m_1 + m_2`, where :math:`r` is the separation at the supernova. For a circular orbit that is half of the total mass; an eccentric orbit is easier to disrupt near periapsis, where the stars move fastest, and harder near apoapsis.

If, like me, you find it hard to understand these concepts from wordy explanations, the interactive plots below are designed to help out! You can interactively explore how the kick angles and magnitudes affect the post-supernova orbit and the Euler angles that describe how the orbit is reoriented. The top view shows where the kick points relative to the orbit and what it does to the orbit's shape, and the 3D view shows how the orbit is reoriented and the Euler angles that COSMIC would store in ``kick_info``.

.. raw:: html

    <link rel="stylesheet" href="../../_static/kick_geometry/kick-geometry.css">

    <div class="kick-geometry">
      <div class="kg-explorer" data-kg-explorer>
        <figure>
          <div class="kg-fig-title">Interactive kick explorer</div>
          <p>The plot shows you a top-down view of the orbit of the binary. You can use the sliders below to adjust the pre-supernova system or change the natal kick and see how it changes the results (below).</p>
          <svg data-kg-svg="top" viewBox="0 0 560 410" role="img" aria-label="Top view of the pre-supernova orbit with the companion at the focus, the exploding star placed by the mean anomaly, its orbital velocity, the in-plane part of the kick at azimuth theta, and the orbit after the kick.">
            <defs>
              <marker id="kg-ah-kick" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" style="fill: var(--kg-kick)"/></marker>
              <marker id="kg-ah-vel" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" style="fill: var(--kg-vel)"/></marker>
              <marker id="kg-ah-ink" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0,0 L10,5 L0,10 z" style="fill: var(--kg-muted)"/></marker>
              <marker id="kg-ah-post" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" style="fill: var(--kg-post)"/></marker>
              <clipPath id="kg-clip-top"><rect x="0" y="0" width="560" height="410"/></clipPath>
              <clipPath id="kg-clip-3d"><rect x="0" y="0" width="560" height="410"/></clipPath>
            </defs>
            <g data-kg-layer="top"></g>
          </svg>
          <figcaption>The solid ellipse is the orbit before the supernova and the dashed curve the orbit after it. The velocity and kick arrows are to scale with one another (and the kick arrow shows only the in-plane projection).</figcaption>
        </figure>

        <div class="kg-controls">
          <div class="kg-sliders">
            <div class="kg-ctl"><label for="kg-ecc">Eccentricity e</label><output data-kg-out="ecc"></output><input id="kg-ecc" data-kg="ecc" type="range" min="0" max="0.9" step="0.01" value="0.5"></div>
            <div class="kg-ctl"><label for="kg-M">Mean anomaly M</label><output data-kg-out="M"></output><input id="kg-M" data-kg="M" type="range" min="0" max="360" step="1" value="60"></div>
            <div class="kg-ctl"><label for="kg-loga">Semi-major axis a</label><output data-kg-out="loga"></output><input id="kg-loga" data-kg="loga" type="range" min="0" max="3" step="0.01" value="1.7"></div>
            <div class="kg-ctl"><label for="kg-theta">Kick azimuth θ</label><output data-kg-out="theta"></output><input id="kg-theta" data-kg="theta" type="range" min="0" max="360" step="0.1" value="60"></div>
            <div class="kg-ctl"><label for="kg-phi">Kick latitude φ</label><output data-kg-out="phi"></output><input id="kg-phi" data-kg="phi" type="range" min="-90" max="90" step="1" value="30"></div>
            <div class="kg-ctl"><label for="kg-vk">Kick speed v<sub>k</sub></label><output data-kg-out="vk"></output><input id="kg-vk" data-kg="vk" type="range" min="0" max="500" step="1" value="100"></div>
            <div class="kg-ctl"><label for="kg-m1">Exploding star m<sub>1</sub></label><output data-kg-out="m1"></output><input id="kg-m1" data-kg="m1" type="range" min="1" max="50" step="0.1" value="12"></div>
            <div class="kg-ctl"><label for="kg-m2">Companion m<sub>2</sub></label><output data-kg-out="m2"></output><input id="kg-m2" data-kg="m2" type="range" min="0.1" max="50" step="0.1" value="8"></div>
            <div class="kg-ctl"><label for="kg-dm">Ejecta Δm</label><output data-kg-out="dm"></output><input id="kg-dm" data-kg="dm" type="range" min="0" max="11.9" step="0.1" value="2"></div>
          </div>
          <div class="kg-buttons" aria-label="Set θ">
            <button type="button" data-kg-set="M">θ = M</button>
            <button type="button" data-kg-set="nu">θ = ν (away from companion)</button>
            <button type="button" data-kg-set="vel">θ = v (along the motion)</button>
            <button type="button" data-kg-set="in">θ = in (toward companion)</button>
          </div>

        <details class="kg-readouts-wrap" open>
        <summary>Results</summary>
        <div class="kg-readouts">
          <dl class="kg-readout">
            <dt class="kg-readout-title kg-span">Immediately prior to supernova</dt>
            <dt>Companion orbital speed</dt><dd data-kg-out="r-v2"></dd>
            <dt>Eccentric anomaly E</dt><dd data-kg-out="r-E"></dd>
            <dt>True anomaly ν</dt><dd data-kg-out="r-nu"></dd>
          </dl>
          <dl class="kg-readout">
            <dt class="kg-readout-title">After the supernova</dt><dd class="kg-state" data-kg-out="r-state"></dd>
            <dt>Ratio of semi-major axes, a<sub>f</sub> / a<sub>i</sub></dt><dd data-kg-out="r-a"></dd>
            <dt>Eccentricity, e<sub>f</sub></dt><dd data-kg-out="r-e"></dd>
            <dt>Centre-of-mass speed</dt><dd data-kg-out="r-vcm"></dd>
            <dt>Remnant runaway speed</dt><dd data-kg-out="r-vsn"></dd>
            <dt>Companion runaway speed</dt><dd data-kg-out="r-vcomp"></dd>
          </dl>
        </div>
        </details>

        <figure>
          <div class="kg-fig-title">How the orbit is reoriented: Euler angles</div>
          <svg class="kg-3d" data-kg-svg="3d" tabindex="0" viewBox="0 0 560 410" role="img" aria-label="Rotatable 3D view of the orbital plane before the supernova and the tilted orbital plane after it, with the line of nodes, the new orbital angular momentum and new periapsis, and arcs marking the Euler angles theta_E, phi_E and psi_E.">
            <g data-kg-layer="3d"></g>
          </svg>
          <div class="kg-buttons" aria-label="3D view">
            <button type="button" data-kg-view="reset">Reset view</button>
            <button type="button" data-kg-view="top">Look down z</button>
            <button type="button" data-kg-view="nodes">Look along the line of nodes</button>
            <label class="kg-check"><input type="checkbox" data-kg="second"> Show the x-axis used at a second supernova</label>
          </div>
          <figcaption>Drag to rotate. The faint disc is the orbital plane before the supernova, with the original x, y and z axes, and the whiter disc is the plane after it. The companion sits at the origin and the open circle marks the exploding star at the moment of the supernova.</figcaption>
        </figure>



        <details class="kg-readouts-wrap" open>
        <summary>Euler angles</summary>
        <div class="kg-readouts">
          <dl class="kg-readout">
            <dt>theta_euler, θ<sub>E</sub></dt><dd data-kg-out="r-thetaE"></dd>
            <dt>phi_euler, φ<sub>E</sub></dt><dd data-kg-out="r-phiE"></dd>
            <dt>psi_euler, ψ<sub>E</sub></dt><dd data-kg-out="r-psiE"></dd>
          </dl>
          <p class="kg-note" data-kg-out="r-euler-note"></p>
        </div>
        </details>

        </div>
      </div>
    </div>

    <script src="../../_static/kick_geometry/kick-physics.js"></script>
    <script src="../../_static/kick_geometry/kick-geometry.js"></script>
