.. _inifile:

##############################
Configuration files for COSMIC
##############################

How to write a configuration file
=================================

The `cosmic-pop` command-line executable cannot run without a configuration file.
Each of the sections below lists inputs for COSMIC's modified version of BSE. Each input has a description of the allowed values and an example of what that section might look like in the INI format.

Here is a link to the most recent stable release version of the default
inifile for COSMIC: `Stable Version INIFILE <https://github.com/COSMIC-PopSynth/COSMIC/blob/master/examples/Params.ini>`_

Here is a link to the unstable development version of the default inifile for COSMIC: `Development Version INIFILE <https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/examples/Params.ini>`_

filters
-------

=======================  ===============================================================
``select_final_state``   The bcm array generally returns the first and last
                         state of the binary system. Since we already
                         save the initial conditions in cosmic-pop, usually
                         we opt to throw out the initial state of the binary
                         from this array and only keep the final state in the
                         bcm array. The allowed values are:

                            ``True`` : bcm array only contains final state of binary

                            ``False`` : bcm array contains multiple binary states 
``binary_state``         Each binary system will end its evolution in one of
                         three states

                            ``0`` : the system is still a binary at the end its 
                            evolution 

                            ``1`` : the system merged before the end of its 
                            evolution

                            ``2`` : the system was disrupted before the end of 
                            its evolution

=======================  ===============================================================

.. code-block:: ini

    [filters]
    ; select_final_state will retain only the final entry of the bcm arry if set to True
    ; default=True
    select_final_state = True

    ; binary_state determines which types of binaries endstates to retain
    ; 0 alive today, 1 merged, 2 disrupted
    ; default=[0,1]
    binary_state = [0,1]

sampling
--------

=======================  ===================================================================================
``sampling_method``      Select which models to use to generate an initial 
                         sample of binary parameters at Zero Age Main Sequence

                            ``independent`` : initialize binaries with 
                            independent parameter distributions for the primary
                            mass, mass ratio, eccentricity, separation, and 
                            binary fraction

                            ``multidim`` : initialize binaries with 
                            multidimensional parameter distributions according to 
                            `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_
``galaxy_component``     There are four star formation histories implemented in 
                         COSMIC, one each for the Milky Way thin disk, thick 
                         disk, and bulge; and one that assumes a delta-function
                         burst of star formation. Their assumptions are as follows:

                            ``ThinDisk`` : constant star formation over 10 Gyr

                            ``ThickDisk`` : 1 Gyr burst of constant star formation
                            11 Gyr in the past

                            ``Bulge`` : 1 Gyr burst of constant star formation
                            10 Gyr in the past

                            ``DeltaBurst`` : single burst of star formation 13.7
                            Gyr in the past

``metallicity``          Single value for the metallicity of the population
                         where solar metallicity is Z = 0.02

=======================  ===================================================================================

.. code-block:: ini

    [sampling]
    ; Specify if you woud like to sample initial conditions via
    ; the independent method (indepedent) or would like to sample
    ; initial conditions follow Moe & Di Stefano (2017) (multidim)
    sampling_method=multidim

    ; Galaxy Components. Options include Bulge, ThinDisk, ThickDisk, and DeltaBurst
    ; Think Star Formation History with which to use as your basis for sampling
    galaxy_component = DeltaBurst

    ; Metallicity of the population of initial binaries
    metallicity = 0.02

[convergence]
-------------

===========================  ===================================================================================
``convergence_params``  
                             A list of parameters you would like to verify have converged
                             to a single distribution shape.
                             Options include: ``mass_1``, ``mass_2``, ``sep``, ``porb``,
                             ``ecc``, ``massc_1``, ``massc_2``, ``rad_1``, ``rad_2``

``convergence_limits``       Specifies limits for parameters included in the convergence
                             params list. If specified, the lower and upper limit must 
                             be specified:  

                                ``convergence_limits = {'mass_1' : [5, 10], 'sep' : [0, 10]}``

``convergence_filter``       Selects the stage of the evolution at which you would like
                             to check for convergence. This will filter for systems that
                             satisfy the final_kstar1 and final_kstar2 selections from
                             the command line call of cosmic-pop

                                ``formation``: computes convergence on binary properties
                                at formation with user-specified final kstars

                                ``1_SN``: computes convergence on binary properties
                                just before the first supernova for the population with
                                user-specified final kstars

                                ``2_SN``: computes convergence on binary properties
                                just before the second supernova for the population with
                                user-specified final kstars

                                ``disruption``: computes convergence on binary properties
                                just before disruption of the population with
                                user-specified final kstars

                                ``final_state``: computes convergence on binary properties
                                after the full evolution specified by the user-supplied evolution time
                                and with the user specified final kstars

                                ``XRB_form``: computes convergence on binary properties
                                at the start of RLO following the first supernova on the population with
                                user-specified final kstars

``match``                    match provides the tolerance for the convergence calculation
                             and is calculated as match = Log\ :sub:`10` (1-convergence)

``bcm_bpp_initCond_filter``  The bcm_bpp_initCond_filter will filter the bcm, bpp, and initCond 
                             DataFrames to only contain the binaries that satisfy the conv DataFrame.

                                ``True``: bcm, bpp, initCond will contain only the binaries which 
                                are in the conv DataFrame

                                ``False``: bcm, bpp, initCond will contain all systems which satisfy the                                 final kstar selection and will **not** be filtered.

                             **bcm_bpp_initCond_filter=False**

===========================  ===================================================================================

.. code-block:: ini

    [convergence]
    ; A list of parameters you would like to verify have converged
    ; to a single distribution shape.
    ; Options include mass_1, mass_2, sep, porb, ecc, massc_1, massc_2
    ; rad_1, rad_2
    convergence_params = [mass_1,mass_2,porb,ecc]

    ; convergence_limits is a dictionary that can contain limits for convergence params
    ; convergence_limits = {"mass_1" : [0, 20], "sep" : [0,5000]}
    convergence_limits = {}

    ; formation computes convergence on binary properties
    ; at formation with user-specified final kstars 

    ; 1_SN computes convergence on binary properties
    ; just before the first supernova for the population with
    ; user-specified final kstars

    ; 2_SN computes convergence on binary properties
    ; just before the second supernova for the population with
    ; user-specified final kstars 

    ; disruption computes convergence on binary properties
    ; just before disruption of the population with
    ; user-specified final kstars 

    ; final_state computes convergence on binary properties
    ; after the full evolution specified by the user-supplied evolution time
    ; and with the user specified final kstars

    ; XRB_form computes convergence on binary properties
    ; at the start of RLO following the first supernova on the population with
    ; user-specified final kstars
    convergence_filter = formation

    ; match provides the tolerance for the convergence calculation
    ; and is calculated as match = log10(1-convergence)
    ; default = -5.0
    match = -5.0

    ; bcm_bpp_initCond_filter filters the bcm/bpp/initCond
    ; arrays to only the binaries that are in the conv array if true
    bcm_bpp_initCond_filter = False

[rand_seed]
-----------

====================  ========================================================
``rand_seed``         Seed used to seed numpy.random.seed
====================  ========================================================

.. code-block:: ini

    [rand_seed]
    ; random seed int
    seed = 42

[bse]
-----

.. note::

    Although this is all one section, we have grouped the
    flags/parameters which get passed to the binary stellar evolution
    code into types. Each group will start with a note to indicate
    the type of parameter or flag.

.. note::

    SAMPLING FLAGS

=======================  =====================================================
``pts1``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         Main Sequence (MS) stars 

                         **pts1 = 0.001** following `Bannerjee+2019 <https://ui.adsabs.harvard.edu/abs/2019arXiv190207718B/abstract>`_
                         
``pts2``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         Giant Branch (GB, CHeB, AGB, HeGB) stars

                         **pts2 = 0.01** following `Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_
``pts3``                 determines the timesteps chosen in each evolution phase as
                         decimal fractions of the time taken in that phase for
                         HG, HeMS stars

                         **pts3 = 0.02** following `Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;
    ;;; SAMPLING FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;

    ; dtp is the timestep (in Myr) for outputting to the bcm array
    ; if dtp=0, will print every timestep (not recommended)
    ; if not set, it will automatically set to dtp=tphsyf (default)
    ;dtp = 1.0

    ; pts1,pts2,pts3 determine the timesteps chosen in each
    ;                 pts1 - MS                  (default=0.001, see Banerjee+ 2019)
    pts1=0.001
    ;                 pts2 - GB, CHeB, AGB, HeGB (default=0.01)
    pts2=0.01
    ;                 pts3 - HG, HeMS            (default=0.02)
    pts3=0.02

.. note::

    WIND FLAGS

=======================  =====================================================
``windflag``             Selects the model for wind mass loss for each star

                            ``0`` : Standard SSE/BSE (`Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_)

                            ``1`` : StarTrack (`Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_)

                            ``2`` : Metallicity dependence for O/B stars and Wolf Rayet stars (`Vink+2001 <http://adsabs.harvard.edu/abs/2001A&amp;A...369..574V>`_, `Vink+2005 <https://ui.adsabs.harvard.edu/abs/2005A%26A...442..587V/abstract>`_)

                            ``3`` : Same as 2, but LBV-like mass loss for giants 
                            and non-degenerate stars beyond the 
                            Humphreys-Davidson limit 

                         **windflag = 3**
``eddlimflag``           Limits the mass-loss rate of low-metallicity stars near 
                         the Eddington limit
                         (see `Grafener+2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...535A..56G/abstract>`_, `Giacobbo+2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.2959G/abstract>`_).

                            ``0`` : does not apply Eddington limit

                            ``1`` : applies Eddington limit

                         **eddlimflag = 0**
``neta``                 Reimers mass-loss coefficent (`Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_).
                         Note: this equation has a typo. There is an extra 
                         :math:`{\eta}` out front; the correct rate is directly proportional
                         to :math:`{\eta}`.
                         See also `Kurdritzki+1978, Section Vb <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1978A%26A....70..227K&link_type=ARTICLE&db_key=AST&high=#page=12>`_ for discussion.

                            ``positive value`` : supplies :math:`{\eta}` to `Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_

                         **neta = 0.5**
``bwind``                Binary enhanced mass loss parameter. 
                         See `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``positive value`` : supplies B\ :sub:`w` to `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_

                         **bwind = 0, inactive for single**
``hewind``               Helium star mass loss parameter: 10\ :sup:`-13` *hewind* L\ :sup:`2/3` gives He star mass-loss. Equivalent to 1 - :math:`{\mu}` in the last equation on `page 19 of SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_.

                         **hewind = 0.5**
``beta``                 Wind velocity factor: v\ :sub:`wind` :sup:`2` goes like *beta*. See `Equation 9 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``negative value`` : StarTrack (`Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_)

                            ``positive value`` : supplies :math:`{\beta}`\ :sub:`w` to `Equation 9 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_

                         **beta = 0.125**
``xi``                   Wind accretion efficiency factor, which gives the fraction 
                         of angular momentum lost via winds from the primary that 
                         transfers to the spin angular momentum of the companion. 
                         Corresponds to :math:`{\mu}`\ :sub:`w` in `Equation 11 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                            ``positive value`` : supplies :math:`{\mu}`\ :sub:`w` in `Equation 11 of Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_ 

                         **xi = 1.0**
``acc2``                 Bondi-Hoyle wind accretion factor where the mean wind accretion rate onto the secondary is proportional to *acc2*. See `Equation 6 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_.

                            ``positive value`` : supplies :math:`{\alpha}`\ :sub:`w` in `Equation 6 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_

                         **acc2 = 1.5**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;
    ;;; WIND FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;

    ; windflag sets the wind prescription
    ; windflag=0: stock BSE; windflag=1: StarTrack 2008
    ; windflag=2: Vink+2001; windflag=3: Vink+2005 (Vink plus LBV winds)
    ; default=3
    windflag=3

    ; neta is the Reimers mass-loss coefficent
    ; for more information, see Kudritzki & Reimers 1978, A&A 70, 227
    ; default=0.5
    neta = 0.5

    ; bwind is the binary enhanced mass loss parameter
    ; bwind it is always inactive for single stars
    ; default=0.0
    bwind = 0.0

    ; hewind is a helium star mass loss factor, between 0 and 1
    ; only applies if windflag=0, otherwise it is overwritten
    ; default=0.5
    hewind = 0.5

    ; beta is wind velocity factor: proportional to vwind^2
    ; beta<0: follows StarTrack 2008; beta=0.125: stock BSE
    ; default=0.125
    beta=0.125

    ; xi is the wind accretion efficiency factor, which gives the fraction of angular momentum lost via winds from the primary that transfers to the spin angular momentum of the companion
    ; default=1.0
    xi=1.0

    ; acc2 sets the Bondi-Hoyle wind accretion factor onto companion
    ; default=1.5
    acc2=1.5

.. note::

    COMMON ENVELOPE FLAGS

=======================  =====================================================
``alpha1``               Common-envelope efficiency parameter which scales the 
                         efficiency of transferring orbital energy to the 
                         envelope. See `Equation 71 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                            ``positive values`` : supplies :math:`{\alpha}` to `Equation 71 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_

                         **alpha1 = 1.0**
``lambdaf``              Binding energy factor for common envelope evolution. 
                         The initial binding energy of the stellar envelope 
                         goes like 1 / :math:`{\lambda}`. See `Equation 69 in Hurley+2002 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                            ``positive values`` : uses variable lambda prescription detailed 
                            in appendix of `Claeys+2014 <https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract>`_

                            ``negative values`` : fixes :math:`{\lambda}` to a value of -1.0* *lambdaf*

                         **lambdaf = 0.5**

``ceflag``               Selects the `de Kool 1990 <https://ui.adsabs.harvard.edu/abs/1990ApJ...358..189D/abstract>`_ 
                         model to set the initial orbital energy using the 
                         total mass of the stars instead of the core masses as 
                         in `Equation 70 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_.

                            ``0`` : Uses the core mass to calculate initial 
                            orbital energy as
                            in `Equation 70 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``1`` : Uses the `de Kool 1990 <https://ui.adsabs.harvard.edu/abs/1990ApJ...358..189D/abstract>`_ 
                            model

                         **ceflag = 0** 
``cekickflag``           Selects which mass and separation values to use when 
                         a supernova occurs during the CE and a kick
                         needs to be applied.

                            ``0`` : uses pre-CE mass and post-CE sep (BSE default)

                            ``1`` : uses pre-CE mass and sep values

                            ``2`` : uses post-CE mass and sep

                         **cekickflag = 2**

``cemergeflag``          Determines whether stars that begin a CE 
                         without a core-envelope boundary automatically lead to 
                         merger in CE. These systems include: 
                         kstars = [0,1,2,7,8,10,11,12].

                            ``0`` : allows the CE to proceed 

                            ``1`` : causes these systems to merge in the CE

                         **cemergeflag = 0**
``cehestarflag``         Uses fitting formulae from `Tauris+2015 <https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2123T/abstract>`_
                         for evolving RLO systems with a helium star donor 
                         and compact object accretor.
                         NOTE: this flag will override choice made by 
                         cekickflag if set

                            ``0`` : does NOT use Tauris+2015 at all

                            ``1`` : uses Tauris+2015 fits for final period only

                            ``2`` : uses Tauris+2015 fits for both final mass and final period

                         **cehestarflag = 0**
``qcflag``               Selects model to determine critical mass ratios for the
                         onset of unstable mass transfer and/or a common envelope
                         during RLO. 
                         NOTE: this is overridden by qcrit_array if any of the 
                         values are non-zero.

                            ``0`` : follows `Section 2.6 of Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                            (Default BSE)

                            ``1`` : same as 0 but with `Hjellming & Webbink 1987 <https://ui.adsabs.harvard.edu/abs/1987ApJ...318..794H/abstract>`_
                            for GB/AGB stars

                            ``2`` : follows `Table 2 of Claeys+2014 <https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract>`_

                            ``3`` : same as 2 but with `Hjellming & Webbink 1987 <https://ui.adsabs.harvard.edu/abs/1987ApJ...318..794H/abstract>`_
                            for GB/AGB stars

                         **qcflag = 2**

``qcrit_array``          Array with length: 16 for user-input values for the 
                         critical mass ratios that govern the onset of unstable
                         mass transfer and a common envelope. Each item is set 
                         individually for its associated kstar, and a value of 
                         0.0 will apply prescription of the qcflag for that kstar.
                         
                         **qcrit_array = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; COMMON ENVELOPE FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; alpha1 is the common-envelope efficiency parameter
    ; default=1.0
    alpha1 = 1.0

    ; lambdaf is the binding energy factor for common envelope evolution
    ; lambdaf>0.0 uses variable lambda prescription in appendix of Claeys+2014
    ; lambdaf<0 uses fixes lambda to a value of -1.0*lambdaf
    ; default=0.5
    lambdaf = 0.5

    ; ceflag=1 used the method from de Kool 1990 for setting the initial orbital energy
    ; ceflag=0 does not use this method (uses the core mass to calculate initial orbital energy)
    ; default=0
    ceflag=0

    ; cekickflag determined the prescription for calling kick.f in comenv.f
    ; 0: default BSE
    ; 1: uses pre-CE mass and sep values
    ; 2: uses post-CE mass and sep
    ; default=0
    cekickflag=0

    ; cemergeflag determines whether stars without a core-envelope boundary automatically lead to merger in CE
    ; cemergeflag=1 turns this on (causes these systems to merge)
    ; default=0
    cemergeflag=0

    ; cehestarflag uses fitting formulae from TLP, 2015, MNRAS, 451 for evolving RLO systems with a helium star donor and compact object accretor
    ; this flag will override choice made by cekickflag if set
    ; 0: off
    ; 1: fits for final period only
    ; 2: fits for both final mass and final period
    ; default=0
    cehestarflag=0

    ; qcflag is an integer flag that sets the model to determine which critical mass ratios to use for the onset of unstable mass transfer and/or a common envelope. NOTE: this is overridden by qcrit_array if any of the values are non-zero.
    ; 0: standard BSE
    ; 1: BSE but with Hjellming & Webbink, 1987, ApJ, 318, 794 GB/AGB stars
    ; 2: following binary_c from Claeys+2014 Table 2
    ; 3: following binary_c from Claeys+2014 Table 2 but with Hjellming & Webbink, 1987, ApJ, 318, 794 GB/AGB stars
    ; default=3
    qcflag=3

    ; qcrit_array is a 16-length array for user-input values for the critical mass ratios that govern the onset of unstable mass transfer and a common envelope
    ; each item is set individually for its associated kstar, and a value of 0.0 will apply prescription of the qcflag for that kstar
    ; default: [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    qcrit_array=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

.. note::

    KICK FLAGS

=======================  =====================================================
``sigma``                Sets the dispersion in the Maxwellian for the 
                         SN kick velocity in km/s

                            ``positive value`` : sets Maxwellian dispersion 

                         **default=265.0**
``bhflag``               Sets the model for how SN kicks are applied to BHs 
                         where bhflag != 0 allows velocity kick at BH formation
                    
                            ``0`` : no BH kicks

                            ``1`` : fallback-modulated kicks following
                            `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_

                            ``2`` : kicks decreased by ratio of BH mass to NS mass
                            (1.44 Msun); conserves linear momentum

                            ``3`` : full strength kick drawn from Maxwellian
                            with dispersion = *sigma* selected above

                         **bhflag = 1**
``ecsn``                 Allows for electron capture SN and sets the 
                         maximum ECSN mass range at the time of SN

                            ``0`` : turns off ECSN

                            ``positive values`` : `BSE (Hurley+2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_ 
                            and `StarTrack (Belczynski+2008) <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_ 
                            use ecsn = 2.25, while `Podsiadlowksi+2004 <https://ui.adsabs.harvard.edu/abs/2004ApJ...612.1044P/abstract>`_ 
                            use ecsn = 2.5

                         **ecsn = 2.5**
``ecsn_mlow``            Sets the low end of the ECSN mass range

                            ``positive values`` : `BSE (Hurley+2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                            use ecsn_mlow = 1.6, while `StarTrack (Belczynski+2008) <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_ 
                            use ecsn_mlow = 1.85, while `Podsiadlowksi+2004 <https://ui.adsabs.harvard.edu/abs/2004ApJ...612.1044P/abstract>`_
                            use ecsn_mlow = 1.4

                         **ecsn_mlow = 1.4**
``sigmadiv``             Sets the modified ECSN kick strength
                         
                         ``positive values`` : divide *sigma* above by *sigmadiv*

                         ``negative values`` : sets the ECSN *sigma* value

                         **sigmadiv = -20.0**
``aic``                  reduces kick strengths for accretion induced collapse SN
                         according to *sigmadiv*

                            ``0`` : AIC SN receive kicks drawn from Maxwellian 
                            with dispersion = *sigma* above

                            ``1`` : sets kick strength according to *sigmadiv* 
                            NOTE: this will applies even if ecsn = 0.0

                         **aic = 1**
``ussn``                 Reduces kicks according to the *sigmadiv* selection
                         for ultra-stripped supernovae which happen whenever
                         a He-star undergoes a CE with a compact companion

                            ``0`` : USSN receive kicks drawn from Maxwellian
                            with dispersion = *sigma* above

                            ``1`` : sets kick strength according to *sigmadiv*

                         **ussn = 0**
``pisn``                 Allows for (pulsational) pair instability supernovae
                         and sets either the model to use or the maximum mass 
                         of the remnant.

                            ``0`` : no pulsational pair instability SN

                            ``-1`` : uses the formulae from `Spera & Mapelli 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.470.4739S/abstract>`_

                            ``-2`` : uses a polynomial fit to `Table 1 in Marchant+2018 <https://ui.adsabs.harvard.edu/abs/2018arXiv181013412M/abstract>`_

                            ``-3`` : uses a polynomial fit to `Table 5 in Woosley 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...878...49W/abstract>`_

                            ``positive values`` : turns on no pulsational pair 
                            instability SN and sets the maximum mass of the allowed
                            remnant

                         **pisn = 45.0**
``bhsigmafrac``          Sets a fractional modification which scales down *sigma*
                         for BHs. This works in addition to whatever is chosen for 
                         *bhflag*, and is applied to *sigma* **before** the *bhflag*
                         prescriptions are applied

                            ``values between [0, 1]`` : reduces *sigma* by *bhsigmafrac*   

                         **bhsigmafrac = 1.0**
``polar_kick_angle``     Sets the opening angle of the SN kick relative to the 
                         pole of the exploding star, where 0 gives strictly polar
                         kicks and 90 gives fully isotropic kicks

                            ``values between [0, 90]`` : sets opening angle for SN kick

                         **polar_kick_angle = 90.0**
``natal_kick_array``     Array of lenght: 6 which takes user input values to fix
                         the SN natal kick, where the array is 
                         formatted as: [vk1, vk2, phi1, phi2, theta1, theta2].
                         NOTE: any numbers outside these ranges will be sampled
                         in the standard ways detailed above.

                            ``vk1, vk2`` : valid on the range [0, inf] 

                            ``phi1, phi2`` : co-lateral polar angles valid from 
                            [-pi/2, pi/2]

                            ``theta1, theta2`` : azimuthal angles valid from 
                            [0, 2pi]

                         **natal_kick_array = [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;
    ;;; KICK FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;

    ; sigma sets is the dispersion in the Maxwellian for the SN kick velocity in km/s
    ; default=265.0
    sigma=265.0

    ; bhflag != 0 allows velocity kick at BH formation
    ; bhflag=0: no BH kicks; bhflag=1: fallback-modulated kicks
    ; bhflag=2: mass-weighted (proportional) kicks; bhflag=3: full NS kicks
    ; default=1
    bhflag=1

    ; ecsn>0 turns on ECSN and also sets the maximum ECSN mass range (at the time of the SN)
    ; stock BSE and StarTrack: ecsn=2.25; Podsiadlowski+2004: ecsn=2.5)
    ; default=2.5
    ecsn=2.5

    ; ecsn_mlow sets the low end of the ECSN mass range
    ; stock BSE:1.6; StarTrack:1.85; Podsiadlowski+2004:1.4)
    ; default=1.4
    ecsn_mlow=1.4

    ; sigmadiv sets the modified ECSN kick
    ; negative values sets the ECSN sigma value, positive values divide sigma above by sigmadiv
    ; default=-20.0
    sigmadiv=-20.0

    ; aic=1 turns on low kicks for accretion induced collapse
    ; works even if ecsn=0
    ; default=1
    aic=1

    ; ussn=1 uses reduced kicks (drawn from the sigmadiv distritbuion) for ultra-stripped supernovae
    ; these happen whenever a He-star undergoes a CE with a compact companion
    ; default=0
    ussn=1

    ; pisn>0 allows for (pulsational) pair instability supernovae
    ; and sets the maximum mass of the remnant
    ; pisn=-1 uses the formulae from Spera+Mapelli 2017 for the mass
    ; pisn=0 turns off (pulsational) pair instability supernovae
    ; default=45
    pisn=45.0

    ; bhsigmafrac sets the fractional modification used for scaling down the sigma for BHs
    ; this works in addition to whatever is chosen for bhflag, and is applied to the sigma beforehand these prescriptions are implemented
    ; default=1.0
    bhsigmafrac = 1.0

    ; polar_kick_angle sets the opening angle of the kick relative to the pole of the exploding star
    ; this can range from 0 (strictly polar kicks) to 90 (fully isotropic kicks)
    ; default=90.0
    polar_kick_angle = 90.0

    ; natal_kick_array is a 6-length array for user-input values for the SN natal kick
    ; formatted as: (vk1, vk2, phi1, phi2, theta1, theta2)
    ; vk is valid on the range [0, inf], phi are the co-lateral polar angles valid from [-pi/2, pi/2], and theta are azimuthal angles [0, 2*pi]
    ; any number outside of these ranges will be sampled in the standard way in kick.f
    ; default=[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]
    natal_kick_array=[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]

.. note::

    REMNANT MASS FLAGS

=======================  =====================================================
``nsflag``               Determines the remnant mass prescription used for
                         NSs and BHs.

                            ``0`` : follows `Section 6 of Hurley+2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_
                            (default BSE)

                            ``1`` : follows `Belczynski+2002 <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..407B/abstract>`_

                            ``2`` : follows `Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                            ``3`` : follows the rapid prescription from `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_

                            ``4`` : delayed prescription from `Fryer+2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...749...91F/abstract>`_

                         **nsflag = 3**
``mxns``                 Sets the boundary between the maximum NS mass
                         and the minimum BH mass

                            ``positive values`` : sets the NS/BH mass bounary

                         **mxns = 2.5**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; REMNANT MASS FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; nsflag determines the remnant mass prescription used
    ; nsflag=0: default BSE; nsflag=1: Belczynski et al. 2002, ApJ, 572, 407
    ; nsflag=2: Belczynski et al. 2008; nsflag=3: rapid prescription (Fryer+ 2012)
    ; nsflag=4: delayed prescription (Fryer+ 2012)
    ; default=3
    nsflag=3

    ; mxns sets the maximum NS mass
    ; default=2.5
    mxns=2.5

.. note::

    REMNANT SPIN FLAGS

=======================  ===============================================================
``bhspinflag``           Uses different prescriptions for BH spin after formation 

                            ``0`` : sets all BH spins to *bhspinmag*

                            ``1`` : draws a random BH spin between 0 and bhspinmag for every BH

                            ``2`` : core-mass dependent BH spin (based on `Belczynski+2017 v1 <https://arxiv.org/abs/1706.07053v1>`_)
                         **bhspinflag = 0**

``bhspinmag``            Sets either the spin of all BHs or the upper limit of the uniform distribution for BH spins 

                            ``values >= 0.0`` : spin or upper limit value
                         **bhspinmag = 0.0**
=======================  ===============================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; REMNANT SPIN FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; bhspinflag uses different prescriptions for BH spin after formation 
    ; bhspinflag=0; sets all BH spins to bhspinmag
    ; bhspinflag=1; draws a random BH spin between 0 and bhspinmag for every BH
    ; bhspinflag=2; core-mass dependent BH spin (based on Belczynski+2017; 1706.07053, v1)
    ; default=0 
    bhspinflag = 0
    
    ; bhspinmag sets either the spin of all BHs or the upper limit of the uniform 
    ; distribution for BH spins 
    ; default=0.0
    bhspinmag=0.0

.. note::

    MASS TRANSFER FLAGS

=======================  =====================================================
``eddfac``               Eddington limit factor for mass transfer. 

                            ``1`` : mass transfer rate is limited by the 
                            Eddington rate following Equation 67 in 
                            `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``values >1`` : permit super-Eddington accretion 
                            up to value of *eddfac*

                         **eddfac = 1.0**

``gamma``                Angular momentum prescriptions for mass lost during RLO
                         at super-Eddington mass transfer rates

                            ``-1`` : assumes the lost material carries away the 
                            specific angular momentum of the primary

                            ``-2`` : assumes material is lost from the system as 
                            if it is a wind from the secondary 

                            ``>0`` : assumes that the lost material takes away a 
                            fraction *gamma* of the orbital angular momentum

                         **gamma = -1.0**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; MASS TRANSFER FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; eddfac is Eddington limit factor for mass transfer
    ; default=1.0
    eddfac=1.0

    ; gamma is the angular momentum factor for mass lost during RLO
    ; gamma=-2: assumes material is lost from the system as if it is a wind from the secondary (for super-Eddington mass transfer rates)
    ; gamma=-1: assumes the lost material carries with is the specific angular momentum of the primary
    ; gamma>0: assumes that the lost material take away a fraction (gamma) of the orbital angular momentum
    ; default=-1
    gamma=-1.0

.. note::

    TIDES FLAGS

=======================  =====================================================
``tflag``                Activates tidal circularisation following 
                         `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``0`` : no tidal circularization

                            ``1`` : activates tidal circularization

                         **tflag = 1**

``ST_tide``              Activates StarTrack setup for tides following
                         `Belczynski+2008 <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                            ``0`` : follows `BSE <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_

                            ``1`` : follows `StarTrack <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_

                         **ST_tide = 1**

``fprimc_array``         controls the scaling factor for convective tides
                         each item is set individually for its associated kstar
                         The releveant equation is `Equation 21 of Hurley+2002 <https://watermark.silverchair.com/329-4-897.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAnAwggJsBgkqhkiG9w0BBwagggJdMIICWQIBADCCAlIGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMYUoYtydpxVKmZePqAgEQgIICI1b5IZldHg9_rX6JacIe-IR042LnNi-4F9DMp-2lm3djjQ8xehKOv5I0VBjSNJfa6n-FErAH7ed1llADY7tMDTvqo1GHKBMDslNku5XDGfmae0sF-Zp5ndeGoZsyqISABLHEbdY4VFl8Uz_6jzAuBjGztnuxVmUh9bKIOaxuDpfB3Mn2xOfP9lcCVkjzQ0JWzr98nQNmVwDkI9bPv98Ab46BjBdGdcBKajCC-sqASjtmAQS2h6SGTTBqyRAyigqXcPtWf3Ye1SbxtL3zag6_Lf01rgCoUCK9eT_pavb5F8vVkUTMWbZQ79DWxn5pfZYi72C7_BtlPoUnS8Gs3wvw18BTIaHTKblwh225DcXuTEh_ngMmRvPEVctvG8tjlr9md-eFK0cEsq0734eGYtnwxeqvFxcWsW6mRbXrFHFsInQK16j6n36XuCimY665l_-HPAuu-lTTlwpMTUR7K1eYMBsco_tp_TdxEipRNvBpaWZX3J0FxPMzi84Y01UvWiW69pxb-LLTpf8aG4YCm9asRFyfDZ9nbSdgrIlCiuzy7QSmkvsHOaTEecmwRimFRycDuIuWLvA_tILmYCIM2KzvqYJSVCQPJH39xEHZG8LbMqImwAVYO3H90qh-90gNrtZn4ofSskcgqxeqfZly9CPfmEevX5s-SlLHMh1N6gdZwenvMC0kTWg_rskbvGiANtuGngD-kKDbunGpYJU_nI7uDnhGtdY#page=5>`_

                            ``positive values`` : sets scaling factor of
                            Equation 21 referenced above

                         **fprimc_array = [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                         2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,
                         2.0/21.0,2.0/21.0]**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;
    ;;; TIDES FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;

    ; tflag=1 activates tidal circularisation
    ; default=1
    tflag=1

    ; ST_tide sets which tidal method to use. 0=Hurley+2002, 1=StarTrack: Belczynski+2008
    ; Note, here startrack method does not use a better integration scheme (yet) but simply
    ; follows similar set up to startrack (including initial vrot, using roche-lobe check
    ; at periastron, and circularisation and synchronisation at start of MT).
    ; default=1
    ST_tide=1

    ; fprimc_array controls the scaling factor for convective tides
    ; each item is set individually for its associated kstar
    ; The releveant equation is Equation 21 from the BSE paper
    ; The default is to send the same coefficient (2/21) as is in the equation
    ; for every kstar
    fprimc_array=[2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]

.. note::

    WHITE DWARF FLAGS

=======================  =====================================================
``ifflag``               Activates the initial-final white dwarf mass relation 
                         from Han+1995 `Equations 3, 4, and 5 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1995MNRAS.272..800H&link_type=ARTICLE&db_key=AST&high=#page=4>`_.

                            ``0`` : no modifications to BSE

                            ``1`` : activates initial-final WD mass relation

                         **ifflag = 0**
``wdflag``               Activates an alternate cooling law found in the description 
                         immediately following `Equation 1 <http://iopscience.iop.org/article/10.1086/374637/pdf#page=3>`_
                         in Hurley & Shara 2003. 
                         Equation 1 gives the BSE default Mestel cooling law.

                            ``0`` : no modifications to BSE

                            ``1`` : activates modified cooling law

                         **wdflag = 1**
``epsnov``               Fraction of accreted matter retained in a nova eruption.
                         This is relevant for accretion onto degenerate objects; 
                         see Section 2.6.6.2 in `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_.

                            ``positive values between [0, 1]`` : retains *epsnov*
                            fraction of accreted matter

                         **epsnov = 0.001**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; WHITE DWARF FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;

    ; ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800
    ; default=0
    ifflag=0

    ; wdflag > 0 uses modified-Mestel cooling for WDs
    ; default=1
    wdflag=1

    ; epsnov is the fraction of accreted matter retained in nova eruptions
    ; default=0.001
    epsnov=0.001

.. note::

    PULSAR FLAGS

=======================  =====================================================
``bdecayfac``            Activates different models for accretion induced field decay; see 
                         `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``0`` : uses an exponential decay

                            ``1`` : uses an inverse decay

                         **bdecayfac = 1**
``bconst``               Sets the magnetic field decay time-scale for pulsars following
                         Section 3 of `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``negative values`` : sets k in Myr from Equation 8 to
                            -1 * *bconst*

                         **bconst = -3000**
``ck``                   Sets the magnetic field decay time-scale for pulsars following
                         Section 3 of `Kiel+2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

                            ``negative values`` : sets :math:`{\tau}`\ :sub:`b` in Myr
                            from Equation 2 to  -1 * *ck* 

                         **ck = -1000**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;
    ;; PULSAR FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;

    ; bdecayfac determines which accretion induced field decay method to 
    ; use from Kiel+2008: 0=exp, 1=inverse
    ; default=1
    bdecayfac=1

    ; bconst is related to magnetic field evolution of pulsars, see Kiel+2008
    ; default=-3000
    bconst=-3000

    ; ck is related to magnetic field evolution of pulsars, see Kiel+2008
    ; default=-1000
    ck=-1000

.. note::

    MIXING VARIABLES

=======================  =====================================================
``rejuv_fac``            Sets the mixing factor in main sequence star collisions. 
                         This is hard coded to 0.1 in the original BSE release
                         and in Equation 80 of `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_
                         but can lead to extended main sequence lifetimes in some cases.

                             ``positive values`` : sets the mixing factor

                         **rejuv_fac = 1.0**
``rejuvflag``            Sets whether to use the orginal prescription for mixing 
                         of main-sequence stars (based on equation 80 of `Hurley+2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_)
                         or whether to use the ratio of the pre-merger He core
                         mass at the base of the giant branch to the merger product's
                         He core mass at the base of the giant branch
                         
                         
                            ``0`` : no modifications to BSE

                            ``1`` : modified mixing times 

                         **rejuvflag = 0**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;
    ;; MIXING VARIABLES ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;

    ; rejuv_fac allows different mixing factors in Equation 80 from the BSE
    ; paper. This was originally hard coded to 0.1, which leads massive 
    ; stars to potentially have extended main sequence lifetimes. 
    rejuv_fac=1.0

    ; rejuvflag toggles between the original BSE prescription for MS mixing and 
    ; lifetimes of stars based on the mass of the MS stars (equation 80) or a
    ; prescription that uses the ratio of helium core mass of the pre-merger stars
    ; at the base of the first ascent of the giant branch to determine relative to the
    ; helium core mass of the merger product at the base of the giant branch
    ; default=0
    rejuvflag=0

.. note::

    MAGNETIC BRAKING FLAGS

=======================  =====================================================
``htpmb``                Activates different models for magnetic braking

                            ``0`` : no modifications to BSE

                            ``1`` : follows `Ivanona and Taam 2003 <https://ui.adsabs.harvard.edu/abs/2003ApJ...599..516I/abstract>`_

                         **htpmb = 1**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; MAGNETIC BRAKING FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ; htpmb allows for different magnetic braking models. 
    ; 0=follows BSE paper Section 2.4
    ; 1=follows Ivanova & Taam 2003 method which kicks in later than the standard
    ; default=1
    htpmb=1


.. note::

    MISC FLAGS

=======================  =====================================================
``ST_cr``                Activates different convective vs radiative boundaries

                            ``0`` : no modifications to BSE

                            ``1`` : follows `StarTrack <https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract>`_
 
                         **ST_cr = 1**
=======================  =====================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;
    ;; MISC FLAGS ;;;
    ;;;;;;;;;;;;;;;;;
    
    ; ST_cr sets which convective/radiative boundary to use
    ; 0=follows BSE paper
    ; 1=follows StarTrack (Belcyznski+2008)
    ; default=1
    ST_cr=1

