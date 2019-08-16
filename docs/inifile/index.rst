.. _inifile:

##############################
Configuration files for COSMIC
##############################

How to write a configuration file
=================================

Here is a link to the most recent stable release version of the default
inifile for COSMIC: `Stable Version INIFILE <https://github.com/COSMIC-PopSynth/COSMIC/blob/master/examples/Params.ini>`_

Here is a link to the unstable development version of the default inifile for COSMIC: `Development Version INIFILE <https://github.com/COSMIC-PopSynth/COSMIC/blob/develop/examples/Params.ini>`_

The `runFixedPop` command-line executable cannot run without a configuration file.
Each of the below sections lists the valid options and a description of what the value should be, and then an example of what that section might look like in the INI format:

[filters]
---------

===================================  =====================================================
``mass_transfer_white_dwarf_to_co``  Do not retain systems which are some point during
                                     their evolution a white dwarf mass transferred
                                     onto a compact object
``select_final_state``               The bcm array generally returns the first and last
                                     state of the binary system. Since we already
                                     save the initial conditions in runFixedPop, usually
                                     we opt to throw out the initial state of the binary
                                     from this array and only keep the final state 
``binary_state``                     Each binary system at the end of its evolution
                                     will be in one of three states. Either it is still
                                     a binary system at the end of the evolution period,
                                     the sytem under went a merger, or the binary
                                     is no longer a binary as it was disrupted
                                     at some point in its evolution (likely due one
                                     of the objects undergoing a Supernova).
``lisa_sources``                     A boolean flag which, if True, will retain
                                     only sources within LISA's frequency
                                     range (porb<1e5 sec)
===================================  =====================================================

.. code-block:: ini

    [filters]
    ; mass_transfer_white_dwarf_to_co gets rid of any binary where this happens
    ; default=False
    mass_transfer_white_dwarf_to_co = False

    ; select_final_state will retain only the final entry of the bcm arry if set to True
    ; default=True
    select_final_state = True

    ; binary_state determines which types of binaries endstates to retain
    ; 0 alive today, 1 merged, 2 disrupted
    ; default=[0,1,2]
    binary_state = [0,1,2]

    ; lisa_sources decides whether to retain only sources within LISA's frequency range (porb<1e5 sec)
    ; default=False
    lisa_sources = False

[convergence]
-------------

====================  ============================================================
``lisa_convergence``  One of the steps performed during `runFixedPop`
                      is to check whether the properities of your evolved binaries
                      have converged to be representative of the underlying
                      distribution. This boolean flag, if True, restricts the
                      systems it checks further to be only systems with
                      porb < 5000 sec 
====================  ============================================================

.. code-block:: ini

    [convergence]
    ; perform convergence over porb < 5000 sec for LISA Galaxies
    ; default=False
    lisa_convergence=False

[rand_seed]
-----------

=============  ============================================================
``rand_seed``  Seed used to seed numpy.random.seed
=============  ============================================================

.. code-block:: ini

    [rand_seed]
    ; random seed int
    seed = 21

[bse]
-----

.. note::

    Although this is all one section, we have grouped the
    flags/parameters which get passed to the binary stellar evolution
    code into types. Each group will start with a note to indicate
    the type of parameter or flag.

.. note::

    SAMPLING FLAGS

========  ============================================================
``pts1``  determines the timesteps chosen in each evolution phase as
          decimal fractions of the time taken in that phase for
          Main Sequence (MS) stars (**default=0.001**, see Banerjee+ 2019)
``pts2``  determines the timesteps chosen in each evolution phase as
          decimal fractions of the time taken in that phase for 
          Giant Branch (GB, CHeB, AGB, HeGB) stars
          (**default=0.01**,)
``pts3``  determines the timesteps chosen in each evolution phase as
          decimal fractions of the time taken in that phase for 
          HG, HeMS stars (**default=0.02**,) 
========  ============================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;
    ;;; SAMPLING FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;

    ; dtp is the timestep (in Myr) for outputting to the bcm array
    ; if dtp=0, will print every timestep (not recommended)
    ; if not set, it will automatically set to dtp=tphsyf (default)
    ;dtp = 1.0

    ; pts1,pts2,pts3 determine the timesteps chosen in each
    ; evolution phase as decimal fractions of the time taken in that phase:
    ;                 pts1 - MS                  (default=0.001, see Banerjee+ 2019)
    pts1=0.001
    ;                 pts2 - GB, CHeB, AGB, HeGB (default=0.01)
    pts2=0.01
    ;                 pts3 - HG, HeMS            (default=0.02)
    pts3=0.02

.. note::

    WIND FLAGS

==============  ============================================================
``windflag``    windflag=0: bse (as outlined in SSE paper),

                windflag=1: StarTrack (`Belczynski et al. 2010 <http://iopscience.iop.org/article/10.1088/0004-637X/714/2/1217/meta>`_)

                windflag=2: Vink (`Vink et al 2001 <http://adsabs.harvard.edu/abs/2001A&amp;A...369..574V>`_)

                windflag=3: Vink+2005 (Vink plus LBV winds)

                **default=3**
``eddlimflag``
                eddlimflag turns on metallicity dependence on winds, affecting the
                mass-loss rate of low-metallicity stars near the Eddington limit
                (see Grafener et al. 2011, Giacobbo et al. 2017)

                eddlimflag = 0 off (**default**)

                eddlimflag = 1 on
``neta``        *neta* is the Reimers mass-loss coefficent.
                `Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_ (due to a typo there's an extra :math:`{\eta}` out front. The rate is directly proportional to :math:`{\eta}`).
                See `Section Vb <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1978A%26A....70..227K&link_type=ARTICLE&db_key=AST&high=#page=12>`_ in Kudritzki R. P., Reimers D., 1978, A&A, 70, 227 for discussion.

                **default=0.5**
``bwind``       *bwind* is the binary enhanced mass loss parameter. See `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.
                **Defaults to 0, inactive for single**
``hewind``      *hewind* is the helium star mass loss parameter. 10\ :sup:`-13` hewind L\ :sup:`2/3` gives He star mass-loss. Equivalent to 1 - :math:`{\mu}` in the last equation on `page 19 of SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_.

                **default=1.0**
``beta``        *beta* is the wind velocity factor. v\ :sub:`wind` :sup:`2` goes like *beta*. See `Equation 9 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                beta<0: follows StarTrack 2008

                beta=0.125: BSE default

                **default=-1.0**
``xi``          *xi* is the wind accretion efficiency factor. It gives the fraction of angular momentum lost via winds from the primary that transfers to the spin angular momentum of the companion. Corresponds to :math:`{\mu}`\ :sub:`w` in `Equation 11 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

                **default=0.5**
``acc2``        *acc2* is the Bondi-Hoyle wind accretion factor. The mean wind accretion rate onto the secondary is proportional to acc2. See `Equation 6 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_.

                **default=1.5**
==============  ============================================================

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
    ; default=1.0
    hewind = 1.0

    ; beta is wind velocity factor: proportional to vwind^2
    ; beta<0: follows StarTrack 2008; beta=0.125: stock BSE
    ; default=-1.0
    beta=-1.0

    ; xi is the wind accretion efficiency factor, which gives the fraction of angular momentum lost via winds from the primary that transfers to the spin angular momentum of the companion
    ; default=0.5
    xi=0.5

    ; acc2 sets the Bondi-Hoyle wind accretion factor onto companion
    ; default=1.5
    acc2=1.5

.. note::

    COMMON ENVELOPE FLAGS

================  ============================================================
``alpha1``        *alpha1* is the common-envelope efficiency parameter. It scales the efficiency of transferring orbital energy to the envelope. See `Equation 71 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                  **default=1.0**
``lambdaf``       *lambda1* is the binding energy factor for common envelope evolution. The initial binding energy of the envelope goes like 1 / :math:`{\lambda}`. See  `Equation 69 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

                  lambdaf=1.0 uses variable lambda prescription written by Onno Pols

                  lambdaf<0 uses fixes lambda to a value of -1.0*lambdaf

                  **default=1.0**

``ceflag``        ceflag=1 used the method from de Kool 1990 for setting the initial orbital energy

                  ceflag=0 does not use this method (uses the core mass to calculate initial orbital energy)

                  **default=0** 
``cekickflag``    cekickflag determined the prescription for calling kick.f in comenv.f
                  0: default BSE

                  1: uses pre-CE mass and sep values

                  2: uses post-CE mass and sep

                  **default=0**

``cemergeflag``   cemergeflag determines whether stars without a core-envelope boundary automatically lead to merger in CE

                  cemergeflag=1 turns this on (causes these systems to merge)

                  **default=0**
``cehestarflag``  cehestarflag uses fitting formulae from TLP, 2015, MNRAS, 451 for evolving RLO systems with a helium star donor and compact object accretor
                  this flag will override choice made by cekickflag if set

                  0: off

                  1: fits for final period only

                  2: fits for both final mass and final period

                  **default=0**
``qcflag``        qcflag is an integer flag that sets the model to determine which critical mass ratios to use for the onset of unstable mass transfer and/or a common envelope. NOTE: this is overridden by qcrit_array if any of the values are non-zero.

                  0: standard BSE

                  1: BSE but with Hjellming & Webbink, 1987, ApJ, 318, 794 GB/AGB stars

                  2: following binary_c from Claeys+2014 Table 2

                  3: following binary_c from Claeys+2014 Table 2 but with Hjellming & Webbink, 1987, ApJ, 318, 794 GB/AGB stars

                  **default=3**

``qcrit_array``   qcrit_array is a 16-length array for user-input values for the critical mass ratios that govern the onset of unstable mass transfer and a common envelope. Each item is set individually for its associated kstar, and a value of 0.0 will apply prescription of the qcflag for that kstar
                  **default: [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]**
================  ============================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; COMMON ENVELOPE FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; alpha1 is the common-envelope efficiency parameter
    ; default=1.0
    alpha1 = 1.0

    ; lambdaf is the binding energy factor for common envelope evolution
    ; lambdaf=1.0 uses variable lambda prescription written by Onno Pols
    ; lambdaf<0 uses fixes lambda to a value of -1.0*lambdaf
    ; default=1.0
    lambdaf = 1.0

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

====================  ==================================================
``sigma``             sigma sets is the dispersion in the Maxwellian for the SN kick velocity in km/s
                      **default=265.0**
``bhflag``            bhflag != 0 allows velocity kick at BH formation

                      bhflag=0: no BH kicks

                      bhflag=1: fallback-modulated kicks

                      bhflag=2: mass-weighted (proportional) kicks

                      bhflag=3: full NS kicks 

                      **default=1**
``ecsn``              ecsn>0 turns on ECSN and also sets the maximum ECSN mass range (at the time of the SN)

                      stock BSE and StarTrack: ecsn=2.25

                      Podsiadlowski+2004: ecsn=2.5)
``ecsn_mlow``         ecsn_mlow sets the low end of the ECSN mass range

                      BSE=1.6

                      Podsiadlowski+2004=1.4

                      StarTrack=1.85
``sigmadiv``          sigmadiv sets the modified ECSN kick
                      negative values sets the ECSN sigma value, positive values divide sigma above by sigmadiv

                      **default=-20.0**
``aic``               aic=1 turns on low kicks for accretion induced collapse works even if ecsn=0

                      **default=1**
``ussn``              ussn=1 uses reduced kicks (drawn from the sigmadiv distritbuion) for ultra-stripped supernovae
                      these happen whenever a He-star undergoes a CE with a compact companion

                      **default=0**
``pisn``              pisn>0 allows for (pulsational) pair instability supernovae
                      and sets the maximum mass of the remnant

                      pisn=-1 uses the formulae from Spera+Mapelli 2017 for the mass

                      pisn=-2 uses a polynomial fit to Table 1 Marchant 2018

                      pisn=-3 uses a polynomial fit to Table 5 in Woosley 2019

                      pisn=0 turns off (pulsational) pair instability supernovae

                      **default=45.0**
``bhsigmafrac``       bhsigmafrac sets the fractional modification used for scaling down the sigma for BHs
                      this works in addition to whatever is chosen for bhflag, and is applied to the sigma beforehand these prescriptions are implemented
                      **default=1.0**
``polar_kick_angle``  polar_kick_angle sets the opening angle of the kick relative to the pole of the exploding star
                      this can range from 0 (strictly polar kicks) to 90 (fully isotropic kicks)
                      **default=90.0**
``natal_kick_array``  natal_kick_array is a 6-length array for user-input values for the SN natal kick
                      formatted as: (vk1, vk2, phi1, phi2, theta1, theta2)
                      vk is valid on the range [0, inf], phi are the co-lateral polar angles valid from [-pi/2, pi/2], and theta are azimuthal angles [0, 2*pi]
                      any number outside of these ranges will be sampled in the standard way in kick.f
                      **default=[-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]**
====================  ==================================================

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

==========  ============================================================
``nsflag``  nsflag determines the remnant mass prescription used

            nsflag=0: default BSE

            nsflag=1: Belczynski et al. 2002, ApJ, 572, 407

            nsflag=2: Belczynski et al. 2008

            nsflag=3: rapid prescription (Fryer+ 2012)

            nsflag=4: delayed prescription (Fryer+ 2012)

            **default=3**
``mxns``    mxns sets the maximum NS mass
            **default=3.0**
==========  ============================================================

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
    ; default=3.0
    mxns=3.0

.. note::

    MASS TRANSFER FLAGS

==========  ============================================================
``eddfac``  eddfac is Eddington limit factor for mass transfer. There is some uncertainty as to whether Eddington limit should be applied.
            In the case of eddfac=1, the mass transfer rate is limited by Eddington rate (Equation (67) in BSE paper).

            Set eddfac >1 to permit some amount of super-Eddington accretion (Section 2.6.6.2 in BSE)

            **default=1.0**

``gamma``   gamma is the angular momentum factor for mass lost during RLO

            gamma=-2: assumes material is lost from the system as if it is a wind from the secondary (for super-Eddington mass transfer rates)
            gamma=-1: assumes the lost material carries with is the specific angular momentum of the primary

            gamma>0: assumes that the lost material take away a fraction (gamma) of the orbital angular momentum

            **default=-2**
==========  ============================================================

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
    ; default=-2
    gamma=-2.0


.. note::

    MISCELLANEOUS FLAGS

================  ============================================================
``tflag``         *tflag* activates tidal circularisation.
                  **default=1**
``ifflag``        *ifflag* activates the initial-final white dwarf mass relation from Han, Podsiadlowski & Eggleton, 1995, MNRAS, 272, 800 `Equations 3, 4, and 5 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1995MNRAS.272..800H&link_type=ARTICLE&db_key=AST&high=#page=4>`_.
                  **default=0**
``wdflag``        *wdflag* activates the alternate cooling law found in the description immediately following `Equation 1 <http://iopscience.iop.org/article/10.1086/374637/pdf#page=3>`_ in Hurley & Shara, 2003, Apj, May 20. Equation 1 gives the default Mestel cooling law (wdflag=0).
                  **default=0**
``epsnov``        *epsnov* is the fraction of accreted matter retained in a nova eruption, set by **default to 0.001**. This is relevant for accretion onto degenerate objects (See Section 2.6.6.2 in BSE paper)
                  **default=0.001**
``bconst``        *bconst* related to magnetic field evolution of pulsars. Implemented by Paul Kiel -- see Section 3 of `Kiel et al. 2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.
                  **default=-3000**
``ck``            *ck* related to magnetic field evolution of pulsars, . Implemented by Paul Kiel -- see Section 3 of `Kiel et al. 2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.
                  **default=-1000**

``fprimc_array``  *fprimc_array* controls the scaling factor for convective tides
                  each item is set individually for its associated kstar
                  The releveant equation is `Equation 21 <https://watermark.silverchair.com/329-4-897.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAnAwggJsBgkqhkiG9w0BBwagggJdMIICWQIBADCCAlIGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMYUoYtydpxVKmZePqAgEQgIICI1b5IZldHg9_rX6JacIe-IR042LnNi-4F9DMp-2lm3djjQ8xehKOv5I0VBjSNJfa6n-FErAH7ed1llADY7tMDTvqo1GHKBMDslNku5XDGfmae0sF-Zp5ndeGoZsyqISABLHEbdY4VFl8Uz_6jzAuBjGztnuxVmUh9bKIOaxuDpfB3Mn2xOfP9lcCVkjzQ0JWzr98nQNmVwDkI9bPv98Ab46BjBdGdcBKajCC-sqASjtmAQS2h6SGTTBqyRAyigqXcPtWf3Ye1SbxtL3zag6_Lf01rgCoUCK9eT_pavb5F8vVkUTMWbZQ79DWxn5pfZYi72C7_BtlPoUnS8Gs3wvw18BTIaHTKblwh225DcXuTEh_ngMmRvPEVctvG8tjlr9md-eFK0cEsq0734eGYtnwxeqvFxcWsW6mRbXrFHFsInQK16j6n36XuCimY665l_-HPAuu-lTTlwpMTUR7K1eYMBsco_tp_TdxEipRNvBpaWZX3J0FxPMzi84Y01UvWiW69pxb-LLTpf8aG4YCm9asRFyfDZ9nbSdgrIlCiuzy7QSmkvsHOaTEecmwRimFRycDuIuWLvA_tILmYCIM2KzvqYJSVCQPJH39xEHZG8LbMqImwAVYO3H90qh-90gNrtZn4ofSskcgqxeqfZly9CPfmEevX5s-SlLHMh1N6gdZwenvMC0kTWg_rskbvGiANtuGngD-kKDbunGpYJU_nI7uDnhGtdY#page=5>`_
                  The default is keep the 2/21 coefficient value as seen in the equation.

                  **default=[2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]**
================  ============================================================

.. code-block:: ini

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; MISCELLANEOUS FLAGS ;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; tflag=1 activates tidal circularisation
    ; default=1
    tflag=1

    ; ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800
    ; default=0
    ifflag=0

    ; wdflag > 0 uses modified-Mestel cooling for WDs
    ; default=0
    wdflag=0

    ; epsnov is the fraction of accreted matter retained in nova eruptions
    ; default=0.001
    epsnov=0.001

    ; bconst is related to magnetic field evolution of pulsars, see Kiel+2008
    ; default=-3000
    bconst=-3000

    ; ck is related to magnetic field evolution of pulsars, see Kiel+2008
    ; default=-1000
    ck=-1000

    ; fprimc_array controls the scaling factor for convective tides
    ; each item is set individually for its associated kstar
    ; The releveant equation is Equation 21 from the BSE paper
    ; The default is to send the same coefficient (2/21) as is in the equation
    ; for every kstar
    fprimc_array=[2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]
