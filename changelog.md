# COSMIC Changelog

## 0.0 - 3.4
See the discussed changes in our previous releases here: https://github.com/COSMIC-PopSynth/COSMIC/releases

## 3.4.1
 - fixed sampling issues for renzo19 and sana12 so that they work with CMC
 - update to Docker ubuntu image
 - minor fixes to documentation to stop warnings/errors which break doc website deployment
 - adds in build wheels for Python 3.9 and 3.10
 - add warning for mlwind prescriptions which assume zsun = 0.019
 - modifications to how m2min and qmin are specified in initial conditions sampling
 - fixed sampling issues from Plummer profiles in cmc.py
 - fixed Eddington limits such that wind accretion is limited to eddington and thus limits RLO accretion if the wind accretion is already super-Eddington
 - fixed errors in Marchant+ PISN prescription

## 3.4.1-3.4.5
 - updated the actions to build properly

## 3.4.6
 - added new Eddington limit prescriptions for BH accretion
 - cleaned up versioning

## 3.4.7
 - Exchanged any print statements for `warnings.warn` to allow silencing
 - Fixed up docs to catch warnings

## 3.4.8
 - removed zsun_wind so that we can match stock BSE exactly and added extra documentation surrounding the winds
 - added a NaN catch in cosmic-pop that throws a warning of NaNs, saves them to a file, and instructs the user to consider changing pts1 since this is the main driver of NaNs so far

## 3.4.9
 - Added `sampling_target == "total_mass"` option to independent sampler so you can target a specific total mass instead of number of binaries
 - Added `trim_extra_samples` parameter to the same function - which trims your samples to get as close as possible to the target total mass
 - Bug fixes: secondaries of single stars are now marked as massless remnants instead of main sequence stars, binfrac=0.0 no longer leads to an infinite loop in sampling

## 3.4.10
 - Bug fixes: `timestep_conditions` for `bcm` arrays now raise errors for invalid columns instead of failing silently
 - Add `teff_1` and `teff_2` as variables that can be used to set `timestep_conditions`
 - Add in `-1` option to turn off Magnetic Braking in htmb 
 - Added `central_bh` and `scale_with_central_bh` as options to the CMC sampler, in order to add central massive black holes to CMC initial conditions


## 3.4.11
 - Added sampling options to ``independent`` sampler to allow for custom power law distributions for ``porb`` and ``q``

## 3.5.0
 - Feature: Added `bpp_columns` and `bcm_columns` parameters to the `evolve()` function to allow users to specify the columns in the bpp and bcm tables
 - Bug fix: Changed `kick.f` to use the Pfahl+02 kick prescription by default instead of Kiel & Hurley 2009, this fixes ejection velocities of secondaries and also changed kick_info to have an extra column

## 3.6.0
 - Overhaul documentation and added debugging environment
 - Feature: Added Disberg+2025 kick prescription as a new choice of `kickflag` (`kickflag=5`). Applies log-normal distribution to regular CCSN, ECSN/USSN still use `sigmadiv` Maxwellian and BH fallback scaling is still applied via `bhflag` and `bhsigmafrac` as with `kickflag=1`

## 3.6.1
 - Add support for single stars in both independent and multidim sampling
 - update documentation

## 3.6.2
 - Add functions to `cosmic.utils` for initC IO that's more efficient (`save_initC`, `load_initC`) by saving
 identical setting columns separately with only one copy - saves ~1kb per binary
 - Make He and CO core masses available in output
 - Change defaults; `qcflag=5` and `eccfac=10`
 - Prevent users from supplying both `qmin` and `m2_min`
 - Prevent users from supplying both `sep` and `porb` when sampling, calculate `porb` from `sep` as necessary
 - Define `binstate=-1` as single stars

## 3.7.0
 - Add `remnantflag=5`: The Mandel & Muller 2020 prescription for remnant masses
 - Add `kickflag=6`: The Mandel & Muller 2020 prescription for natal kicks. This can be tuned with `mm_mu_ns` and `mm_mu_bh`
 - Add `remnantflag=6`: The Maltsev+25 and Willcox+25 prescription for remnant masses. This can be tuned with `maltsev_mode` and `maltsev_fallback`
 - (Docs) Any tutorial that uses a plain BSEDict now uses one that is drawn from the cosmic-settings.json file to avoid missing changes and saves us updating each docs page


## 3.7.1
 - Bug fix [[#729](https://github.com/COSMIC-PopSynth/COSMIC/issues/729)]: ensure disruptions are logged with ``evol_type==11`` when SNe occur during CEs
 - Bug fix [[#725](https://github.com/COSMIC-PopSynth/COSMIC/issues/725)]: set ``tb=sep=0`` for cases where AIC caused a WD to explode and leave behind no remnant (``evolve_type==9``)
 - Bug fix [[#724](https://github.com/COSMIC-PopSynth/COSMIC/issues/724)]: remove bug where ``kstar=15`` was assigned its previous epoch mass after merging during a CE and the merger product goes SN

## 3.7.2
This release contains _several_ fixes to how CO core masses/remnant masses are handled. It also adds a new PISN prescription, windflag and LBV winds flag.

- Fixes:
    - Update ``mc_co`` and ``mc_he`` after adjusting ``mc = mcmax`` in ``hrdiag.f`` for stripped stars. This can be a fairly significant change, up to ~2 Msun. (Used to actually be up to 10 Msun because we added HeMS core mass growth)
    - Make ``evolv2.f`` actually log PISN (they get set to ``kw = 15`` so were missed previously)
    - Change ``loop = 40000`` (up from 20,000), otherwise 70 Msun stars run out of loops with new mass modifier
    - Changed all existing ``pisn`` prescriptions to use ``mc_tot`` instead of ``mcbagb`` for conditions
    - Update ``mc = mt`` for **ALL** remnants instead of just for certain flags
    - Only apply conversion from baryonic to gravitational mass once. Previous code **APPLIED IT TWICE!!** for core-collapse BHs! Create subroutine ``baryonic_to_gravitational_mass``, which is called throughout
    - Ensure stellar type changes during hrdiag are logged to the bpp before SN occurs

- Additions/changes:
    - Linearly grow the core mass of HeMS stars from 0 to the value at the start of HeHG and check for core-collapse for all He stars (not just HeHG and HeGB)
    - Add ``dt_mass_modifier`` argument to ``Evolve.evolve()``. This lets you specify a list of tuples ``(m_low, m_high, mod)``, which multiplies ``pts1/pts2/pts3`` by ``mod`` for stars with primary star (only care about primary, conservative) ZAMS masses between ``m_low <= m_ZAMS < m_high``
    - Added a ``windflag = -1`` option to turn off stellar winds entirely
    - Added ``pisn = -4`` for Renzo+22/Hendriks+23 prescription for PPISN mass loss
        - With options ``ppi_co_shift`` and ``ppi_extral_ml``
    - Add new setting ``LBV_flag`` which allows one to turn off LBV winds, use Hurley+2000, or use Belcyznski+2008
    - Change the default LBV winds to Hurley

- Code cleanup:
    - ``assign_remnant`` no longer takes ``mc_tot`` as a parameter, just get it from the common block
    - WD chandrasekhar check use ``mc_co`` instead of ``mc``
    - Remnant flags 0-4 inclusive use ``mc_co`` instead of ``mc``
    - Got rid of ``mcx`` in ``assign_remnant`` in favour of clearer ``m_proto`` and ``m_FeNi`` to match the papers
    - [Very minor] Fryer Rapid was using <= instead of < everywhere

- Documentation:
    - Start new settings gallery in the documentation
    - Tag settings/options with the version they were added in the docs page and auto link them to release