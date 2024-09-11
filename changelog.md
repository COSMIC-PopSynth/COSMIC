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