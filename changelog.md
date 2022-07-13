# COSMIC Changelog

## 0.0 - 3.4
See the discussed changes in our previous releases here: https://github.com/COSMIC-PopSynth/COSMIC/releases

## 3.4.1
 - minor fixes to documentation to stop warnings/errors which break doc website deployment
 - adds in build wheels for Python 3.9 and 3.10
 - add warning for mlwind prescriptions which assume zsun = 0.019
 - modifications to how m2min and qmin are specified in initial conditions sampling
 - fixed sampling issues from Plummer profiles in cmc.py
 - fixed Eddington limits such that wind accretion is limited to eddington and thus limits RLO accretion if the wind accretion is already super-Eddington
 - fixed errors in Marchant+ PISN prescription
