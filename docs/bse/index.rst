.. _bse:

#############################################
Binary Stellar Evolution (BSE) Flag Deep Dive
#############################################

BSE Dict contains all the bse flags for different prescriptions.

****
neta
****

*neta* is the Reimers mass-loss coefficent.
`Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_ (due to a typo there's an extra :math:`{\eta}` out front. The rate is directly proportional to :math:`{\eta}`). See `Section Vb <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1978A%26A....70..227K&link_type=ARTICLE&db_key=AST&high=#page=12>`_ in Kudritzki R. P., Reimers D., 1978, A&A, 70, 227 for discussion.

**Defaults to 0.5**.


*****
bwind
*****

*bwind* is the binary enhanced mass loss parameter. See `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

**Defaults to 0, inactive for single**.

******
hewind
******

*hewind* is the helium star mass loss parameter. 10\ :sup:`-13` hewind L\ :sup:`2/3` gives He star mass-loss. Equivalent to 1 - :math:`{\mu}` in the last equation on `page 19 of SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_.

**Defaults to 1**.

******
alpha1
******

*alpha1* is the common-envelope efficiency parameter. It scales the efficiency of transferring orbital energy to the envelope. See `Equation 71 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

**Defaults to 1**.

******
lambda
******

*lambda* is the binding energy factor for common envelope evolution. The initial binding energy of the envelope goes like 1 / :math:`{\lambda}`. See  `Equation 69 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_.

**Defaults to 1.0**.

**********
cekickflag
**********
*cekickflag*: prescription for calling what masses and sep to feed into kick.f
from comenv.f

0: default BSE

1: uses pre-CE mass and sep values

2: uses post-CE mass and sep

**Defaults to 0**.

***********
cemergeflag
***********

*cemergeflag*: whether stars without core-envelope boundary automatiically lead to merger in CE

0: off

1: on

**Defaults to 0**.

************
cehestarflag
************

*cehestarflag*: uses fitting formulae from TLP, 2015, MNRAS, 451 for evolving systems involving a mass-transferring He-star witha compact object. This flag will override choice made by cekickflag if set

0: off

1: fits for final period only

2: fits for both final mass and final period

**Defaults to 0**.

*****
tflag
*****

*tflag*>0 activates tidal circularization.

**Defaults to 1**.

******
ifflag
******

*ifflag*>0 activates the initial-final white dwarf mass relation from Han, Podsiadlowski & Eggleton, 1995, MNRAS, 272, 800 `Equations 3, 4, and 5 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1995MNRAS.272..800H&link_type=ARTICLE&db_key=AST&high=#page=4>`_

**Defaults to 0**.

******
wdflag
******

*wdflag*>0 activates the alternate cooling law found in the description immediately following `Equation 1 <http://iopscience.iop.org/article/10.1086/374637/pdf#page=3>`_ in Hurley & Shara, 2003, Apj, May 20. Equation 1 gives the default Mestel cooling law (wdflag=0).

**Defaults to 0**.

******
bhflag
******

*bhflag* > 0 allows velocity kick at BH formation (0).

bhflag=0: no BH kicks;

bhflag=1: fallback-modulated kicks

bhflag=2: mass-weighted (proportional) kicks

**Defaults to 3**.

******
nsflag
******

nsflag=0: default BSE;

nsflag=1: Belczynski et al. 2002, ApJ, 572, 407

nsflag=2: Belczynski et al. 2008;

nsflag=3: rapid prescription (Fryer+ 2012)

nsflag=4: delayed prescription (Fryer+ 2012)

**Defaults to 3**.

****
mxns
****

*mxns* is the maximum neutron star mass.

**Defaults to 1.8 if *nsflag*=0, and 3.0 if *nsflag*=1 (Note: nsflag defaults to 1)**.

****************
pts1, pts2, pts3 
****************

These parameters determine the timesteps chosen in each
evolution phase as decimal fractions of the time taken in that phase:
                 pts1 - MS                  (0.001, see Banerjee+ 2019)
**pts1 defaults to 0.001**.
                 pts2 - GB, CHeB, AGB, HeGB (0.01)
**pts2 defaults to 0.01**.
                 pts3 - HG, HeMS            (0.02)
**pts3 defaults to 0.02**.

****
ppsn
****

Pair Instability supernova flag

0: off

1: on, set to M=45 Msun

**Defaults to 1**.

*****
ecsnp
*****

*ecsnp*>0 turns on ECSN and also sets the maximum ECSN mass range 
(mass at the time of the SN)

BSE/StarTrack=2.25

Podsiadlowski+2004=2.5

**Defaults to 2.5**.

*********
ecsn_mlow
*********

ecsn_mlow sets the low end of the ECSN mass range 

BSE=1.6

Podsiadlowski+2004=1.4

StarTrack=1.85

**Defaults to 1.4**.

***
aic
***

aic=0: off

aic=1 includes AIC low kicks even if ecsnp=0

**Defaults to 1**.

*****
sigma
*****

sigma is the dispersion in the Maxwellian for the SN kick speed (265 km/s)

**Defaults to 265.0 km/s**.

********
sigmadiv
********

sigmadiv sets the ECSN kick, negative values sets the ECSN sigma value 
to sigmadiv and positive divides sigma above by sigmadiv

**Defaults to -20.0 km/s**.

***********
bhsigmafrac
***********

bhsigmafrac sets the fractional modification to sigma for BHs (0<=bhsigmafrac<=1.0)

**Defaults to 1.0**.

****************
polar_kick_angle
****************

*polar_kick_angle* sets the opening angle of the kick relative to 
the pole of the exploding star, value of 90 degrees is 
the default for isotropic kicks (0<=polar_kick_angle<=90.0)

**Defaults to 90.0 degrees**.

****
beta
****

*beta* is the wind velocity factor. v\ :sub:`wind` :sup:`2` goes like *beta*. See `Equation 9 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

beta<0: follows StarTrack 2008

beta=0.125: BSE default

**Defaults to 0.125**.

**
xi
**

*xi* is the wind accretion efficiency factor. It gives the fraction of angular momentum lost via winds from the primary that transfers to the spin angular momentum of the companion. Corresponds to :math:`{\mu}`\ :sub:`w` in `Equation 11 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

**Defaults to 0.5**.

****
acc2
****

*acc2* is the Bondi-Hoyle wind accretion factor. The mean wind accretion rate onto the secondary is proportional to acc2. See `Equation 6 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_.

**Defaults to 1.5**.

******
epsnov
******

*epsnov*  is the fraction of accreted matter retained in a nova eruption, set by **default to 0.001**. This is relevant for accretion onto degenerate objects (See Section 2.6.6.2 in BSE paper)

******
eddfac
******

*eddfac* is Eddington limit factor for mass transfer. There is some uncertainty as to whether Eddington limit should be applied.

In the case of eddfac=1, the mass transfer rate is limited by Eddington rate (Equation (67) in BSE paper).

Set eddfac >1 to permit some amount of super-Eddington accretion (Section 2.6.6.2 in BSE)

**default to 1.0**

*****
gamma
*****

*gamma* is the angular momentum factor for mass lost during RLO 

gamma=-2: assumes material is lost from the system as if it is a wind 
from the secondary (for super-Eddington mass transfer rates)

gamma=-1: assumes the lost material carries with is the specific angular
momentum of the primary

gamma>0: assumes that the lost material take away a fraction (gamma) of`
the orbital angular momentum

**Defaults to -2.0**.


*************
bconst and ck
*************

*bconst* and *CK* both pertain to neutron star (pulsar) evolution. Implemented by Paul Kiel -- see Section 3 of `Kiel et al. 2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

**Defaults to -3000 and -1000**.

********
windflag
********

*windflag* sets which wind prescription is to be used.

0=bse (as outlined in SSE paper),

1=StarTrack (`Belczynski et al. 2010 <http://iopscience.iop.org/article/10.1088/0004-637X/714/2/1217/meta>`_)

2=Vink (`Vink et al 2001 <http://adsabs.harvard.edu/abs/2001A&amp;A...369..574V>`_)

windflag=3: Vink+2005 (Vink plus LBV winds)

**Defaults to 3**.

****************
natal_kick_array
****************

6-length array for user-input values for the SN natal kick of format:
(vk1, vk2, phi1, phi2, theta1, theta2)
vk is valid on the range [0, inf], phi are the co-lateral polar angles valid from [-pi/2, pi/2],
and theta are azimuthal angles [0, 2*pi]; any number outside of these ranges
will be sampled in the standard way in kick.f

**Defaults to [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0]**.

***********
qcrit_array
***********

16-length array for user-input values for the critical mass ratios that govern 
the onset of unstable mass transfer and a common envelope
Each item is set individually for its associated kstar, and a value of 0.0
will apply the standard BSE prescription for that kstar

**Defaults to [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]**
