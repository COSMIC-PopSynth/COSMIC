.. _bse:

#############################################
Binary Stellar Evolution (BSE) Flag Deep Dive
#############################################

BSE Dict contains all the bse flags for different prescriptions.

****
neta
****

*neta* is the Reimers mass-loss coefficent.
`Equation 106 SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_ (due to a typo there's an extra :math:`{\eta}` out front. The rate is directly proportional to :math:`{\eta}`). See `Section Vb <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1978A%26A....70..227K&link_type=ARTICLE&db_key=AST&high=#page=12>`_ in Kudritzki R. P., Reimers D., 1978, A&A, 70, 227 for discussion. **Defaults to 0.5**.


*****
bwind
*****

*bwind* is the binary enhanced mass loss parameter. See `Equation 12 BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_. **Defaults to 0, inactive for single**.

******
hewind
******

*hewind* is the helium star mass loss parameter. 10\ :sup:`-13` hewind L\ :sup:`2/3` gives He star mass-loss. Equivalent to 1 - :math:`{\mu}` in the last equation on `page 19 of SSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2000MNRAS.315..543H&link_type=ARTICLE&db_key=AST&high=#page=19>`_. **Defaults to 1**.

******
alpha1
******

*alpha1* is the common-envelope efficiency parameter. It scales the efficiency of transferring orbital energy to the envelope. See `Equation 71 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_. **Defaults to 1**.

******
lambda
******

*lambda* is the binding energy factor for common envelope evolution. The initial binding energy of the envelope goes like 1 / :math:`{\lambda}`. See  `Equation 69 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=11>`_. **Defaults to 0.5**.

******
ceflag
******

Code commments say "*ceflag>0 activates the spin-energy correction in common-envelope evolution*". This isn't true. The code only cares whether or not ceflag == 3. See lines 83-96 in comenv.f. **Defaults to 0**.

*****
tflag
*****

*tflag*>0 activates tidal circularization. **Defaults to 1**.

******
ifflag
******

*ifflag*>0 activates the initial-final white dwarf mass relation from Han, Podsiadlowski & Eggleton, 1995, MNRAS, 272, 800 `Equations 3, 4, and 5 <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=1995MNRAS.272..800H&link_type=ARTICLE&db_key=AST&high=#page=4>`_ **Defaults to 0**.

******
wdflag
******

*wdflag*>0 activates the alternate cooling law found in the description immediately following `Equation 1 <http://iopscience.iop.org/article/10.1086/374637/pdf#page=3>`_ in Hurley & Shara, 2003, Apj, May 20. Equation 1 gives the default Mestel cooling law (wdflag=0). **Defaults to 0**.

******
nsflag
******

Documentation says "*nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1).*"
FIXME: Looking at the code, this isn't true. There are at least 3 prescriptions. Completing this document will require taking a closer look. More to come. **Defaults to 1**.

******
bhflag
******

*bhflag*>0 applies supernova kicks with Maxwellian dispersion velocity *sigma* (see below). Note: This is terrible nomenclature, as *bhflag* isn't at all analogous to *nsflag* (above). Might be less confusing to merge *bhflag* and *sigma* into a single parameter; the exception being a constant kick prescription, but how often is that used? **Defaults to 0.**

*****
sigma
*****

*sigma* is the Maxwellian dispersion velocity for the supernova kick. **Defaults to 190 km/s**.

****
mxns
****

*mxns* is the maximum neutron star mass. **Defaults to 1.8 if *nsflag*=0, and 3.0 if *nsflag*=1 (Note: nsflag defaults to 1)**.

****************
pts1, pts2, pts3 
****************

*pts1, pts2, pts3* scale the size of the time step depending on the ktype. FIXME: No documentation available. Further investigation required. More to come

****
beta
****

*beta* is the wind velocity factor. v\ :sub:`wind` :sup:`2` goes like *beta*. See `Equation 9 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_. **Defaults to 1/8**.

**
xi
**

*xi* is the wind accretion efficiency factor. It gives the fraction of angular momentum lost via winds from the primary that transfers to the spin angular momentum of the companion. Corresponds to :math:`{\mu}`\ :sub:`w` in `Equation 11 of BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=3>`_.

****
acc2
****

*acc2* is the Bondi-Hoyle wind accretion factor. The mean wind accretion rate onto the secondary is proportional to acc2. See `Equation 6 in BSE <http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2002MNRAS.329..897H&link_type=ARTICLE&db_key=AST&high=#page=2>`_.

******
epsnov
******

*epsnov*  is the fraction of accreted matter retained in a nova eruption, set by **default to 0.001**. This is relevant for accretion onto degenerate objects (See Section 2.6.6.2 in BSE paper)

******
eddfac
******

*eddfac* is Eddington limit factor for mass transfer, set by **default to 1.0**. In the case of eddfac=1, the mass transfer rate is limited by Eddington rate (Equation (67) in BSE paper). There is some uncertainty as to whether Eddington limit should be applied. Set eddfac >1 to permit some amount of super-Eddington accretion (Section 2.6.6.2 in BSE)

*****
gamma
*****

*gamma* is the angular momentum factor for mass lost during Roche lobe overflow. For super-Eddington mass transfer rates, for gamma = -2.0, and for novae systems, assume that material is lost from the system as if a wind from the secondary. If gamma = -1.0 then assume the lost material carries with it the specific angular momentum of the primary and for all gamma > 0.0 assume that it takes away a fraction gamma of the orbital angular momentum. **Default value is -1.0**.

*************
bconst and ck
*************

*bconst* and *CK* both pertain to neutron star (pulsar) evolution. Implemented by Paul Kiel -- see Section 3 of `Kiel et al. 2008 <https://academic.oup.com/mnras/article/388/1/393/1013977>`_.

******
merger
******

Not sure what this is, but ask Sourav.

********
windflag
********

*windflag* sets which wind prescription is to be used. 0=bse (as outlined in SSE paper), 1=StarTrack (`Belczynski et al. 2010 <http://iopscience.iop.org/article/10.1088/0004-637X/714/2/1217/meta>`_), 2=Vink (`Vink et al 2001 <http://adsabs.harvard.edu/abs/2001A&amp;A...369..574V>`_)
