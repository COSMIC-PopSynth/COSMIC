***
      SUBROUTINE evolv2(kstar,mass,tb,ecc,z,tphysf,
     \ dtp,mass0,rad,lumin,massc,radc,
     \ menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
     \ bhspin,tphys,zpars,bkick,kick_info,
     \ bpp_index_out,bcm_index_out,kick_info_out)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      INCLUDE 'checkstate.h'
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Adapted from Aarseth's code 21st September 1996.
* Fully revised on 27th November 1996 to remove vestiges of N-body code and
* incorporate corrections.
* Fully revised on 1st April 1998 to include new stellar evolution formulae
* and associated binary evolution changes.
* Fully revised on 4th July 1998 to include eccentricity, tidal
* circularization, wind accretion, velocity kicks for supernovae and all
* associated orbital momentum changes.
*
***
*
* See Tout et al., 1997, MNRAS, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
*
* Reference for the stellar evolution formulae is Hurley, Pols & Tout,
* 2000, MNRAS, 315, 543 (SSE paper).
* Reference for the binary evolution algorithm is Hurley, Tout & Pols,
* 2002, MNRAS, 329, 897 (BSE paper).
*
***
*
* March 2001 *
* Changes since version 3, i.e. since production of Paper3:
*
* 1) The Eddington limit flag (on/off) has been replaced by an
*    Eddington limit multiplicative factor (eddfac). So if you
*    want to neglect the Eddington limit you would set eddfac
*    to a large value.
*
* 2) To determine whether material transferred during RLOF forms
*    an accretion disk around the secondary or hits the secondary
*    in a direct stream we calculate a minimum radial distance, rmin,
*    of the mass stream from the secondary. This is taken from eq.(1)
*    of Ulrich & Burger (1976, ApJ, 206, 509) which they fitted to
*    the calculations of Lubow & Shu (1974, ApJ, 198, 383).
*    If rmin is less than the radius of the secondary then an
*    accretion disk is not formed.
*    Note that the formula for rmin given by Ulrich & Burger is valid
*    for all q whereas that given by Nelemans et al. (2001, A&A,
*    submitted) in their eq.(6) is only valid for q < 1 where
*    they define q = Mdonor/Maccretor, i.e. DD systems.
*
* 3) The changes to orbital and spin angular momentum owing to
*    RLOF mass transfer have been improved, and an new input option
*    now exists.
*    When mass is lost from the system during RLOF there are now
*    three choices as to how the orbital angular momentum is
*    affected: a) the lost material carries with it a fraction
*    gamma of the orbital angular momentum, i.e.
*    dJorb = gamma*dm*a^2*omega_orb; b) the material carries with it
*    the specific angular momentum of the primary, i.e.
*    dJorb = dm*a_1^2*omega_orb; or c) assume the material is lost
*    from the system as if a wind from the secondary, i.e.
*    dJorb = dm*a_2^2*omega_orb.
*    The parameter gamma is an input option.
*    Choice c) is used if the mass transfer is super-Eddington
*    or the system is experiencing novae eruptions.
*    In all other cases choice a) is used if gamma > 0.0, b) if
*    gamma = -1.0 and c) is used if gamma = -2.0.
*    The primary spin angular momentum is reduced by an amount
*    dm1*r_1^2*omega_1 when an amount of mass dm1 is transferred
*    from the primary.
*    If the secondary accretes through a disk then its spin
*    angular momentum is altered by assuming that the material
*    falls onto the star from the inner edge of a Keplerian
*    disk and that the system is in a steady state, i.e.
*    an amount dm2*SQRT(G*m_2*r_2).
*    If there is no accretion disk then we calculate the angular
*    momentum of the transferred material by using the radius at
*    at which the disk would have formed (rdisk = 1.7*rmin, see
*    Ulrich & Burger 1976) if allowed, i.e. the angular momentum
*    of the inner Lagrangian point, and add this directly to
*    the secondary, i.e. an amount dm2*SQRT(G*m_2*rdisk).
*    Total angular momentum is conserved in this model.
*
* 4) Now using q_crit = 3.0 for MS-MS Roche systems (previously we
*    had nothing). This corresponds roughly to R proportional to M^5
*    which should be true for the majority of the MS (varies from
*    (M^17 -> M^2). If q > q_crit then contact occurs.
*    For CHeB primaries we also take q_crit = 3.0 and allow
*    common-envelope to occur if this is exceeded.
*
* 5) The value of lambda used in calculations of the envelope binding
*    energy for giants in common-envelope is now variable (see function
*    in zfuncs). The lambda function has been fitted by Onno to detailed
*    models ... he will write about this soon!
*
* 6) Note that eq.42 in the paper is missing a SQRT around the
*    MR^2/a^5 part. This needs to be corrected in any code update
*    paper with a thanks to Jeremy Sepinsky (student at NorthWestern).
*    It is ok in the code.
*
* March 2003 *
* New input options added:
*
*    ifflag - for the mass of a WD you can choose to use the mass that
*             results from the evolution algorithm (basically a competition
*             between core-mass growth and envelope mass-loss) or use the IFMR
*             proposed by Han, Podsiadlowski & Eggleton, 1995, MNRAS, 272, 800
*             [>0 activates HPE IFMR].
*
*    wdflag - for the cooling of WDs you can choose to use either the standard
*             Mestel cooling law (see SSE paper) or a modified-Mestel law that
*             is better matched to detailed models (provided by Brad Hansen
*             ... see Hurley & Shara, 2003, ApJ, May 20, in press)
*             [>0 activates modified-Mestel].
*
*    bhflag - choose whether or not black holes should get velocity kicks
*             at formation
*             [0= no kick; >0 kick].
*
*    remnantflag - for the mass of neutron stars and black holes you can use either
*             the SSE prescription or the prescription presented by
*             Belczynski et al. 2002, ApJ, 572, 407 who found that SSE was
*             underestimating the masses of these stars. In either case you also
*             need to set the maximum NS mass (mxns) for the prescription
*             [0= SSE, mxns=1.8; >0 Belczynski, mxns=3.0].
*
* Sept 2004 *
* Input options added/changed:
*
*    ceflag - set to 1 this uses de Kool (or Podsiadlowski) CE prescription,
*             other options, such as Yungelson, could be added as well.
*
*    hewind - factor to control the amount of He star mass-loss, i.e.
*             1.0e-13*hewind*L^(2/3) gives He star mass-loss.
*
*
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
*
      INTEGER loop,iter,intpol,k,ip,jp,j1,j2
      INTEGER bcm_index_out, bpp_index_out
      INTEGER kcomp1,kcomp2,formation(2)
      PARAMETER(loop=20000)
      INTEGER kstar(2),kw,kst,kw1,kw2,kmin,kmax
      INTEGER kstar1_bpp,kstar2_bpp
*
      REAL*8 km,km0,tphys,tphys0,dtm0,tphys00,tphysfhold
      REAL*8 tphysf,dtp,tsave,dtp_original
      REAL*8 aj(2),aj0(2),epoch(2),tms(2),tbgb(2),tkh(2),dtmi(2)
      REAL*8 mass0(2),mass(2),massc(2),menv(2),mass00(2),mcxx(2)
      REAL*8 mass1_bpp,mass2_bpp
      REAL*8 rad(2),rol(2),rol0(2),rdot(2),radc(2),renv(2),radx(2)
      REAL*8 lumin(2),k2str(2),q(2),dms(2),dmr(2),dmt(2)
      REAL*8 dml,vorb2,vwind2,omv2,ivsqm,lacc,kick_info(2,18)
      REAL*8 bkick(20)
      REAL*8 kick_info_out(2,18)
      REAL*8 sep,dr,tb,dme,tdyn,taum,dm1,dm2,dmchk,qc,dt,pd,rlperi
      REAL*8 m1ce,m2ce,mch,tmsnew,dm22,mew
      PARAMETER(mch=1.44d0)
      REAL*8 yeardy,yearsc,aursun
      PARAMETER(yeardy=365.24d0,aursun=214.95d0,yearsc=3.1557d+07)
      REAL*8 acc1,tiny
      PARAMETER(acc1=3.920659d+08,tiny=1.0d-14)
      REAL*8 ecc,ecc1,tc,tcirc,ttid,ecc2,omecc2,sqome2,sqome3,sqome5
      REAL*8 f1,f2,f3,f4,f5,f,raa2,raa6,eqspin,rg2,tcqr,gammadisc
      REAL*8 k3,mr23yr,twopi
      PARAMETER(k3=0.21d0,mr23yr=0.4311d0)
      REAL*8 jspin(2),ospin(2),jorb,oorb,jspbru,ospbru
      REAL*8 bhspin(2)
      REAL*8 delet,delet1,dspint(2),djspint(2),djtx(2)
      REAL*8 dtj,djorb,djgr,djmb,djt,djtt,rmin,rdisk
      REAL*8 etaBH,maxspinBH
*
      INTEGER pulsar
      INTEGER mergemsp,merge_mem,notamerger,binstate,mergertype
      REAL*8 fallback,sigmahold
      REAL*8 vk,u1,u2,s,Kconst,betahold,convradcomp(2),teff(2)
      REAL*8 B_0(2),bacc(2),tacc(2),xip,xihold
      REAL*8 deltam1_bcm,deltam2_bcm,b01_bcm,b02_bcm
      REAL*8 B(2),Bbot,omdot,b_mdot,b_mdot_lim,evolve_type
      COMMON /fall/fallback
      REAL ran3
      EXTERNAL ran3
*

*
      REAL*8 z,tm,tn,m0,mt,rm,lum,mc,rc,me,re,k2,age,dtm,dtr
      REAL*8 tscls(20),lums(10),GB(10),zpars(20)
      REAL*8 zero,ngtv,ngtv2,mt2,rrl1,rrl2,mcx,teff1,teff2
      REAL*8 mass1i,mass2i,tbi,ecci
      LOGICAL coel,com,prec,inttry,change,snova,sgl
      LOGICAL supedd,novae,disk,inspiral
      LOGICAL iplot,isave
      REAL*8 rl,mlwind,vrotf,corerd,f_fac
      EXTERNAL rl,mlwind,vrotf,corerd
*
      REAL*8 kw3,wsun,wx
      PARAMETER(kw3=619.2d0,wsun=9.46d+07,wx=9.46d+08)
      LOGICAL output
*
      REAL*8 qc_fixed
      LOGICAL switchedCE,disrupt

Cf2py intent(in) kstar
Cf2py intent(in) mass
Cf2py intent(in) tb
Cf2py intent(in) ecc
Cf2py intent(in) z
Cf2py intent(in) tphysf
Cf2py intent(in) dtp
Cf2py intent(in) mass0
Cf2py intent(in) rad
Cf2py intent(in) lumin
Cf2py intent(in) massc
Cf2py intent(in) radc
Cf2py intent(in) menv
Cf2py intent(in) renv
Cf2py intent(in) ospin
Cf2py intent(in) B_0
Cf2py intent(in) bacc
Cf2py intent(in) tacc
Cf2py intent(in) epoch
Cf2py intent(in) tms
Cf2py intent(in) bhspin
Cf2py intent(in) tphys
Cf2py intent(in) zpars
Cf2py intent(in) bkick
Cf2py intent(in) kick_info
Cf2py intent(out) bpp_index_out
Cf2py intent(out) bcm_index_out
Cf2py intent(out) kick_info_out

      if(using_cmc.eq.0)then
              CALL instar
      endif

*
* Save the initial state.
*

*      CE2flag = 0
      kstar1_bpp = 0
      kstar2_bpp = 0

      mass1_bpp = 0.d0
      mass2_bpp = 0.d0

      mass1i = mass(1)
      mass2i = mass(2)
      tbi = tb
      ecci = ecc
*
      zero = 0.d0
      ngtv = -1.d0
      ngtv2 = -2.d0
      twopi = 2.d0*ACOS(-1.d0)

      Mbh_initial = 0.d0


* disrupt tracks if system get disrupted by a SN during the common
* envelope
      disrupt = .false.
* value for bcm[ii,37] which tracks binary state; 0 for binary, 1 for merger, 2 for disrupted
      binstate = 0
* value for bcm[ii,38] which tracks merger types; only set when binstate is 1
* the logic is to combine kstar values of merged objects. so 1313 or 0809.
      mergertype = -1
*Captures original sigma so after ECSN we can reset it.
      sigmahold = sigma
*memory for wind mass loss factor
      betahold = beta

** SET PULSAR VALUES HERE**
* PDK
      pulsar = 1 ! allows for pulsar physics; not a flag since we heart pulsars
      Kconst = 2.5d-49
      Bbot = 5e+7 !100.d0 or ~d+07.
      b_mdot_lim = -1.0e-11 !limiting accretion induced field decay with mdot as a proxy for
*                           accretion temperature and number of impurities.
      xihold = xi
      xip = 1 !Modifies NS wind ang. mom. accretion.
      formation(1) = 0 !helps determine formation channel of interesting systems.
      formation(2) = 0
      pisn_track(1) = 0 !tracks whether a PISN occurred for each
component.
      pisn_track(2) = 0
      merger = -1 !used in CMC to track systems that dynamically merge
      notamerger = 0 !if 0 you reset the merger NS product to new factory settings else you don't.
      mergemsp = 1 !if set to 1 any NS that merges with another star where the NS is an MSP stays an MSP...
      merge_mem = 0
      output = .false.  ! .true. turns on, .false. turns off.
*                       WARNING: can fill up the output file very quickly.
*                       With N=2e6 .stdout was 3.2 GB in 6 mins. If needed you can
*                       be more selective with outputting, but must add this yourself!
*      if(id1_pass.eq.1)then!.and.tphysf.gt.17.50d0.and.
**     &   tphysf.lt.17.6d0)then
*         output = .true.
*      endif
*
* Initialize the parameters.
*

      if(using_cmc.eq.0)then
          bcm_index_out = 0
          bpp_index_out = 0
          kick_info_out = 0.d0
      endif


*
* Set the seed for the random number generator.
*
*      idum1 = INT(sep*100)
      if(idum1.gt.0.and.using_cmc.eq.0) idum1 = -idum1

*
* Set the collision matrix.
*
      if(using_cmc.eq.0)then
          CALL zcnsts(z,zpars)
      endif

      kmin = 1
      kmax = 2
      sgl = .false.
      mt2 = MIN(mass(1),mass(2))
      kst = 0
      iter = 0 ! PDK addition, just incase you bail out before 4 loop (usually from multiple calls for one evolve)
*
      if(mt2.lt.tiny.or.tb.le.0.d0)then
         sgl = .true.
         if(mt2.lt.tiny)then
            mt2 = 0.d0
            if(mass(1).lt.tiny)then
               if(tphys.lt.tiny)then
                  mass0(1) = 0.01d0
                  mass(1) = mass0(1)
                  kst = 1
               else
                  kmin = 2
                  lumin(1) = 1.0d-10
                  rad(1) = 1.0d-10
                  massc(1) = 0.d0
                  dmt(1) = 0.d0
                  dmr(1) = 0.d0
               endif
               ospin(1) = 1.0d-10
               jspin(1) = 1.0d-10
               dtmi(1) = 5.0d+04 !5Gyr, added for use of single stars in evolv2.f
            else
               if(tphys.lt.tiny)then
                  mass0(2) = 0.01d0
                  mass(2) = mass0(2)
                  kst = 2
               else
                  kmax = 1
                  lumin(2) = 1.0d-10
                  rad(2) = 1.0d-10
                  massc(2) = 0.d0
                  dmt(2) = 0.d0
                  dmr(2) = 0.d0
               endif
               ospin(2) = 1.0d-10
               jspin(2) = 1.0d-10
               dtmi(2) = 5.0d+04 !5Gyr, added for use of single stars in evolv2.f
            endif
         endif
         ecc = -1.d0
         tb = 0.d0
         sep = 1.0d+10
         oorb = 0.d0
         jorb = 0.d0
         if(kstar(1).ne.14.d0.or.using_cmc.eq.0) bhspin(1) = 0.d0
         if(kstar(2).ne.14.d0.or.using_cmc.eq.0) bhspin(2) = 0.d0
         if(ospin(1).lt.0.0) ospin(1) = 1.0d-10
         if(ospin(2).lt.0.0) ospin(2) = 1.0d-10
         q(1) = 1.0d+10
         q(2) = 1.0d+10
         rol(1) = 1.0d+10
         rol(2) = 1.0d+10
      else
         tb = tb/yeardy
         sep = aursun*(tb*tb*(mass(1) + mass(2)))**(1.d0/3.d0)
         oorb = twopi/tb
         jorb = mass(1)*mass(2)/(mass(1)+mass(2))
     &          *SQRT(1.d0-ecc*ecc)*sep*sep*oorb
         if(ospin(1).lt.0.d0) ospin(1) = oorb
         if(ospin(2).lt.0.d0) ospin(2) = oorb
      endif
*
      do 500 , k = kmin,kmax
         age = tphys - epoch(k)
         mc = massc(k)
         rc = radc(k)
         CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(k),age,mass(k),tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kstar(k),mc,rc,me,re,k2,bhspin(k),k)
         aj(k) = age
         epoch(k) = tphys - age
         rad(k) = rm
         lumin(k) = lum
         teff(k) = 1000.d0*((1130.d0*lumin(k)/
     &                    (rad(k)**2.d0))**(1.d0/4.d0))
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
         if(tphys.lt.tiny.and.ospin(k).le.0.001d0)then
            ospin(k) = 45.35d0*vrotf(mass(k),ST_tide)/rm
         endif
         jspin(k) = ospin(k)*(k2*rm*rm*(mass(k)-mc)+k3*rc*rc*mc)
         if(.not.sgl)then
            q(k) = mass(k)/mass(3-k)
            rol(k) = rl(q(k))*sep*(1.d0-ecc)
         endif
         rol0(k) = rol(k)
         dmr(k) = 0.d0
         dmt(k) = 0.d0
         djspint(k) = 0.d0
         dtmi(k) = 1.0d+06
         B(k) = 0.d0 !PK
         if(kstar(k).ne.13)then
            bacc(k) = 0.d0
            tacc(k) = 0.d0
         endif
*
 500  continue
*
      if(output) write(*,*)'Init:',mass(1),mass(2),massc(1),massc(2),
     & rad(1),rad(2),kstar(1),kstar(2),sep,ospin(1),ospin(2),jspin(1),
     & jspin(2),sigma,eddfac,z,id1_pass,id2_pass,tphysf,tphys,iter,tsave
*
      if(mt2.lt.tiny)then
         sep = 0.d0
         if(kst.gt.0)then
            mass0(kst) = 0.d0
            mass(kst) = 0.d0
            kmin = 3 - kst
            kmax = kmin
         endif
      endif
*
* On the first entry the previous timestep is zero to prevent mass loss.
*
      dtm = 0.d0
      delet = 0.d0
      djorb = 0.d0
*
* Setup variables which control the output (if it is required).
*
      ip = 0
      jp = 0

      dtp_original = dtp

      tsave = tphys
      isave = .true.
      iplot = .false.
      if(dtp.le.0.d0)then
         iplot = .true.
         isave = .false.
         tsave = tphysf
      elseif(dtp.gt.tphysf)then
         isave = .false.
         tsave = tphysf
      endif
      if(tphys.ge.tphysf) goto 140

 4    iter = 0
      intpol = 0
      inttry = .false.
      change = .false.
      prec = .false.
      snova = .false.
      coel = .false.
      com = .false.
      inspiral = .false.
      tphys0 = tphys
      ecc1 = ecc
      j1 = 1
      j2 = 2
      if(kstar(1).ge.10.and.kstar(1).le.14) dtmi(1) = 0.01d0
      if(kstar(2).ge.10.and.kstar(2).le.14) dtmi(2) = 0.01d0
      dm1 = 0.d0
      dm2 = 0.d0
*
 5    kw1 = kstar(1)
      kw2 = kstar(2)

      if(xi.eq.0)then
         ospin(1) = 1.0d-10
         ospin(2) = 1.0d-10
      endif

*
      dt = 1.0d+06*dtm
      eqspin = 0.d0
      djtt = 0.d0
*
      if(output) write(*,*)'1st in 5: ',tphys,dt,kw1,kw2,
     & mass(1),mass(2),intpol,iter
*
      if(intpol.eq.0.and.ABS(dtm).gt.tiny.and..not.sgl)then
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 501 , k = 1,2
*
* Determine the eddington limit for the accretor (3-k)
* Just in case the wind mass loss rates are *very* high
*
            dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(3-k)
* For BHs, follow Marchant et al. 2017
            if(kstar(3-k).eq.14)then
               maxspinBH = 6.d0**(1.d0/2.d0) * Mbh_initial
               if(mass(3-k).lt.maxspinBH)then
                  etaBH = 1.d0 -
     &    (1.d0 - (mass(3-k)/(3.d0*Mbh_initial))**(2.d0))**(1.d0/2.d0)
               else
                  etaBH = 0.42
               endif
               dme = 1.04e-3*eddfac*(1.d0/(1.d0 + zpars(11)))
     &    *(1.d0/etaBH)*rad(3-k)
            endif

*
* Calculate wind mass loss from the previous timestep.
*
            !check for kstar added by PA
            if(neta.gt.tiny .and. kstar(k)<15)then
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
*
* Calculate how much of wind mass loss from companion will be
* accreted (Boffin & Jorissen, A&A 1988, 205, 155).
*
               if(beta.lt.0.d0)then !PK. following startrack
                  beta = 0.125
                  if(kstar(k).le.1)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.1.4d0)then
                        beta = 0.5
                     else
                        beta = 7.d0*((mass(k)-1.4d0)/(120.d0-1.4d0))
     &                         + 0.5d0
                     endif
                  elseif(kstar(k).ge.7.and.kstar(k).le.9)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.10.d0)then
                        beta = 0.125
                     else
                        beta = 7.d0*((mass(k)-10.d0)/(120.d0-10.d0))
     &                               + 0.125d0
                     endif
                  endif
               endif
               vwind2 = 2.d0*beta*acc1*mass(k)/rad(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),0.8d0*dmr(k))
*
* Apply Eddington limit just in case
*
               if(dt.gt.0)then
                  dmt(3-k) = MIN(dmt(3-k),dme)
               endif
               beta = betahold
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 501     continue

*
* Diagnostic for Symbiotic-type stars.
*
         !check for kstar added by PA
         if(neta.gt.tiny .and. kstar(j2)<15)then
            lacc = 3.14d+07*mass(j2)*dmt(j2)/rad(j2)
            lacc = lacc/lumin(j1)
         endif
*
* Calculate orbital angular momentum change due to wind mass loss.
*
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) +
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))*
     &           sep*sep*sqome2*oorb/(mass(1)+mass(2))**2
         delet = ecc*(dmt(1)*(0.5d0/mass(1) + 1.d0/(mass(1)+mass(2))) +
     &                dmt(2)*(0.5d0/mass(2) + 1.d0/(mass(1)+mass(2))))
*
* For very close systems include angular momentum loss owing to
* gravitational radiation.
*
         if(sep.le.100000.d0.and.grflag.eq.1)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr
            delet = delet + delet1
         endif
*
         do 502 , k = 1,2

* Evaluate convective/radiative limits for a variety of stars as based
* on the work of Belczynski et al. (2008). As option.
* Note only certain stars of type k = 0, 1, 2, 4, 5, 6, 9 **double check
* 0, 1 have range of 0.35-Mms,conv. Mms,conv is function of metalicity.
* 2 & 4 have a temperature dependence. Convective if Teff < 10**3.73
* 3, 5 & 6 are giants and have no radiative envelopes.
* 9's envelope is convective when M < Mhe,conv, Mhe,conv = 3.0Msun.
*
            if(ST_cr.le.0)then
               if(kstar(k).le.1)then
                  convradcomp(k) = 1.25d0
               else
                  convradcomp(k) = 99999999.d0
               endif
            else
               if(kstar(k).le.1)then
*                 Main sequence mass limit varying with metallicity
                  convradcomp(k) = 1.25d0
                  if(z.gt.0.001d0.and.z.lt.zsun)then
                     convradcomp(k) = 0.747d0 + 55.73d0*z - 1532*z*z
                  elseif(z.le.0.001d0)then
                     convradcomp(k) = 0.8d0
                  endif
               elseif(kstar(k).eq.2.or.kstar(k).eq.4)then
*                 H-rich HG and CHeB temperature limit
                  convradcomp(k) = 10.d0**3.73d0
               elseif(kstar(k).le.6)then
*                 No limit or all other giant stars
*                 (compare this large value to mass).
                  convradcomp(k) = 99999999.d0
               elseif(kstar(k).eq.9)then
*                 Mass limit for evolved He star
                  convradcomp(k) = 3.d0
               endif
            endif
*
* Calculate change in the intrinsic spin of the star.
*
*
* Pulsar wind spinup. PK.
*
            if(kstar(k).eq.13.and.pulsar.gt.0)then
* Place propeller stuff here eventually.
               if(xip.eq.1)then
* Modify spin angular momentum accretion according to Ruffert (1999)
* or magnetic field strength limit tested (but by no means established - in fact it
* is likely to change, especially when propeller and ablation are implemented)
* in Kiel et al. (2008).
                  if(B(k).gt.0.d0)then
                     xi = MIN(1.d0,(0.01d0*(2.0d+11/B(k)))+0.01d0)
                     xi = 0.01
                  else
                     xi = 0.01
                  endif
               endif
               djtx(k) = (2.d0/3.d0)*xi*dmt(k)*rad(3-k)*rad(3-k)*
     &                   ospin(3-k)
               djspint(k) = (2.d0/3.d0)*(dmr(k)*rad(k)*rad(k)*ospin(k))-
     &                      djtx(k)
               xi = xihold
               if(output) write(*,*)'502 1:',k,djspint(k),djtx(k),
     & dmt(k),rad(3-k),ospin(3-k),xi,xihold
            else
               djtx(k) = (2.d0/3.d0)*xi*dmt(k)*
     &                     rad(3-k)*rad(3-k)*ospin(3-k)
               djspint(k) = (2.d0/3.d0)*
     &                       (dmr(k)*rad(k)*rad(k)*ospin(k)) - djtx(k)
            endif
*
* Include magnetic braking for stars that have appreciable convective
* envelopes. This includes MS stars with M < 1.25, HG stars near the GB
* and giants. MB is not allowed for fully convective MS stars.
*
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
*               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
*               djspint(k) = djspint(k) + djmb
*
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
*     &              menv(k).gt.0.0d0)then
*              if (ospin(k) .le. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &              (ospin(k)/wsun)**3.0d0
*              if (ospin(k) .gt. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &             (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
*              djspint(k) = djspint(k) + djmb
            djmb = 0.d0
            if(htpmb.eq.0)then
* HTP02 method
               if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
                  djmb = 5.83d-16*menv(k)*
     &                   (rad(k)*ospin(k))**3/mass(k)
                  djspint(k) = djspint(k) + djmb
               endif
            else
               if(ST_cr.le.0.and.mass(k).gt.0.35d0.and.
     &            kstar(k).lt.10.and.menv(k).gt.0.0d0)then
* Ivanova & Taam (2002) method
                  if(ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                (ospin(k)/wsun)**3.0d0
                  if(ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                  djspint(k) = djspint(k) + djmb
               elseif(ST_cr.gt.0.and.menv(k).gt.0.d0.and.
     &            ((kstar(k).le.1.and.mass(k).gt.0.35d0.and.
     &              mass(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.2.and.teff(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.4.and.teff(k).le.convradcomp(k)).or.
     &              ((kstar(k).eq.3).or.(kstar(k).eq.5).or.
     &              (kstar(k).eq.6))))then
* Ivanova & Taam (2002) method
                  if(ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &               (ospin(k)/wsun)**3.0d0
                  if(ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &               (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                  djspint(k) = djspint(k) + djmb
               endif
            endif
* Limit to a 3% angular momentum change for the star owing to MB.
* This is found to work best with the maximum iteration of 20000,
* i.e. does not create an excessive number of iterations, while not
* affecting the evolution outcome when compared with a 2% restriction.
*
            if(djmb.gt.tiny)then
               dtj = 0.03d0*jspin(k)/ABS(djmb)
               dt = MIN(dt,dtj)
               if(output) write(*,*)'mb1:',tphys,dt,djmb,djt
            endif
*
            if(kstar(k).eq.13.and.pulsar.gt.0)then
*
* NS(pulsar) magnetic braking. PK.
* No contact binary NS system.
*
               if(bdecayfac.eq.0)then
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny)then
                     B(k) = B_0(k)*exp(-CK*bacc(k)) + Bbot
                  else
                     B(k) = B_0(k)*
     &                       EXP(-(tphys-epoch(k)-tacc(k))/bconst)*
     &                       exp(-CK*bacc(k)) + Bbot
                  endif
               else
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny.and.
     &                bacc(k).eq.0.d0)then
                     B(k) = B_0(k) + Bbot
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny)then
                     B(k) = B_0(k)/(1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  elseif(bacc(k).eq.0.d0)then
                     B(k) = B_0(k)*EXP(-(tphys-epoch(k)-tacc(k))/bconst)
     &                       + Bbot
                  else
                     B(k) = B_0(k)*
     &                       EXP(-(tphys-epoch(k)-tacc(k))/bconst)/
     &                       (1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  endif
               endif
               omdot = Kconst*B(k)*B(k)*ospin(k)**3
               djmb = 0.4d0*mass(k)*rad(k)*rad(k)*omdot
               djspint(k) = djspint(k) + djmb
               if(output) write(*,*)'502 2:',k,djspint(k),djmb
               if(djmb.gt.tiny)then
                  dtj = 0.1d0*(jspin(k)/ABS(djmb))
                  dt = MIN(dt,dtj)
               endif
            endif
*
* Calculate circularization, orbital shrinkage and spin up.
*
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (rad(k)/sep)**2
               raa6 = raa2**3
*
* Hut's polynomials.
*
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if(ST_cr.le.0.and.
     &            ((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7))then
*
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).
*
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*rad(k)*rad(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(ST_cr.gt.0.and.
     &                ((kstar(k).le.1.and.mass(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.2.and.teff(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.4.and.teff(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.9.and.mass(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.7).or.(kstar(k).eq.8)))then
*
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).
*
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*rad(k)*rad(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k) !for startrack comparison this should be =0.31622777d0 (sqrt(0.1)).
               elseif(kstar(k).le.9)then
*
* Convective damping (Hut, 1981, A&A, 99, 126).
*
* In BSE paper Equation 30, the default scaling coefficient is 2./21
* the fprimc_array kstar dependent array that is fed in
* keeps this same coefficient by default but allows user to
* specify their own
*
                  tc = mr23yr*(menv(k)*renv(k)*(rad(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
                  f = MIN(1.d0,(ttid/(2.d0*tc))**2)
                  !kstar(k) to kstar(k)+1 by PA for kstar(k)=0
                  tcqr = fprimc_array(kstar(k)+1)*
     &                 f*q(3-k)*raa6*menv(k)/
     &                 (tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               elseif(ST_tide.le.0)then
*
* Degenerate damping (Campbell, 1984, MNRAS, 207, 433)
*
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
*
* Circularization.
*
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1
*
* Spin up of star.
*
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
               if(output) write(*,*)'502 3:',k,dspint(k),tcqr
*
* Calculate the equilibrium spin at which no angular momentum
* can be transferred.
*
               eqspin = oorb*f2/(sqome3*f5)
*
* Calculate angular momentum change for the star owing to tides.
*
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               if(kstar(k).le.6.or.ABS(djt)/jspin(k).gt.0.1d0)then
                  djtt = djtt + djt
               endif
            endif
 502     continue
*
* Limit to 2% orbital angular momentum change.
*
         djtt = djtt + djorb
         if(ABS(djtt).gt.tiny)then
            dtj = 0.02d0*jorb/ABS(djtt)
            dt = MIN(dt,dtj)
         endif
         dtm = dt/1.0d+06
         if(output) write(*,*)'bin lim orb ang mom:',tphys,dt,kstar(1)
*
      elseif(ABS(dtm).gt.tiny.and.sgl)then
         do 503 , k = kmin,kmax
            !check for kstar added by PA
            if(neta.gt.tiny .and. kstar(k)<15)then
               rlperi = 0.d0
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
            else
               dmr(k) = 0.d0
            endif
            dmt(k) = 0.d0
            djspint(k) = (2.d0/3.d0)*dmr(k)*rad(k)*rad(k)*ospin(k)
            if(output) write(*,*)'503 1:',k,djspint(k)
*
* Evaluate convective/radiative limits for a variety of stars as based
* on the work of Belczynski et al. (2008). As option.
* Note only certain stars of type k = 0, 1, 2, 4, 5, 6, 9 **double check
* 0, 1 have range of 0.35-Mms,conv. Mms,conv is function of metalicity.
* 2 & 4 have a temperature dependence. Convective if Teff < 10**3.73
* 3, 5 & 6 are giants and have no radiative envelopes.
* 9's envelope is convective when M < Mhe,conv, Mhe,conv = 3.0Msun.
*
            if(ST_cr.le.0)then
               if(kstar(k).le.1)then
                  convradcomp(k) = 1.25d0
               else
                  convradcomp(k) = 99999999.d0
               endif
            else
               if(kstar(k).le.1)then
*                 Main sequence mass limit varying with metallicity
                  convradcomp(k) = 1.25d0
                  if(z.gt.0.001d0.and.z.lt.zsun)then
                     convradcomp(k) = 0.747d0 + 55.73d0*z - 1532*z*z
                  elseif(z.le.0.001d0)then
                     convradcomp(k) = 0.8d0
                  endif
               elseif(kstar(k).eq.2.or.kstar(k).eq.4)then
*                 H-rich HG and CHeB temperature limit
                  convradcomp(k) = 10.d0**3.73d0
               elseif(kstar(k).le.6)then
*                 No limit or all other giant stars
*                 (compare this large value to mass).
                  convradcomp(k) = 99999999.d0
               elseif(kstar(k).eq.9)then
*                 Mass limit for evolved He star
                  convradcomp(k) = 3.d0
               endif
            endif
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
*               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
*               djspint(k) = djspint(k) + djmb
*             if(mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
*     &              menv(k).gt.0.0d0)then
*                if (ospin(k) .le. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &              (ospin(k)/wsun)**3.0d0
*                if (ospin(k) .gt. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &             (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
*                djspint(k) = djspint(k) + djmb
*
*MB
            djmb = 0.d0
            if(htpmb.eq.0)then
               if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
                  djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
                  djspint(k) = djspint(k) + djmb
               endif
            elseif(htpmb.eq.1)then
               if(ST_cr.le.0.and.
     &            mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
     &            menv(k).gt.0.0d0)then
* Old BSE convective/radiative divide...
                   if(ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                 (ospin(k)/wsun)**3.0d0
                   if(ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &             (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                   djspint(k) = djspint(k) + djmb
               elseif(ST_cr.gt.0.and.(menv(k).gt.0.d0.and.
     &            ((kstar(k).le.1.and.mass(k).gt.0.35d0.and.
     &              mass(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.2.and.teff(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.4.and.teff(k).le.convradcomp(k)).or.
     &              ((kstar(k).eq.3).or.(kstar(k).eq.5).or.
     &              (kstar(k).eq.6)))))then
* MB given in Ivanova & Taam (2002)
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
*     &              menv(k).gt.0.0d0)then
                  if(ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                (ospin(k)/wsun)**3.0d0
                  if(ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &               (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                  djspint(k) = djspint(k) + djmb
               endif
            endif
*                if(output) write(*,*)'503 2:',k,djspint(k),djmb
            if(djmb.gt.tiny)then
               dtj = 0.03d0*jspin(k)/ABS(djmb)
               dt = MIN(dt,dtj)
            endif
            if(kstar(k).eq.13.and.pulsar.gt.0)then
*
* NS magnetic braking. PK.
* Single NS evolution.
*
               if(bdecayfac.eq.0)then
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny)then
                     B(k) = B_0(k)*exp(-CK*bacc(k)) + Bbot
                  else
                     B(k) = B_0(k)*
     &                       EXP(-(tphys-epoch(k)-tacc(k))/bconst)*
     &                       exp(-CK*bacc(k)) + Bbot
                  endif
               else
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny.and.
     &                bacc(k).eq.0.d0)then
                     B(k) = B_0(k) + Bbot
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny)then
                     B(k) = B_0(k)/(1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  elseif(bacc(k).eq.0.d0)then
                     B(k) = B_0(k)*EXP(-(tphys-epoch(k)-tacc(k))/bconst)
     &                       + Bbot
                  else
                     B(k) = B_0(k)*
     &                       EXP(-(tphys-epoch(k)-tacc(k))/bconst)/
     &                       (1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  endif
               endif
               omdot = Kconst*B(k)*B(k)*ospin(k)**3
               djmb = 0.4d0*mass(k)*rad(k)*rad(k)*omdot
               djspint(k) = djspint(k) + djmb
               if(output) write(*,*)'503 2:',k,djmb,djspint(k)
* Consider update of time-stepping due to dj, i.e. dt = dj/(dj/dt).
               if(djmb.gt.tiny)then
                  dtj = 0.1d0*(jspin(k)/ABS(djmb))
                  dt = MIN(dt,dtj)
               endif
            endif
 503     continue
         dtm = dt/1.0d+06
      endif
*
      do 504 , k = kmin,kmax
*
         dms(k) = (dmr(k) - dmt(k))*dt
         if(kstar(k).lt.10)then
            dml = mass(k) - massc(k)
            if(dml.lt.dms(k))then
               dml = MAX(dml,2.d0*tiny)
               dtm = (dml/dms(k))*dtm
               if(k.eq.2) dms(1) = dms(1)*dml/dms(2)
               dms(k) = dml
               dt = 1.0d+06*dtm
            endif
*
* Limit to 1% mass loss.
*
            if(dms(k).gt.0.01d0*mass(k))then
               dtm = 0.01d0*mass(k)*dtm/dms(k)
               if(k.eq.2) dms(1) = dms(1)*0.01d0*mass(2)/dms(2)
               dms(k) = 0.01d0*mass(k)
               dt = 1.0d+06*dtm
            endif
         endif
         if(output) write(*,*)'after 1p ml: ',tphys,k,dt,dtm
*
 504  continue
*
* Update mass and intrinsic spin (checking that the star is not spun
* past the equilibrium) and reset epoch for a MS (and possibly a HG) star.
*
      do 505 , k = kmin,kmax
*
         if(eqspin.gt.0.d0.and.ABS(dspint(k)).gt.tiny)then
            if(intpol.eq.0)then
               if(dspint(k).ge.0.d0)then
                  dspint(k) = MIN(dspint(k),(eqspin-ospin(k))/dt)
               else
                  dspint(k) = MAX(dspint(k),(eqspin-ospin(k))/dt)
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt
               djspint(k) = djspint(k) - djt
               if(output) write(*,*)'505: ',k,djt,djspint(k),jspin(k),dt
            endif
         endif
*
         if(output) write(*,*)'505 1:',tphys,k,kstar(k),djspint(k),
     & djspint(k)*dt,jspin(k),intpol
         jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k)*dt)
         if(output) write(*,*)'505 2:',tphys,k,kstar(k),djspint(k),
     & djspint(k)*dt,jspin(k)
*
* Ensure that the star does not spin up beyond break-up.
*
         ospbru = twopi*SQRT(mass(k)*aursun**3/rad(k)**3)
         jspbru = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &             k3*massc(k)*radc(k)*radc(k))*ospbru
         if((jspin(k).gt.jspbru.or.
     &      (jspin(k).eq.1.0d-10.and.intpol.gt.0)).and.
     &      ABS(dtm).gt.tiny)then !PDK add check for jspin intpol issues here.
* If rapidly spinning star (generally NS) pushed over the edge (spins up so
* much so that jspin would have become negative according to bse, so it sets
* it to 1.d0-10) in intpol then give it the stars maximum spin.
            mew = 1.d0
            if(djtx(k).gt.0.d0.and.jspin(k).gt.1.d0-10)then !dont go in here if jspin hit wall.
               mew = MIN(mew,(jspin(k) - jspbru)/djtx(k))
            endif
            jspin(k) = jspbru
* If excess material should not be accreted, activate next line.
*           dms(k) = dms(k) + (1.d0 - mew)*dmt(k)*dt
         endif
* Update mass
         if(ABS(dms(k)).gt.tiny)then
            mass(k) = mass(k) - dms(k)
            if(kstar(k).le.2.or.kstar(k).eq.7)then
               m0 = mass0(k)
               mass0(k) = mass(k)
               CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,
     &                   lums,GB,zpars)
               if(kstar(k).eq.2)then
                  if(GB(9).lt.massc(k).or.m0.gt.zpars(3))then
                     mass0(k) = m0
                  else
                     epoch(k) = tm + (tscls(1) - tm)*(aj(k)-tms(k))/
     &                               (tbgb(k) - tms(k))
                     epoch(k) = tphys - epoch(k)
                  endif
               else
                  epoch(k) = tphys - aj(k)*tm/tms(k)
               endif
            endif
* Update NS magnetic field owing to accretion, as a function of mass accreted. PK.
            if(kstar(k).eq.13.and.pulsar.gt.0)then
               if(dms(k).lt.0.d0)then !negative dms is mass gained.
* When propeller ev. include .not.prop here...
                  b_mdot = dms(k)/dt
                  if(b_mdot_lim.gt.0.d0.and.b_mdot.gt.b_mdot_lim)then
                     bacc(k) = bacc(k) - dms(k)
                     tacc(k) = tacc(k) + dtm
                  elseif(b_mdot_lim.le.0.d0)then
                     bacc(k) = bacc(k) - dms(k)
                     tacc(k) = tacc(k) + dtm
                  endif
               endif
            endif
         endif
*
 505  continue
*
      if(.not.sgl)then
*
         ecc1 = ecc1 - delet*dt
         ecc = MAX(ecc1,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
*
         jorb = jorb - djorb*dt
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      endif
*
* Advance the time.
*
      if(intpol.eq.0)then
         tphys0 = tphys
         dtm0 = dtm
      endif
      tphys = tphys + dtm

      if(output) write(*,*)'time upd:',tsave, tphys,dtm,
     &                      kstar(1),kstar(2),ecc
*
      do 6 , k = kmin,kmax
*
* Acquire stellar parameters (M, R, L, Mc & K*) at apparent evolution age.
*
         age = tphys - epoch(k)
         aj0(k) = age
         kw = kstar(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
         if(intpol.eq.0) mcxx(k) = mc
         if(intpol.gt.0) mc = mcxx(k)
         mass00(k) = m0
*
* Masses over 100Msun should probably not be trusted in the
* evolution formulae.
*
         if(mt.gt.100.d0)then
*            WRITE(99,*)' MASS EXCEEDED ',mass1i,mass2i,tbi,ecci,mt,
*     & tphysf,id1_pass,id2_pass
*            goto 140
         endif
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2,bhspin(k),k)
*
         if(kw.ne.15)then
            ospin(k) = jspin(k)/(k2*(mt-mc)*rm*rm+k3*mc*rc*rc)
         endif
*
* At this point there may have been a supernova.
*
         if((kw.ne.kstar(k).and.kstar(k).le.12.and.
     &      (kw.eq.13.or.kw.eq.14)).or.(ABS(merger).ge.20.d0))then
            if(formation(k).ne.11) formation(k) = 1
            if(kw.eq.13.and.ecsn.gt.0.d0)then
               if(kstar(k).le.6)then
                  if(mass0(k).le.zpars(5))then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation(k) = 2
                  endif
               elseif(kstar(k).ge.7.and.kstar(k).le.9)then
                  if(mass(k).gt.ecsn_mlow.and.mass(k).le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation(k) = 2
                  endif
               elseif(formation(k).eq.11)then
* MIC
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 5
               elseif(kstar(k).ge.10.or.kstar(k).eq.12)then
* AIC formation, will never happen here but...
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 4
               elseif(merger.ge.20.d0)then
                  sigma = merger
                  fallback = 0.d0
                  if(merger.ge.200.d0)then!estimate CC SN
                     formation(k) = 1
                  else
                     formation(k) = 5
                  endif
               elseif(merger.le.-20.d0)then
                  sigma = ABS(merger)
                  fallback = 0.d0
                  if(merger.ge.200.d0)then!estimate CC SN
*Sourav:Possible bug in the line above. merger should really be sigms!!
*                  if(sigma.ge.200.d0)then!estimate CC SN
                     formation(k) = 1
                  else
                     formation(k) = 5
                  endif
               endif
            elseif(kw.eq.13.and.aic.gt.0)then
               if(formation(k).eq.11)then
* MIC
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 5
               elseif(kstar(k).ge.10.or.kstar(k).eq.12)then
* AIC formation, will never happen here but...
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 4
               endif
            elseif(ABS(merger).ge.20.d0)then
               sigma = ABS(merger)
               fallback = 0.d0
               if(merger.ge.200.d0)then!estimate CC SN
                  formation(k) = 1
               else
                  formation(k) = 5
               endif
            endif
            if(sgl)then
               evolve_type = 14.d0 + FLOAT(k)
               teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
               teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
               if(B_0(1).eq.0.d0)then !PK.
                  b01_bcm = 0.d0
               elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                  b01_bcm = B_0(1)
               else
                  b01_bcm = B(1)
               endif
               if(B_0(2).eq.0.d0)then
                  b02_bcm = 0.d0
               elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                  b02_bcm = B_0(2)
               else
                  b02_bcm = B(2)
               endif
               CALL writetab(jp,tphys,evolve_type,
     &                      mass(1),mass(2),kstar(1),kstar(2),
     &                      sep,tb,ecc,rrl1,rrl2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      teff1,teff2,radc(1),radc(2),
     &                      menv(1),menv(2),renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin(1),bhspin(2),
     &                      deltam1_bcm,deltam2_bcm,formation(1),
     &                      formation(2),binstate,mergertype,'bpp')
               CALL kick(kw,mass(k),mt,0.d0,0.d0,-1.d0,0.d0,vk,k,
     &                  0.d0,fallback,sigmahold,kick_info,disrupt,bkick)

               sigma = sigmahold !reset sigma after possible ECSN kick dist. Remove this if u want some kick link to the intial pulsar values...
* set kick values for the bcm array

            else
               evolve_type = 14.d0 + FLOAT(k)
               teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
               teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
               if(B_0(1).eq.0.d0)then !PK.
                  b01_bcm = 0.d0
               elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                  b01_bcm = B_0(1)
               else
                  b01_bcm = B(1)
               endif
               if(B_0(2).eq.0.d0)then
                  b02_bcm = 0.d0
               elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                  b02_bcm = B_0(2)
               else
                  b02_bcm = B(2)
               endif

               CALL writetab(jp,tphys,evolve_type,
     &                       mass(1),mass(2),kstar(1),kstar(2),
     &                       sep,tb,ecc,rrl1,rrl2,
     &                       aj(1),aj(2),tms(1),tms(2),
     &                       massc(1),massc(2),rad(1),rad(2),
     &                       mass0(1),mass0(2),lumin(1),lumin(2),
     &                       teff1,teff2,radc(1),radc(2),
     &                       menv(1),menv(2),renv(1),renv(2),
     &                       ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                       bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                       epoch(2),bhspin(1),bhspin(2),
     &                       deltam1_bcm,deltam2_bcm,formation(1),
     &                       formation(2),binstate,mergertype,'bpp')

               CALL kick(kw,mass(k),mt,mass(3-k),ecc,sep,jorb,vk,k,
     &              rad(3-k),fallback,sigmahold,kick_info,disrupt,bkick)
               sigma = sigmahold !reset sigma after possible ECSN kick dist. Remove this if u want some kick link to the intial pulsar values...
* set kick values for the bcm array
               if(mass(3-k).lt.0.d0)then
                  if(kstar(3-k).lt.0.d0) mt = mt-mass(3-k) !ignore TZ object
                  if(kw.eq.13.and.mt.gt.mxns) kw = 14
                  CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
                  kstar(k) = kw
                  mass(k) = mt
                  epoch(k) = tphys - age
                  kstar(3-k) = 15
                  mass(3-k) = 0.d0
                  coel = .true.
                  binstate = 1
                  goto 135
               endif
               if(ecc.gt.1.d0)then
                  kstar(k) = kw
                  mass(k) = mt
                  epoch(k) = tphys - age
                  goto 135
               endif
               tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
               oorb = twopi/tb
            endif

            merger = -1.d0
            snova = .true.
         endif
*
         if(kw.ne.kstar(k))then
            change = .true.
            mass(k) = mt
            dtmi(k) = 0.01d0
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
         endif

*
*
* Force new NS or BH to have a birth spin peirod and magnetic field.
*
         if(kstar(k).eq.13.or.kstar(k).eq.14)then
* re-calculate the Eddington limit of the newly formed CO
            dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(k)
* For BHs, follow Marchant et al. 2017
            if(kstar(k).eq.14)then
               maxspinBH = 6.d0**(1.d0/2.d0) * Mbh_initial
               if(mass(k).lt.maxspinBH)then
                  etaBH = 1.d0 -
     &    (1.d0 - (mass(k)/(3.d0*Mbh_initial))**(2.d0))**(1.d0/2.d0)
               else
                  etaBH = 0.42
               endif
               dme = 1.04e-3*eddfac*(1.d0/(1.d0 + zpars(11)))
     &    *(1.d0/etaBH)*rad(k)
            endif
            if(tphys-epoch(k).lt.tiny)then
               if(kstar(k).eq.13.and.pulsar.gt.0)then
*                  write(93,*)'birth start: ',tphys,k,B_0(k),ospin(k)
*                  CALL FLUSH(93)
 170              u1 = ran3(idum1)
                  if(u1.ge.1.d0) goto 170
                  u2 = ran3(idum1)
                  s = sqrt(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
                  s = 0.7d0*s - 0.6d0
                  if(s.ge.0.013d0.or.s.le.-1.5d0) goto 170
                  ospin(k) = (twopi*yearsc)/(10.d0**s)
 174              u1 = ran3(idum1)
                  if(u1.ge.1.d0) goto 174
                  u2 = ran3(idum1)
                  s = sqrt(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
                  s = 0.68d0*s + 12.6d0
                  if(s.lt.11.5d0.or.s.gt.13.8d0) goto 174
                  B_0(k) = 10.d0**s
                  bacc(k) = 0.d0 ! If it has been a NS before reset
                  tacc(k) = 0.d0
                  if((merger.le.-2.d0.and.merger.gt.-20.d0).or.
     &                merge_mem.eq.1)then
* Reset as MSP.
 175                 u1 = ran3(idum1)
                     u2 = ran3(idum1)
                     if(u1.gt.0.9999d0) u1 = 0.9999d0
                     if(u2.gt.1.d0) u2 = 1.d0
                     s = SQRT(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
*                  s = 0.7d0*s - 0.6d0
                     s = 0.5d0*s - 2.25d0
*                  if(s.ge.0.013.or.s.le.-1.5) goto 173
*                  if(s.ge.-2.0457d0.or.s.le.-2.53d0) goto 175
                     if(s.ge.-1.6457d0.or.s.le.-2.53d0) goto 175
                     ospin(k) = (twopi*yearsc)/(10.d0**s)!have commented this out to keeps same spin
*                  write(*,*)'P=',s
 176                 u1 = ran3(idum1)
                     u2 = ran3(idum1)
                     if(u1.gt.0.9999d0) u1 = 0.9999d0
                     if(u2.gt.1.d0) u2 = 1.d0
                     s = SQRT(-2.d0*LOG(1.d0-u1))*COS(twopi*u2)
                     s = 0.3d0*s + 8.38d0
                     if(s.lt.8.d0.or.s.gt.8.778d0) goto 176
                     B_0(k) = 10.d0**s
                  endif
               else
                  ospin(k) = 2.0d+08
               endif
               bacc(k) = 0.d0 !make sure if its been a NS before its now a new one...
               tacc(k) = 0.d0
               jspin(k) = k3*rc*rc*mc*ospin(k)
               if(output) write(*,*)'SN: ',k,kstar(k),ospin(k),sigma
               sigma = sigmahold !reset sigma after possible ECSN kick dist.
            endif
         endif
         merge_mem = 0
*
* Set radius derivative for later interpolation.
*
         if(ABS(dtm).gt.tiny)then
            rdot(k) = ABS(rm - rad(k))/dtm
         else
            rdot(k) = 0.d0
         endif
*
*     Base new time scale for changes in radius & mass on stellar type.
*
         dt = dtmi(k)
         CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
         if(output) write(*,*)'post deltat:',tphys,dt,dtr,kw,
     & age,intpol,iter,k,kmin,kmax
*
* Choose minimum of time-scale and remaining interval.
*
         dtmi(k) = MIN(dt,dtr)
*
* Save relevent solar quantities.
*
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         lumin(k) = lum
         teff(k) = 1000.d0*((1130.d0*lumin(k)/
     &                    (rad(k)**2.d0))**(1.d0/4.d0))
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
 6    continue
*
      if(.not.sgl)then
*
* Determine the mass ratios.
*
         do 506 , k = 1,2
            q(k) = mass(k)/mass(3-k)
 506     continue
*
* Determine the Roche lobe radii and adjust the radius derivative.
*
         do 507 , k = 1,2
            rol(k) = rl(q(k))*sep*(1.d0-ecc)
            if(ABS(dtm).gt.tiny)then
               rdot(k) = rdot(k) + (rol(k) - rol0(k))/dtm
               rol0(k) = rol(k)
            endif
 507     continue
      else
         do 508 , k = kmin,kmax
            rol(k) = 10000.d0*rad(k)
 508     continue
      endif
*
      if((tphys.lt.tiny.and.ABS(dtm).lt.tiny.and.
     &    (mass2i.lt.0.1d0.or..not.sgl)).or.snova)then
          evolve_type = 1.d0
          rrl1 = rad(1)/rol(1)
          rrl2 = rad(2)/rol(2)
          teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
          teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
          if(B_0(1).eq.0.d0)then !PK.
             b01_bcm = 0.d0
          elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
             b01_bcm = B_0(1)
          else
             b01_bcm = B(1)
          endif
          if(B_0(2).eq.0.d0)then
             b02_bcm = 0.d0
          elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
             b02_bcm = B_0(2)
          else
             b02_bcm = B(2)
          endif

          CALL writetab(jp,tphys,evolve_type,
     &                  mass(1),mass(2),kstar(1),kstar(2),sep,
     &                  tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bpp')
         if(snova)then
            bpp(jp,11) = 2.0
            dtm = 0.d0
            goto 4
         endif
      endif
      if(check_dtp.eq.1)then
          CALL checkstate(dtp,dtp_original,tsave,tphys,tphysf,
     &                      iplot,isave,binstate,evolve_type,
     &                      mass(1),mass(2),kstar(1),kstar(2),sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      radc(1),radc(2),menv(1),menv(2),
     &                      renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),
     &                      tacc(1),tacc(2),epoch(1),epoch(2),
     &                      bhspin(1),bhspin(2),teff1,teff2)
      endif
*
      if((isave.and.tphys.ge.tsave).or.iplot)then
         if(sgl.or.(rad(1).lt.rol(1).and.rad(2).lt.rol(2)).
     &      or.tphys.lt.tiny)then
            if(B_0(1).eq.0.d0)then !PK.
               b01_bcm = 0.d0
            elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
               b01_bcm = B_0(1)
            else
               b01_bcm = B(1)
            endif
            if(B_0(2).eq.0.d0)then
               b02_bcm = 0.d0
            elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
               b02_bcm = B_0(2)
            else
               b02_bcm = B(2)
            endif
            teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
            teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
            rrl1 = rad(1)/rol(1)
            rrl2 = rad(2)/rol(2)
            deltam1_bcm = dmt(1) - dmr(1)
            deltam2_bcm = dmt(2) - dmr(2)
* Check if PISN occurred, and if so overwrite formation
            if(pisn_track(1).ne.0) formation(1) = pisn_track(1)
            if(pisn_track(2).ne.0) formation(2) = pisn_track(2)

            CALL writetab(ip,tphys,evolve_type,
     &                    mass(1),mass(2),kstar(1),kstar(2),
     &                    sep,tb,ecc,rrl1,rrl2,
     &                    aj(1),aj(2),tms(1),tms(2),
     &                    massc(1),massc(2),rad(1),rad(2),
     &                    mass0(1),mass0(2),lumin(1),lumin(2),
     &                    teff1,teff2,radc(1),radc(2),
     &                    menv(1),menv(2),renv(1),renv(2),
     &                    ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                    bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                    epoch(2),bhspin(1),bhspin(2),
     &                    deltam1_bcm,deltam2_bcm,formation(1),
     &                    formation(2),binstate,mergertype,'bcm')
            if(isave) tsave = tsave + dtp
            if(output) write(*,*)'bcm1',kstar(1),kstar(2),mass(1),
     & mass(2),rad(1),rad(2),ospin(1),ospin(2),jspin(1)
*     & mass(2),rad(1),rad(2),ospin(1),ospin(2),b01_bcm,b02_bcm,jspin(1)
         endif
      endif
*
* If not interpolating set the next timestep.
*
      if(intpol.eq.0)then
*         WRITE(*,*)'you should see this to advance the time'
         if(output) write(*,*)'nxt t, prior:',tphys,dtm,dtmi(1),dtmi(2)
         dtm = MAX(1.0d-07*tphys,MIN(dtmi(1),dtmi(2)))
         dtm = MIN(dtm,tsave-tphys)
         if(output) write(*,*)'nxt t, after:',tphys,dtm,dtmi(1),dtmi(2)
         if(iter.eq.0) dtm0 = dtm
      endif
      if(sgl) goto 98
*
* Set j1 to the donor - the primary
* and j2 to the accretor - the secondary.
*
      if(intpol.eq.0)then
         if(rad(1)/rol(1).ge.rad(2)/rol(2))then
            j1 = 1
            j2 = 2
         else
            j1 = 2
            j2 = 1
         endif
      endif
*
* Test whether Roche lobe overflow has begun.
*
      if(rad(j1).gt.rol(j1))then
*
* Interpolate back until the primary is just filling its Roche lobe.
*
         if(rad(j1).ge.1.002d0*rol(j1))then
            if(intpol.eq.0) tphys00 = tphys
            intpol = intpol + 1
            if(iter.eq.0) goto 7
            if(inttry) goto 7
            if(intpol.ge.100)then
               WRITE(99,*)' INTPOL EXCEEDED ',mass1i,mass2i,tbi,ecci
               goto 140
            endif
            dr = rad(j1) - 1.001d0*rol(j1)
            if(ABS(rdot(j1)).lt.tiny.or.prec)then
               goto 7
            endif
            dtm = -dr/ABS(rdot(j1))
            if(ABS(tphys0-tphys).gt.tiny) dtm = MAX(dtm,tphys0-tphys)
            if(kstar(1).ne.kw1)then
               kstar(1) = kw1
               mass0(1) = mass00(1)
               epoch(1) = tphys - aj0(1)
            endif
            if(kstar(2).ne.kw2)then
               kstar(2) = kw2
               mass0(2) = mass00(2)
               epoch(2) = tphys - aj0(2)
            endif
            change = .false.
         else
*
* Enter Roche lobe overflow
*
            if(tphys.ge.tphysf) goto 140
            goto 7
         endif
      else
*
* Check if already interpolating.
*
         if(intpol.gt.0)then
            intpol = intpol + 1
            if(intpol.ge.80)then
               inttry = .true.
            endif
            if(ABS(rdot(j1)).lt.tiny)then
               prec = .true.
               dtm = 1.0d-07*tphys
            else
               dr = rad(j1) - 1.001d0*rol(j1)
               dtm = -dr/ABS(rdot(j1))
            endif
            if((tphys+dtm).ge.tphys00)then
*
* If this occurs then most likely the star is a high mass type 4
* where the radius can change very sharply or possibly there is a
* discontinuity in the radius as a function of time and HRDIAG
* needs to be checked!
*
               dtm = 0.5d0*(tphys00 - tphys0)
               dtm = MAX(dtm,1.0d-10)
               prec = .true.
            endif
            tphys0 = tphys
         endif
      endif
*
* Check for collision at periastron.
*
      pd = sep*(1.d0 - ecc)
      if(pd.lt.(rad(1)+rad(2)).and.intpol.eq.0) goto 130
*
* Go back for the next step or interpolation.
*
 98   continue
      if(tphys.ge.tphysf.and.intpol.eq.0) goto 140
      if(change)then
         change = .false.
         evolve_type = 2.d0
         rrl1 = rad(1)/rol(1)
         rrl2 = rad(2)/rol(2)
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
      endif
*
      iter = iter + 1
*
      if(iter.ge.loop)then
         WRITE(99,*)' MAXIMUM ITER EXCEEDED ',mass1i,mass2i,tbi,ecci,
     & tphysf,id1_pass,id2_pass
         goto 140
      endif
      goto 5
*
* Set the nuclear timescale in years and slow-down factor.
*
 7    km0 = dtm0*1.0d+03/tb
      if(km0.lt.tiny) km0 = 0.5d0
      
* Check for collision at periastron for a stable RLOF 
      pd = sep*(1.d0 - ecc)
      if(pd.lt.(rad(1)+rad(2))) goto 130
      
*      
* Force co-rotation of primary and orbit to ensure that the tides do not
* lead to unstable Roche (not currently used).
*
*     if(ospin(j1).gt.1.05d0*oorb)then
*        ospin(j1) = oorb
*        jspin(j1) = (k2str(j1)*rad(j1)*rad(j1)*(mass(j1)-massc(j1))+
*    &                k3*radc(j1)*radc(j1)*massc(j1))*ospin(j1)
*     endif
*
      sep = sep*(1-ecc)
      ecc = 0.d0
      tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
      oorb = twopi/tb
      jorb = (mass(1)*mass(2)/(mass(1)+mass(2)))*
     &        sqrt(1.d0-ecc*ecc)*sep*sep*oorb
* Note: this will modify pulsar spin period
      if(kstar(1).lt.13)then
         ospin(1) = oorb
         jspin(1) = (k2str(1)*rad(1)*rad(1)*(mass(1)-massc(1))+
     &               k3*radc(1)*radc(1)*massc(1))*ospin(1)
      endif
      if(kstar(2).lt.13)then
         ospin(2) = oorb
         jspin(2) = (k2str(2)*rad(2)*rad(2)*(mass(2)-massc(2))+
     &               k3*radc(2)*radc(2)*massc(2))*ospin(2)
      endif
      km0 = dtm0*1.0d+03/tb
      if(km0.lt.tiny) km0 = 0.5d0
      rol(1) = rl(q(1))*sep !okay like this 'cos sep is peri. dist.
      rol(2) = rl(q(2))*sep
      iter = 0
      coel = .false.
      change = .false.
      radx(j1) = MAX(radc(j1),rol(j1))
      radx(j2) = rad(j2)
*
      evolve_type = 3.0
      rrl1 = rad(1)/rol(1)
      rrl2 = rad(2)/rol(2)
      teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
      teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
      if(B_0(1).eq.0.d0)then !PK.
         b01_bcm = 0.d0
      elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
         b01_bcm = B_0(1)
      else
         b01_bcm = B(1)
      endif
      if(B_0(2).eq.0.d0)then
         b02_bcm = 0.d0
      elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
         b02_bcm = B_0(2)
      else
         b02_bcm = B(2)
      endif

      CALL writetab(jp,tphys,evolve_type,
     &              mass(1),mass(2),kstar(1),kstar(2),sep,
     &              tb,ecc,rrl1,rrl2,
     &              aj(1),aj(2),tms(1),tms(2),
     &              massc(1),massc(2),rad(1),rad(2),
     &              mass0(1),mass0(2),lumin(1),lumin(2),
     &              teff1,teff2,radc(1),radc(2),
     &              menv(1),menv(2),renv(1),renv(2),
     &              ospin(1),ospin(2),b01_bcm,b02_bcm,
     &              bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &              epoch(2),bhspin(1),bhspin(2),
     &              deltam1_bcm,deltam2_bcm,formation(1),
     &              formation(2),binstate,mergertype,'bpp')
*
      if(check_dtp.eq.1)then
          CALL checkstate(dtp,dtp_original,tsave,tphys,tphysf,
     &                      iplot,isave,binstate,evolve_type,
     &                      mass(1),mass(2),kstar(1),kstar(2),sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      radc(1),radc(2),menv(1),menv(2),
     &                      renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),
     &                      tacc(1),tacc(2),epoch(1),epoch(2),
     &                      bhspin(1),bhspin(2),teff1,teff2)
      endif

      if(iplot.and.tphys.gt.tiny)then
          if(B_0(1).eq.0.d0)then !PK.
              b01_bcm = 0.d0
          elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
              b01_bcm = B_0(1)
          else
              b01_bcm = B(1)
          endif
          if(B_0(2).eq.0.d0)then
              b02_bcm = 0.d0
          elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
              b02_bcm = B_0(2)
          else
              b02_bcm = B(2)
          endif
          teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
          teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
          rrl1 = rad(1)/rol(1)
          rrl2 = rad(2)/rol(2)
          deltam1_bcm = 0.0
          deltam2_bcm = 0.0
* Check if PISN occurred, and if so overwrite formation
          if(pisn_track(1).ne.0) formation(1) = pisn_track(1)
          if(pisn_track(2).ne.0) formation(2) = pisn_track(2)
          CALL writetab(ip,tphys,evolve_type,
     &                  mass(1),mass(2),kstar(1),kstar(2),
     &                  sep,tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bcm')
         if(output) write(*,*)'bcm2:',kstar(1),kstar(2),mass(1),
     & mass(2),rad(1),rad(2),ospin(1),ospin(2),jspin(1)
*     & mass(2),rad(1),rad(2),ospin(1),ospin(2),b01_bcm,b02_bcm,jspin(1)
      endif
*
* Eddington limit for accretion on to the secondary in one orbit.
*
 8    dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb
* For BHs, follow Marchant et al. 2017
      if(kstar(j2).eq.14)then
         maxspinBH = 6.d0**(1.d0/2.d0) * Mbh_initial
         if(mass(j2).lt.maxspinBH)then
            etaBH = 1.d0 -
     &    (1.d0 - (mass(j2)/(3.d0*Mbh_initial))**(2.d0))**(1.d0/2.d0)
         else
            etaBH = 0.42
         endif
         dme = 1.04e-3*eddfac*(1.d0/(1.d0 + zpars(11)))
     &    *(1.d0/etaBH)*rad(j2)*tb
      endif

      supedd = .false.
      novae = .false.
      disk = .false.
*
* Determine whether the transferred material forms an accretion
* disk around the secondary or hits the secondary in a direct
* stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509)
* fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383).
*
*     if(kstar(j2).ge.10) disk = .true.
      rmin = 0.0425d0*sep*(q(j2)*(1.d0+q(j2)))**(1.d0/4.d0)
      if(rmin.gt.rad(j2)) disk = .true.
*
* Kelvin-Helmholtz time from the modified classical expression.
*
      do 13 , k = 1,2
         tkh(k) = 3.13d+07*mass(k)/(rad(k)*lumin(k))
         if(kstar(k).le.1.or.kstar(k).eq.7.or.kstar(k).ge.10)then
            tkh(k) = tkh(k)*mass(k)
         else
            tkh(k) = tkh(k)*(mass(k) - massc(k))
         endif
 13   continue
*
* Dynamical timescale for the primary.
*
      tdyn = 5.05d-05*SQRT(rad(j1)**3/mass(j1))

*
* Set default qcrit values and identify special cases.
*

      if(qcflag.eq.0)then
*
* Use the BSE equations and defaults:
* If q1 = m_donor/m_acc > qc then common envelope
*

         if(kstar(j1).eq.0)then
            qc = 0.695
         elseif(kstar(j1).eq.1)then
            qc = 3.d0
         elseif(kstar(j1).eq.2)then
            qc = 4.d0
         elseif(kstar(j1).eq.3)then
            qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0
         elseif(kstar(j1).eq.4)then
            qc = 3.d0
         elseif(kstar(j1).eq.5)then
            qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0
         elseif(kstar(j1).eq.6)then
            qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0
         elseif(kstar(j1).eq.7)then
            qc = 3.0d0
         elseif(kstar(j1).eq.8)then
            qc = 0.784d0
         elseif(kstar(j1).eq.9)then
            qc = 0.784d0
         elseif(kstar(j1).ge.10)then
            qc = 0.628
         endif
      elseif(qcflag.eq.1)then
*
* Use the BSE equations and Hjellming & Webbink, 1987, ApJ, 318, 794.
* for the GB/AGB:
* If q1 = m_donor/m_acc > qc then common envelope
*
         if(kstar(j1).eq.0)then
            qc = 0.695
         elseif(kstar(j1).eq.1)then
            qc = 3.d0
         elseif(kstar(j1).eq.2)then
            qc = 4.d0
         elseif(kstar(j1).eq.3)then
            qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.4)then
            qc = 3.d0
         elseif(kstar(j1).eq.5)then
            qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.6)then
            qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.7)then
            qc = 3.0d0
         elseif(kstar(j1).eq.8)then
            qc = 0.784d0
         elseif(kstar(j1).eq.9)then
            qc = 0.784d0
         elseif(kstar(j1).ge.10)then
            qc = 0.628
         endif

      elseif(qcflag.eq.2)then
*
* Use the binary_c prescriptions taken from Claeys+2014 Table 2
* If q1 = m_donor/m_acc > qc then common envelope
*
         if(kstar(j2).lt.10)then
             if(kstar(j1).eq.0)then
                qc = 0.695
             elseif(kstar(j1).eq.1)then
                qc = 1.6d0
             elseif(kstar(j1).eq.2)then
                qc = 4.d0
             elseif(kstar(j1).eq.3)then
               qc = (1.67d0-zpars(7)+
     &               2.d0*(massc(j1)/mass(j1))**5)/2.13d0
             elseif(kstar(j1).eq.4)then
                qc = 3.d0
             elseif(kstar(j1).eq.5)then
               qc = (1.67d0-zpars(7)+
     &               2.d0*(massc(j1)/mass(j1))**5)/2.13d0
             elseif(kstar(j1).eq.6)then
               qc = (1.67d0-zpars(7)+
     &               2.d0*(massc(j1)/mass(j1))**5)/2.13d0
             elseif(kstar(j1).eq.7)then
                qc = 3.0d0
             elseif(kstar(j1).eq.8)then
                qc = 4.0d0
             elseif(kstar(j1).eq.9)then
                qc = 0.784d0
             elseif(kstar(j1).ge.10)then
                qc = 3.0d0
             endif
         elseif(kstar(j2).ge.10)then
             if(kstar(j1).eq.0)then
                qc = 1.0d0
             elseif(kstar(j1).eq.1)then
                qc = 1.0d0
             elseif(kstar(j1).eq.2)then
                qc = 4.7619d0
             elseif(kstar(j1).eq.3)then
                qc = 1.15d0
             elseif(kstar(j1).eq.4)then
                qc = 3.d0
             elseif(kstar(j1).eq.5)then
                qc = 1.15d0
             elseif(kstar(j1).eq.6)then
                qc = 1.15d0
             elseif(kstar(j1).eq.7)then
                qc = 3.0d0
             elseif(kstar(j1).eq.8)then
               qc = 4.7619d0
             elseif(kstar(j1).eq.9)then
                qc = 1.15d0
             elseif(kstar(j1).ge.10)then
                qc = 0.625d0
             endif
         endif
      elseif(qcflag.eq.3)then
*
* Use the binary_c prescriptions taken from Claeys+2014 Table 2
* but w/ Hjellming & Webbink for GB/AGB
* If q1 = m_donor/m_acc > qc then common envelope
*
         if(kstar(j2).lt.10)then
             if(kstar(j1).eq.0)then
                qc = 0.695
             elseif(kstar(j1).eq.1)then
                qc = 1.6d0
             elseif(kstar(j1).eq.2)then
                qc = 4.d0
             elseif(kstar(j1).eq.3)then
               qc = 0.362 +
     &              1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
             elseif(kstar(j1).eq.4)then
                qc = 3.d0
             elseif(kstar(j1).eq.5)then
               qc = 0.362 +
     &              1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
             elseif(kstar(j1).eq.6)then
               qc = 0.362 +
     &              1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
             elseif(kstar(j1).eq.7)then
                qc = 3.0d0
             elseif(kstar(j1).eq.8)then
                qc = 4.0d0
             elseif(kstar(j1).eq.9)then
                qc = 0.784d0
             elseif(kstar(j1).ge.10)then
                qc = 3.0d0
             endif
         elseif(kstar(j2).ge.10)then
             if(kstar(j1).eq.0)then
                qc = 1.0d0
             elseif(kstar(j1).eq.1)then
                qc = 1.0d0
             elseif(kstar(j1).eq.2)then
                qc = 4.7619d0
             elseif(kstar(j1).eq.3)then
                qc = 1.15d0
             elseif(kstar(j1).eq.4)then
                qc = 3.d0
             elseif(kstar(j1).eq.5)then
                qc = 1.15d0
             elseif(kstar(j1).eq.6)then
                qc = 1.15d0
             elseif(kstar(j1).eq.7)then
                qc = 3.0d0
             elseif(kstar(j1).eq.8)then
               qc = 4.7619d0
             elseif(kstar(j1).eq.9)then
                qc = 1.15d0
             elseif(kstar(j1).ge.10)then
                qc = 0.625d0
             endif
         endif
      elseif(qcflag.eq.4)then
*
* Use the StarTrack prescriptions taken from Belczynski+2008
* section 5.1
* If q1 = m_donor/m_acc > qc then common envelope
*
* Note: this does not do the WD mass ratios exactly the same.
* We do the standard 0.628 mass ratio from BSE since we don't
* compute the mass transfer rates the same
         if(kstar(j1).eq.0)then
            qc = 3.d0
         elseif(kstar(j1).eq.1)then
            qc = 3.d0
         elseif(kstar(j1).eq.2)then
            qc = 3.d0
         elseif(kstar(j1).eq.3)then
            qc = 3.d0
         elseif(kstar(j1).eq.4)then
            qc = 3.d0
         elseif(kstar(j1).eq.5)then
            qc = 3.d0
         elseif(kstar(j1).eq.6)then
            qc = 3.d0
         elseif(kstar(j1).eq.7)then
            qc = 1.7d0
         elseif(kstar(j1).eq.8)then
            qc = 3.5
         elseif(kstar(j1).eq.9)then
            qc = 3.5
         elseif(kstar(j1).ge.10)then
            qc = 0.628
         endif
      elseif(qcflag.eq.5)then
*
* Use the COMPAS prescriptions taken from Neijssel+2020,
* section 2.3
* We convert from radial response to qcrit for MS and HG,
* which assumes conservative mass transfer
* Stable MT is always assumed for stripped stars
* Assume standard qcrit from BSE for kstar>=10
* If q1 = m_donor/m_acc > qc then common envelope
*
         if(kstar(j1).eq.0)then
            qc = 1.717d0
         elseif(kstar(j1).eq.1)then
            qc = 1.717d0
         elseif(kstar(j1).eq.2)then
            qc = 3.825d0
         elseif(kstar(j1).eq.3)then
           qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.4)then
            qc = 3.d0
         elseif(kstar(j1).eq.5)then
           qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.6)then
           qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
         elseif(kstar(j1).eq.7)then
            qc = 1000.d0
         elseif(kstar(j1).eq.8)then
            qc = 1000.d0
         elseif(kstar(j1).eq.9)then
            qc = 1000.d0
         elseif(kstar(j1).ge.10)then
            qc = 0.628
         endif
      endif
*
* Allow for manually overriding qcrit values with fixed
* values supplied from ini file.
*
      qc_fixed = qcrit_array(kstar(j1)+1)
      if(qc_fixed.ne.0)then
         qc = qc_fixed
      endif

      if((kstar(j1).le.1.or.kstar(j1).eq.7).and.q(j1).gt.qc)then
*
* KB 20 jul 23: This will lead to a merger of the two MS stars.
* The result is always a single star placed in *2.
*
         CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(kstar(j2).le.1)then
*
* Assume perfect accretion efficiency
*
            dm2 = dm1
            mass(j2) = mass(j2) + dm2
*
* Rejuvenate if the star is still on the main sequence.
*
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,
     &                tscls,lums,GB,zpars)
* If the star has no convective core then the effective age decreases,
* otherwise it will become younger still.
            if(mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0)then
               aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm2)/mass(j2)
            else
               aj(j2) = tmsnew/tms(j2)*aj(j2)
            endif
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.6)then
*
* Add all the material to the giant's envelope.
*
            dm2 = dm1
            mass(j2) = mass(j2) + dm2
            if(kstar(j2).eq.2)then
               mass0(j2) = mass(j2)
               CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                   lums,GB,zpars)
               aj(j2) = tmsnew + tscls(1)*(aj(j2)-tms(j2))/tbgb(j2)
               epoch(j2) = tphys - aj(j2)
            endif
         elseif(kstar(j2).le.12)then
*
* Form a new giant envelope.
*
            dm2 = dm1
            kst = ktype(kstar(j1),kstar(j2))
            if(kst.gt.100) kst = kst - 100
            if(kst.eq.4)then
               aj(j2) = aj(j2)/tms(j2)
               massc(j2) = mass(j2)
            endif
*
* Check for planets or low-mass WDs.
*
            if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &         (kstar(j2).ge.11.and.mass(j2).lt.0.1d0))then
               kst = kstar(j1)
               mass(j1) = mass(j2) + dm2
               mass(j2) = 0.d0
            else
               mass(j2) = mass(j2) + dm2
               CALL gntage(massc(j2),mass(j2),kst,zpars,
     &                     mass0(j2),aj(j2))
               epoch(j2) = tphys - aj(j2)
            endif
            kstar(j2) = kst
         else
*
* The neutron star or black hole simply accretes at the Eddington rate.
*
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true.
            mass(j2) = mass(j2) + dm2
         endif
         coel = .true.
         binstate = 1
         if(mass(j2).gt.0.d0)then
            mass(j1) = 0.d0
            kstar(j1) = 15
         else
            kstar(j1) = kstar(j2)
            kstar(j2) = 15
         endif
         goto 135
*
*KB added RRLO_1 flag to send RRLO_1 > 10 into CE for post-MS binaries
*
      elseif(((kstar(j1).eq.3.or.kstar(j1).eq.5.or.kstar(j1).eq.6.or.
     &        kstar(j1).eq.8.or.kstar(j1).eq.9)
     &        .and.(q(j1).gt.qc.or.radx(j1).le.radc(j1))).or.
     &       (kstar(j1).eq.2.and.q(j1).gt.qc).or.
     &       (kstar(j1).eq.4.and.q(j1).gt.qc))then
*
* Common-envelope evolution.
*
         m1ce = mass(j1)
         m2ce = mass(j2)
         kcomp1 = kstar(j1) !PDK
         kcomp2 = kstar(j2)
         if(j1.eq.2)then
             switchedCE = .true.
         else
             switchedCE = .false.
         endif
         evolve_type = 7.d0
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),
     &                 kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')

         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel,j1,j2,
     &               vk,kick_info,formation(j1),formation(j2),sigmahold,
     &               bhspin(j1),bhspin(j2),binstate,mergertype,
     &               jp,tphys,switchedCE,rad,tms,evolve_type,disrupt,
     &               lumin,B_0,bacc,tacc,epoch,menv,renv,bkick,
     &               deltam1_bcm,deltam2_bcm)
         if(j1.eq.2.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif
         if(j1.eq.1.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif
         com = .true.
         if(com.and..not.coel.and..not.disrupt)then
* if it went through common envelope
* did not disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system is still in binary
            binstate = 0
            mergertype = -1
         elseif(com.and..not.coel.and.disrupt)then
* if it went through common envelope
* and did disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system should be marked as disrupted
            binstate = 2
            mergertype = -1
         endif
         if(binstate.eq.1.d0)then
             sep = 0.d0
             tb = 0.d0
         elseif(binstate.eq.2.d0)then
             sep = -1.d0
             tb = -1.d0
         endif
         if(j1.eq.2.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif
         if(j1.eq.1.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif

*
         evolve_type = 8.0
         
         
         mc = massc(1)
         rc = radc(1)
         CALL star(kstar(1),mass0(1),mass(1),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(1),aj(1),mass(1),tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kstar(1),mc,rc,me,re,k2,bhspin(1),1)
     
         rad(1) = rm
         lumin(1) = lum  
         massc(1) = mc
         radc(1) = rc
         menv(1) = me
         renv(1) = re
         
         
         mc = massc(2)
         rc = radc(2)
         CALL star(kstar(2),mass0(2),mass(2),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(2),aj(2),mass(2),tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kstar(2),mc,rc,me,re,k2,bhspin(2),2)
     
         rad(2) = rm
         lumin(2) = lum  
         massc(2) = mc
         radc(2) = rc
         menv(2) = me
         renv(2) = re

         mass1_bpp = mass(1)
         mass2_bpp = mass(2)
         if(kstar(1).eq.15) mass1_bpp = mass0(1)
         if(kstar(2).eq.15) mass2_bpp = mass0(2)
         rrl1 = rad(1)/rol(1)
         rrl2 = rad(2)/rol(2)
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass1_bpp,mass2_bpp,
     &                 kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
*
         epoch(j1) = tphys - aj(j1)
         com = .false.
         if(coel)then
            com = .true.
            goto 135
         endif
         epoch(j2) = tphys - aj(j2)
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
*
* Next step should be made without changing the time.
*
         dm1 = m1ce - mass(j1)
         dm2 = mass(j2) - m2ce
         dm22 = dm2
         dtm = 0.d0
*
* Reset orbital parameters as separation may have changed.
*
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      elseif(kstar(j1).ge.10.and.kstar(j1).le.12.and.
     &       q(j1).gt.qc)then
*
* Dynamic transfer from a white dwarf.  Secondary will have KW > 9.
*
         CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(eddfac.lt.10.d0)then
*           dm2 = MIN(dme*taum/tb,dm1)
            if(wd_mass_lim.eq.1)then
               dm2 = MIN(dme*taum/tb,dm1)
            elseif(wd_mass_lim.eq.0)then
               dm2 = dm1
            endif
            if(dm2.lt.dm1) supedd = .true.
         else
            dm2 = dm1
         endif
         mass(j2) = mass(j2) + dm2
*
         if(kstar(j1).eq.10.and.kstar(j2).eq.10)then
*
* Assume the energy released by ignition of the triple-alpha reaction
* is enough to destroy the star.
*
            kstar(j2) = 15
            mass(j2) = 0.d0
         elseif(kstar(j1).eq.10.and.kstar(j2).gt.10)then
*
* Should be helium overflowing onto a CO or ONe core in which case the
* helium swells up to form a giant envelope so a HeGB star is formed.
* Allowance for the rare case of CO or ONe flowing onto He is made.
*
            kst = 9
            if(kstar(j2).eq.10) massc(j2) = dm2
            CALL gntage(massc(j2),mass(j2),kst,zpars,mass0(j2),aj(j2))
            kstar(j2) = kst
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.12)then
            mass0(j2) = mass(j2)
            if(kstar(j1).eq.12.and.kstar(j2).eq.11)then
*
* Mixture of ONe and CO will result in an ONe product.
*
               kstar(j2) = 12
            endif
            if(kstar(j1).eq.11.and.kstar(j2).eq.11)then !PK.
*
* Following Startrack of merger of two CO will result in an ONe or NS product
* (depending upon mass, if kstar = 12 and mass > Mecsn then hrdiag will update
* accordingly).
*
               kstar(j2) = 12
            endif
            formation(j2) = 11
         endif
         kstar(j1) = 15
         mass(j1) = 0.d0
*
* Might be a supernova that destroys the system.
*
         if(kstar(j2).le.11.and.mass(j2).gt.mch)then
            kstar(j2) = 15
            mass(j2) = 0.d0
         endif
         coel = .true.
         binstate = 1
         goto 135
      elseif(kstar(j1).eq.13)then
*
* Gamma ray burster?
*
         CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         kstar(j2) = 14
         coel = .true.
         binstate = 1
         goto 135
      elseif(kstar(j1).eq.14)then
*
* Both stars are black holes.  Let them merge quietly.
*
         CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         coel = .true.
         binstate = 1
         goto 135
      else
*
* Mass transfer in one Kepler orbit.
*
*
* KB: adding in stable mass transfer factor from
*     eqs 10-11 of Claeys+2014 as don_lim flag instead
*     of qcflag (10/12/20)
*
         if(don_lim.eq.-2)then
            if(q(j1).gt.1)then
               f_fac=1000.d0
            else
               f_fac=1000*q(j1)*
     &               EXP(-0.5d0*(-LOG(q(j1))/0.15d0)**2)
               if(f_fac.lt.1)then
                  f_fac=1
               endif
            endif
            dm1 = f_fac*3.0d-06*tb*(LOG(rad(j1)/rol(j1))**3)*
     &            MIN(mass(j1),5.d0)**2
         elseif(don_lim.eq.-1)then
            dm1 = 3.0d-06*tb*(LOG(rad(j1)/rol(j1))**3)*
     &            MIN(mass(j1),5.d0)**2
         endif
*         dm1 = MIN(dm1, mass(j1)*tb/tkh(j1))
         if(kstar(j1).eq.2)then
            mew = (mass(j1) - massc(j1))/mass(j1)
            dm1 = MAX(mew,0.01d0)*dm1
         elseif(kstar(j1).ge.10)then
*           dm1 = dm1*1.0d+03/MAX(rad(j1),1.0d-04)
            dm1 = dm1*1.0d+03*mass(j1)/MAX(rad(j1),1.0d-04)
         endif
         kst = kstar(j2)
*
* Possibly mass transfer needs to be reduced if primary is rotating
* faster than the orbit (not currently implemented).
*
*        spnfac = MIN(3.d0,MAX(ospin(j1)/oorb,1.d0))
*        dm1 = dm1/spnfac**2
*
* Limit mass transfer to the thermal rate for all fusing stars
* and to the dynamical rate for all others.
*
         if(kstar(j1).ge.0.and.kstar(j1).le.9)then
***
* JH_temp ... this may be good for HG RLOF??
*           if(kstar(j1).eq.2)then
*              mew = rad(j1)/rol(j1) - 1.d0
*              mew = 2.d0*mew
*              dm1 = dm1*10.d0**mew
*           endif
***
            dm1 = MIN(dm1,mass(j1)*tb/tkh(j1))
         endif
*
* Now we will initiate a CE for giants or merger for MS
* if RRLO_1 > 10. Previously a MS + MS & q > c conditional
* existed here, but it is completely defunct at this point
*

         if(rad(j1).gt.10.d0*rol(j1).and.kstar(j1).le.1.and.
     &      kstar(j2).le.1)then
*
* Allow the stars to merge with the product in *1.
*
            CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
            m1ce = mass(j1)
            m2ce = mass(j2)
            if((kstar(1).ge.10.and.kstar(1).le.12).and.
     &         (kstar(2).ge.10.and.kstar(2).le.12))then
               formation(1) = 11
               formation(2) = 11
            endif
            CALL mix(mass0,mass,aj,kstar,zpars,bhspin)
            dm1 = m1ce - mass(j1)
            dm2 = mass(j2) - m2ce
*
* Next step should be made without changing the time.
*
            dtm = 0.d0
            epoch(1) = tphys - aj(1)
            coel = .true.
            binstate = 1
            goto 135
         endif
         if(kstar(j1).gt.9)then
            dm1 = MIN(dm1,mass(j1)*tb/tdyn)
         endif
*
* Calculate wind mass loss from the stars during one orbit.
*
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 14 , k = 1,2
*
* Determine the eddington limit for the accretor (3-k)
* Just in case the wind mass loss rates are *very* high
*
            dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(3-k)*tb
* For BHs, follow Marchant et al. 2017
            if(kstar(3-k).eq.14)then
               maxspinBH = 6.d0**(1.d0/2.d0) * Mbh_initial
               if(mass(3-k).lt.maxspinBH)then
                  etaBH = 1.d0 -
     &    (1.d0 - (mass(3-k)/(3.d0*Mbh_initial))**(2.d0))**(1.d0/2.d0)
               else
                  etaBH = 0.42
               endif
               dme = 1.04e-3*eddfac*(1.d0/(1.d0 + zpars(11)))
     &    *(1.d0/etaBH)*rad(3-k)*tb
            endif

            !check for kstar added by PA
            if(neta.gt.tiny .and. kstar(k)<15)then
               if(beta.lt.0.d0)then !PK. following startrack
                  beta = 0.125
                  if(kstar(k).le.1)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.1.4d0)then
                        beta = 0.5
                     else
                        beta = 7.d0*((mass(k)-1.4d0)/(120.d0-1.4d0))
     &                         + 0.5d0
                     endif
                  elseif(kstar(k).ge.7.and.kstar(k).le.9)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.10.d0)then
                        beta = 0.125
                     else
                        beta = 7.d0*((mass(k)-10.d0)/(120.d0-10.d0))
     &                               + 0.125d0
                     endif
                  endif
               endif
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),radx(k),
     &                         mass(k),massc(k),rlperi,z)
               vwind2 = 2.d0*beta*acc1*mass(k)/radx(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),dmr(k))
*
* Apply Eddington limit
*
               dmt(3-k) = MIN(dmt(3-k),dme/tb)
               if(dmt(3-k).eq.dme/tb) supedd = .true.
               beta = betahold
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 14      continue
*
         do 15 , k = 1,2
            dms(k) = (dmr(k)-dmt(k))*tb
 15      continue
*
* Increase time-scale to relative mass loss of 0.5% but not more than twice.
* KM is the number of orbits for the timestep.
*
         km = MIN(2.d0*km0,5.0d-03/
     &            MAX(ABS(dm1+dms(j1))/mass(j1),dms(j2)/mass(j2)))
         km0 = km
*
*       Modify time-step & mass loss terms by speed-up factor.
*
         dt = km*tb
         dtm = dt/1.0d+06
*
* Take the stellar evolution timestep into account but don't let it
* be overly restrictive for long lived phases.
*

* NOTE: This can cause NaN values! If you get NaNs, check what happens
*       when you increase `loop` (see PR #647) - TW
*
         if(iter.le.loop) dtm = MIN(dtm,dtmi(1),dtmi(2))

         dtm = MIN(dtm,tsave-tphys)
         dt = dtm*1.0d+06
         km = dt/tb
*
* Decide between accreted mass by secondary and/or system mass loss.
*
         taum = mass(j2)/dm1*tb


*
* KB 4/Jan/21: adding in acc_lim flags
* MZ 30/Mar/21: slight adjustment to flag designations
* acc_lim >= 0: fraction of donor mass loss accreted; supersceded by
*              eddington limits and/or novae on WDs.
* acc_lim = -1: standard BSE w/ MS/HG/CHeB assumed to have thermal
*              limit of 10*tkh while giants are unlimited. If the
*              accretor is He-rich (>=7), limit is 10*tkh if accretor
*              is He-MS/HeHG/HeAGB and unlimited if donor is H-rich
* acc_lim = -2: same as acc_lim = -1, but for 1*tkh instead
* acc_lim = -3: accretor kw=0-6 have thermal limit of 10*tkh. If the
*               acretor is He-rich (>=7), limit is 10*tkh if donor
*               is He-MS/HeHG/HeAGB and unlimited if donor is H-rich
* acc_lim = -4: accretor kw=0-6 have thermal limit of tkh. If the
*               accretor is He-rich (>=7), limit is 1*tkh if donor
*               is He-MS/HeHG/HeAGB and unlimited if donor is H-rich
*
* Note that acc_lim >= 0 means that the thermal limit is not considered
*


         if(kstar(j2).le.2.or.kstar(j2).eq.4)then
            if(acc_lim.eq.-1.or.acc_lim.eq.-3)then
               dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
            elseif(acc_lim.eq.-2.or.acc_lim.eq.-4)then
               dm2 = MIN(1.d0,taum/tkh(j2))*dm1
            elseif(acc_lim.ge.0.d0)then
               dm2 = acc_lim*dm1
            endif
         elseif(kstar(j2).ge.7.and.kstar(j2).le.9)then
*
* Naked helium star secondary swells up to a core helium burning star
* or SAGB star unless the primary is also a helium star.
*
            if(kstar(j1).ge.7)then
               if(acc_lim.eq.-1.or.acc_lim.eq.-3)then
                  dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
               elseif(acc_lim.eq.-2.or.acc_lim.eq.-4)then
                  dm2 = MIN(1.d0,taum/tkh(j2))*dm1
               elseif(acc_lim.ge.0.d0)then
                  dm2 = acc_lim*dm1
               endif
            else
               if(acc_lim.lt.0.d0)then
                  dm2 = dm1
               elseif(acc_lim.ge.0.d0)then
                  dm2 = acc_lim*dm1
               endif
               dmchk = dm2 - 1.05d0*dms(j2)
               if(dmchk.gt.0.d0.and.dm2/mass(j2).gt.1.0d-04)then
                  kst = MIN(6,2*kstar(j2)-10)
                  if(kst.eq.4)then
                     aj(j2) = aj(j2)/tms(j2)
                     mcx = mass(j2)
                  else
                     mcx = massc(j2)
                  endif
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(mcx,mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
               endif
            endif
         elseif(kstar(j1).le.6.and.
     &           (kstar(j2).ge.10.and.kstar(j2).le.12))then
*
* White dwarf secondary.
*
            if(dm1/tb.lt.2.71d-07)then
               if(dm1/tb.lt.1.03d-07)then
*
* Accrete until a nova explosion blows away most of the accreted material.
*
                  novae = .true.
                  if(acc_lim.lt.0.d0)then
                     dm2 = MIN(dm1,dme)
                     if(dm2.lt.dm1) supedd = .true.
                  elseif(acc_lim.ge.0.d0)then
                     dm2 = MIN(dm2,acc_lim*dm1)
                     if(dm2.lt.acc_lim*dm1) supedd = .true.
                  endif
                  dm22 = epsnov*dm2
               else
*
* Steady burning at the surface
*
                  if(acc_lim.lt.0.d0)then
                     dm2 = dm1
                  elseif(acc_lim.ge.0.d0)then
                     dm2 = acc_lim*dm1
                  endif
               endif
            else
*
* Make a new giant envelope.
*
               if(acc_lim.lt.0.d0)then
                  dm2 = dm1
               elseif(acc_lim.ge.0.d0)then
                  dm2 = MIN(dm2,acc_lim*dm1)
               endif
*
* Check for planets or low-mass WDs.
*
               if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &            (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
                  kst = kstar(j2)
               else
                  kst = MIN(6,3*kstar(j2)-27)
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(massc(j2),mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
*
               endif
*
            endif
         elseif(kstar(j2).eq.3.or.kstar(j2).eq.5.or.kstar(j2).eq.6)then
* We have a giant w/ kstar(j2) = 3,5,6
*
            if(acc_lim.eq.-1.or.acc_lim.eq.-2)then
               dm2 = dm1
            elseif(acc_lim.eq.-3)then
               dm2 = MIN(1.d0,10*taum/tkh(j2))*dm1
            elseif(acc_lim.eq.-4)then
               dm2 = MIN(1.d0,taum/tkh(j2))*dm1
            elseif(acc_lim.ge.0.d0)then
               dm2 = MIN(dm2,acc_lim*dm1)
            endif

         endif

*
* Impose the Eddington limit.
*
         if(kstar(j2).ge.10)then
            if(acc_lim.lt.0.d0)then
*
* If there is wind accretion the total amount of mass change is
* dms(j2) = dmr(j2) - dmt(j2), where dmt(j2) is the accretion
* from the companion. We should limit to the Eddington limit minus
* the amount of accretion that is already coming in from Winds
*
               dm2 = MIN(dm1,dme + dms(j2))
*
* If we already hit supereddington wind accretion, don't add
* any more mass through RLO
*
               if(supedd.eqv..true.) dm2 = 0.d0
               if(dm2.lt.dm1) supedd = .true.
            elseif(acc_lim.ge.0.d0)then
*
* If there is wind accretion the total amount of mass change is
* dms(j2) = dmr(j2) - dmt(j2), where dmt(j2) is the accretion
* from the companion. We should limit to the Eddington limit minus
* the amount of accretion that is already coming in from Winds
*
               dm2 = MIN(acc_lim*dm1,dme + dms(j2))
*
* If we already hit supereddington wind accretion, don't add
* any more mass through RLO
*
               if(supedd.eqv..true.) dm2 = 0.d0
               if(dm2.lt.acc_lim*dm1) supedd = .true.
            endif

*
* Can add pulsar propeller evolution here if need be. PK.
*
         endif


         if(.not.novae) dm22 = dm2
*
         if(kst.ge.10.and.kst.le.12)then
            mt2 = mass(j2) + km*(dm22 - dms(j2))
            if(kstar(j1).le.10.and.kst.eq.10.and.mt2.ge.0.7d0)then
*
* HeWD can only accrete helium-rich material up to a mass of 0.7 when
* it is destroyed in a possible Type 1a SN.
*
               mass(j1) = mass(j1) - km*(dm1 + dms(j1))
               mass(j2) = 0.d0
               kstar(j2) = 15
               goto 135
            elseif(kstar(j1).le.10.and.kst.ge.11)then
*
* CO and ONeWDs accrete helium-rich material until the accumulated
* material exceeds a mass of 0.15 when it ignites. For a COWD with
* mass less than 0.95 the system will be destroyed as an ELD in a
* possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs
* will survive with all the material converted to ONe (JH 30/09/99).
*
** Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).
*
               if((mt2-mass0(j2)).ge.0.15d0)then
                  if(kst.eq.11)then
                     mass(j1) = mass(j1) - km*(dm1 + dms(j1))
                     mass(j2) = 0.d0
                     kstar(j2) = 15
                     goto 135
                  endif
                  mass0(j2) = mt2
               endif
            else
               mass0(j2) = mt2
            endif
*
* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
* the white dwarf in a supernova. If the WD is ONe then a neutron star
* will survive the supernova and we let HRDIAG take care of this when
* the stars are next updated.
*
            if(kst.eq.10.or.kst.eq.11)then
               if(mt2.ge.mch)then
                  dm1 = mch - mass(j2) + km*dms(j2)
                  mass(j1) = mass(j1) - dm1 - km*dms(j1)
                  mass(j2) = 0.d0
                  kstar(j2) = 15
                  goto 135
               endif
            endif
         endif
*
*       Modify mass loss terms by speed-up factor.
*
         dm1 = km*dm1
         dm2 = km*dm2
         dm22 = km*dm22
         dme = km*dme

*
* Calculate orbital angular momentum change due to system mass loss.
*
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) +
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))/
     &           (mass(1)+mass(2))**2
         djorb = djorb*dt
*
* For super-Eddington mass transfer rates, for gamma == -2,
* and for novae systems, assume that material is lost from
* the system as if a wind from the secondary.
*
         if(supedd.or.novae.or.int(gamma).eq.-2)then
            djorb = djorb + (dm1 - dm22)*mass(j1)*mass(j1)/
     &              (mass(1)+mass(2))**2
* If gamma == -3: Assume mass is lost through the outer 
* Lagrangian point, forming a circumbinary disk. See
* Zapartas+17 Eq. 9 and Artymowicz & Lubow (1994).
* Set rmin=2a, which is consistent with A&L94 simulations.
         elseif(int(gamma).eq.-3)then
            gammadisc = (mass(1) + mass(2))**2
     &                / (mass(1) * mass(2)) * sqrt(2.d0)
            djorb = djorb + gammadisc * (dm1 - dm2)
* If gamma == -1 then assume the lost material carries with it
* the specific angular momentum of the primary and for all
         elseif(int(gamma).eq.-1) then
            djorb = djorb + (dm1 - dm2)*mass(j2)*mass(j2)/
     &              (mass(1)+mass(2))**2
* gamma > 0.0 assume that it takes away a fraction gamma of
* the orbital angular momentum.
         elseif(gamma.ge.0.d0)then
            djorb = djorb + gamma*(dm1 - dm2)
         endif
*
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = djorb*sep*sep*sqome2*oorb
         delet = 0.d0
*
* For very close systems include angular momentum loss mechanisms.
*
         if(sep.le.100000.d0.and.grflag.eq.1)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr*dt
            delet = delet + delet1*dt
         endif
*
         do 602 , k = 1,2
*
            dms(k) = km*dms(k)
*            WRITE(*,*)dme/tb,dms(j2)/tb/km,dmt(j2),dms(j1)/tb/km,dmr(j1)

            if(kstar(k).lt.10) dms(k) = MIN(dms(k),mass(k) - massc(k))
*
* Calculate change in the intrinsic spin of the star.
*
* Spin up via wind accretion. NSs are allowed to vary the angular momentum accretion efficiency. PK.
            if(kstar(k).eq.13.and.xip.eq.1.and.pulsar.gt.0)then
               if(B(k).gt.0.d0)then
                  xi = MIN(1.d0,(0.01d0*(2.0d+11/B(k)))+0.01d0)
                  xi = 0.01
               else
                  xi = 0.01
               endif
            else
               xi = xihold
            endif
            djspint(k) = (2.d0/3.d0)*(dmr(k)*radx(k)*radx(k)*ospin(k) -
     &                   xi*dmt(k)*radx(3-k)*radx(3-k)*ospin(3-k))
            djspint(k) = djspint(k)*dt
*
* Evaluate convective/radiative limits for a variety of stars as based
* on the work of Belczynski et al. (2008). As option.
* Note only certain stars of type k = 0, 1, 2, 4, 5, 6, 9 **double check
* 0, 1 have range of 0.35-Mms,conv. Mms,conv is function of metalicity.
* 2 & 4 have a temperature dependence. Convective if Teff < 10**3.73
* 3, 5 & 6 are giants and have no radiative envelopes.
* 9's envelope is convective when M < Mhe,conv, Mhe,conv = 3.0Msun.
*
            if(ST_cr.le.0)then
               if(kstar(k).le.1)then
                  convradcomp(k) = 1.25d0
               else
                  convradcomp(k) = 9999999999.d0 !doesn't work well with startrack_tide = .false.
               endif
            else
               if(kstar(k).le.1)then
*                 Main sequence mass limit varying with metallicity
                  convradcomp(k) = 1.25d0
                  if(z.gt.0.001d0.and.z.lt.zsun)then
                     convradcomp(k) = 0.747d0 + 55.73d0*z - 1532*z*z
                  elseif(z.le.0.001d0)then
                     convradcomp(k) = 0.8d0
                  endif
               elseif(kstar(k).eq.2.or.kstar(k).eq.4)then
*                 H-rich HG and CHeB temperature limit
                  convradcomp(k) = 10.d0**3.73d0
               elseif(kstar(k).le.6)then
*                 No limit or all other giant stars
*                 (compare this large value to mass).
                  convradcomp(k) = 9999999999.d0
               elseif(kstar(k).eq.9)then
*                 Mass limit for evolved He star
                  convradcomp(k) = 3.d0
               endif
            endif
*
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
*               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
*               djspint(k) = djspint(k) + djmb*dt
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
*     &              menv(k).gt.0.0d0)then
*                if (ospin(k) .le. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &              (ospin(k)/wsun)**3.0d0
*                if (ospin(k) .gt. wx) djmb = kw3 * rad(k)**4.0d0 *
*     &             (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
*                djspint(k) = djspint(k) + djmb
*            endif
* MB
            djmb = 0.d0
            if(htpmb.eq.0)then
* HTP02 method
               if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
                  djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
                  djspint(k) = djspint(k) + djmb*dt
               endif
            elseif(htpmb.eq.1)then
               if(ST_cr.le.0.and.
     &            (mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
     &            menv(k).gt.0.0d0))then
* Ivanova & Taam (2003)
                   if (ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                 (ospin(k)/wsun)**3.0d0
                   if (ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                   djspint(k) = djspint(k) + djmb
               elseif(ST_cr.gt.0.and.(menv(k).gt.0.d0.and.
     &            ((kstar(k).le.1.and.mass(k).gt.0.35d0.and.
     &              mass(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.2.and.teff(k).le.convradcomp(k)).or.
     &              (kstar(k).eq.4.and.teff(k).le.convradcomp(k)).or.
     &              ((kstar(k).eq.3).or.(kstar(k).eq.5).or.
     &              (kstar(k).eq.6)))))then
* MB given in Ivanova & Taam (2003)
*            if(mass(k).gt.0.35d0.and.kstar(k).lt.10.and.
*     &              menv(k).gt.0.0d0)then
                  if(ospin(k).le.wx) djmb = kw3 * rad(k)**4.0d0 *
     &                (ospin(k)/wsun)**3.0d0
                  if(ospin(k).gt.wx) djmb = kw3 * rad(k)**4.0d0 *
     &               (ospin(k)/wsun)**1.3d0 * (wx/wsun)**1.7d0
                  djspint(k) = djspint(k) + djmb
               endif
            endif
            if(kstar(k).eq.13.and.pulsar.gt.0)then
               if(bdecayfac.eq.0)then
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).le.tiny)then
                     B(k) = B_0(k)*exp(-CK*bacc(k)) + Bbot
                  else
                     B(k) = B_0(k)*
     &                      exp(-(tphys-epoch(k)-tacc(k))/bconst)*
     &                      exp(-CK*bacc(k)) + Bbot
                  endif
               else
                  if(B_0(k).eq.0.d0)then
                     B(k) = 0.d0
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny.and.
     &                bacc(k).eq.0.d0)then
                     B(k) = B_0(k) + Bbot
                  elseif((tphys-epoch(k)-tacc(k)).lt.tiny)then
                     B(k) = B_0(k)/(1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  elseif(bacc(k).eq.0.d0)then
                     B(k) = B_0(k)*EXP(-(tphys-epoch(k)-tacc(k))/bconst)
     &                       + Bbot
                  else
                     B(k) = B_0(k)*
     &                       EXP(-(tphys-epoch(k)-tacc(k))/bconst)/
     &                       (1.d0 + (bacc(k)/1.0d-6)) + Bbot
                  endif
               endif
            endif
*
 602     continue
*
* Adjust the spin angular momentum of each star owing to mass transfer
* and conserve total angular momentum.
*
         djt = dm1*radx(j1)*radx(j1)*ospin(j1)
         djspint(j1) = djspint(j1) + djt
         djorb = djorb - djt
         if(disk)then
*
* Alter spin of the degenerate secondary by assuming that material
* falls onto the star from the inner edge of a Keplerian accretion
* disk and that the system is in a steady state.
*
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*radx(j2))
            djspint(j2) = djspint(j2) - djt
            djorb = djorb + djt
*
         else
*
* No accretion disk.
* Calculate the angular momentum of the transferred material by
* using the radius of the disk (see Ulrich & Burger) that would
* have formed if allowed.
*
            rdisk = 1.7d0*rmin
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*rdisk)
            djspint(j2) = djspint(j2) - djt
            djorb = djorb + djt
*
         endif
         djtx(2) = djt
*
* Adjust the secondary spin if a nova eruption has occurred.
*
         if(novae)then
            djt = (dm2 - dm22)*radx(j2)*radx(j2)*ospin(j2)
            djspint(j2) = djspint(j2) + djt
            djtx(2) = djtx(2) - djt
         endif
*
* Calculate circularization, orbital shrinkage and spin up.
*
         do 603 , k = 1,2
*
* Evaluate convective/radiative limits for a variety of stars as based
* on the work of Belczynski et al. (2008). As option.
* Note only certain stars of type k = 0, 1, 2, 4, 5, 6, 9 **double check
* 0, 1 have range of 0.35-Mms,conv. Mms,conv is function of metalicity.
* 2 & 4 have a temperature dependence. Convective if Teff < 10**3.73
* 3, 5 & 6 are giants and have no radiative envelopes.
* 9's envelope is convective when M < Mhe,conv, Mhe,conv = 3.0Msun.
*
            if(ST_cr.le.0)then
               if(kstar(k).le.1)then
                  convradcomp(k) = 1.25d0
               else
                  convradcomp(k) = 99999999.d0
               endif
            else
               if(kstar(k).le.1)then
*                 Main sequence mass limit varying with metallicity
                  convradcomp(k) = 1.25d0
                  if(z.gt.0.001d0.and.z.lt.zsun)then
                     convradcomp(k) = 0.747d0 + 55.73d0*z - 1532*z*z
                  elseif(z.le.0.001d0)then
                     convradcomp(k) = 0.8d0
                  endif
               elseif(kstar(k).eq.2.or.kstar(k).eq.4)then
*                 H-rich HG and CHeB temperature limit
                  convradcomp(k) = 10.d0**3.73d0
               elseif(kstar(k).le.6)then
*                 No limit or all other giant stars
*                 (compare this large value to mass).
                  convradcomp(k) = 99999999.d0
               elseif(kstar(k).eq.9)then
*                 Mass limit for evolved He star
                  convradcomp(k) = 3.d0
               endif
            endif
*
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (radx(k)/sep)**2
               raa6 = raa2**3
*
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if(ST_cr.le.0.and.
     &            ((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7))then
*
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).
*
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*radx(k)*radx(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(ST_cr.gt.0.and.
     &                ((kstar(k).le.1.and.mass(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.2.and.teff(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.4.and.teff(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.9.and.mass(k).gt.convradcomp(k)).or.
     &           (kstar(k).eq.7).or.(kstar(k).eq.8)))then
*
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).
*
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*rad(k)*rad(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(kstar(k).le.9)then
* Convective damping

* In BSE paper Equation 30, the default scaling coefficient is 2./21
* the fprimc_array kstar dependent array that is fed in
* keeps this same coefficient by default but allows user to
* specify their own
*
                  renv(k) = MIN(renv(k),radx(k)-radc(k))
                  renv(k) = MAX(renv(k),1.0d-10)
                  tc = mr23yr*(menv(k)*renv(k)*(radx(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
                  f = MIN(1.d0,(ttid/(2.d0*tc))**2)
                  !kstar(k) to kstar(k)+1 by PA for kstar(k)=0
                  tcqr = fprimc_array(kstar(k)+1)*
     &                 f*q(3-k)*raa6*menv(k)/
     &                 (tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               elseif(ST_tide.le.0)then
*Degenerate damping
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1*dt
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
               eqspin = oorb*f2/(sqome3*f5)
               if(dt.gt.0.d0)then
                  if(dspint(k).ge.0.d0)then
                     dspint(k) = MIN(dt*dspint(k),eqspin-ospin(k))/dt
                  else
                     dspint(k) = MAX(dt*dspint(k),eqspin-ospin(k))/dt
                  endif
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt*dt
               djspint(k) = djspint(k) - djt*dt
*
            endif
*
            jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k))
*
* Ensure that the star does not spin up beyond break-up, and transfer
* the excess angular momentum back to the orbit.
*
            ospbru = twopi*SQRT(mass(k)*aursun**3/radx(k)**3)
            jspbru = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*ospbru
            if(jspin(k).gt.jspbru)then
               mew = 1.d0
               if(djtx(2).gt.0.d0)then
                  mew = MIN(mew,(jspin(k) - jspbru)/djtx(2))
               endif
               djorb = djorb - (jspin(k) - jspbru)
               jspin(k) = jspbru
* If excess material should not be accreted, activate next line.
*              dm22 = (1.d0 - mew)*dm22
            endif
*
 603     continue
*
* Update the masses.
*
         kstar(j2) = kst
         mass(j1) = mass(j1) - dm1 - dms(j1)
         if(kstar(j1).le.1.or.kstar(j1).eq.7) mass0(j1) = mass(j1)
         mass(j2) = mass(j2) + dm22 - dms(j2)
         if(kstar(j2).le.1.or.kstar(j2).eq.7) mass0(j2) = mass(j2)
         if(kstar(j2).eq.13.and.pulsar.gt.0)then
* if propeller add and if(.not.prop)then here... PK.
            b_mdot = (dm2-dms(j2))/dt
            if(b_mdot_lim.gt.0.d0.and.b_mdot.gt.b_mdot_lim)then
               bacc(j2) = bacc(j2) + (dm2-dms(j2))
               tacc(j2) = tacc(j2) + dtm
            elseif(b_mdot_lim.le.0.d0)then
               bacc(j2) = bacc(j2) + (dm2-dms(j2))
               tacc(j2) = tacc(j2) + dtm
            endif
         endif
*
* For a HG star check if the initial mass can be reduced.
*
         if(kstar(j1).eq.2.and.mass0(j1).le.zpars(3))then
            m0 = mass0(j1)
            mass0(j1) = mass(j1)
            CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j1))then
               mass0(j1) = m0
            endif
         endif
         if(kstar(j2).eq.2.and.mass0(j2).le.zpars(3))then
            m0 = mass0(j2)
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j2))then
               mass0(j2) = m0
            endif
         endif
*
         ecc = ecc - delet
         ecc = MAX(ecc,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
*
* Ensure that Jorb does not become negative which could happen if the
* primary overfills its Roche lobe initially. In this case we simply
* allow contact to occur.
*
* NOTE: there is a possibility that if contact occurs and we're at/past
* the current tphysf, then the binary will spiral in and we don't catch
* it.  If that happens, we could return an impossibly small binary which
* will break energy conservation in CMC.  Flag this here
         if((jorb-djorb).lt.1.d0) inspiral=.true.
         jorb = MAX(1.d0,jorb - djorb)
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
*
      endif
*
* Always rejuvenate the secondary and age the primary if they are on
* the main sequence.
*
      if(kstar(j1).le.2.or.kstar(j1).eq.7)then
         CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j1).eq.2)then
            aj(j1) = tmsnew + (tscls(1) - tmsnew)*(aj(j1)-tms(j1))/
     &                        (tbgb(j1) - tms(j1))
         else
            aj(j1) = tmsnew/tms(j1)*aj(j1)
         endif
         epoch(j1) = tphys - aj(j1)
      endif
*
      if(kstar(j2).le.2.or.kstar(j2).eq.7)then
         CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j2).eq.2)then
            aj(j2) = tmsnew + (tscls(1) - tmsnew)*(aj(j2)-tms(j2))/
     &                        (tbgb(j2) - tms(j2))
         elseif((mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0).
     &           and.kstar(j2).ne.7)then
            aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm22)/mass(j2)
         else
            aj(j2) = tmsnew/tms(j2)*aj(j2)
         endif
         epoch(j2) = tphys - aj(j2)
      endif
*
* Obtain the stellar parameters for the next step.
*
      tphys = tphys + dtm
      do 90 , k = 1,2
         age = tphys - epoch(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
*
* Masses over 100Msun should probably not be trusted in the
* evolution formulae.
*
         if(mt.gt.100.d0)then
*            WRITE(99,*)' MASS EXCEEDED ',mass1i,mass2i,tbi,ecci,mt,
*     & tphysf,id1_pass,id2_pass
*            goto 140
         endif
         kw = kstar(k)
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2,bhspin(k),k)
         pd = sep*(1.d0 - ecc)
         if(pd.lt.(rad(1)+rad(2))) goto 130


     
*
* Check for a supernova and correct the semi-major axis if so.
*
         if(kw.ne.kstar(k).and.kstar(k).le.12.and.
     &      (kw.eq.13.or.kw.eq.14))then
            dms(k) = mass(k) - mt
            if(formation(k).ne.11) formation(k) = 1
            if(kw.eq.13.and.ecsn.gt.0.d0)then
               if(kstar(k).le.6)then
                  if(mass0(k).le.zpars(5))then
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation(k) = 2
                  endif
               elseif(kstar(k).ge.7.and.kstar(k).le.9)then
                  if(mass(k).gt.ecsn_mlow.and.mass(k).le.ecsn)then
* BSE orgi: 1.6-2.25, Pod: 1.4-2.5, StarTrack: 1.83-2.25 (all in Msun)
                     if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                        sigma = -sigmahold/sigmadiv
                     elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                        sigma = sigmadiv
                     endif
                     formation(k) = 2
                  endif
               elseif(formation(k).eq.11)then
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 5
               elseif(kstar(k).ge.10.or.kstar(k).eq.12)then
* AIC formation, will never happen here but...
                  if(sigma.gt.0.d0.and.sigmadiv.gt.0.d0)then
                     sigma = -sigmahold/sigmadiv
                  elseif(sigma.gt.0.d0.and.sigmadiv.lt.0.d0)then
                     sigma = sigmadiv
                  endif
                  formation(k) = 4
               endif
            endif

            evolve_type = 14.d0 + FLOAT(k)
            teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
            teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
            if(B_0(1).eq.0.d0)then !PK.
               b01_bcm = 0.d0
            elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
               b01_bcm = B_0(1)
            else
               b01_bcm = B(1)
            endif
            if(B_0(2).eq.0.d0)then
               b02_bcm = 0.d0
            elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
               b02_bcm = B_0(2)
            else
               b02_bcm = B(2)
            endif
            CALL writetab(jp,tphys,evolve_type,
     &                    mass(1),mass(2),kstar(1),kstar(2),
     &                    sep,tb,ecc,rrl1,rrl2,
     &                    aj(1),aj(2),tms(1),tms(2),
     &                    massc(1),massc(2),rad(1),rad(2),
     &                    mass0(1),mass0(2),lumin(1),lumin(2),
     &                    teff1,teff2,radc(1),radc(2),
     &                    menv(1),menv(2),renv(1),renv(2),
     &                    ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                    bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                    epoch(2),bhspin(1),bhspin(2),
     &                    deltam1_bcm,deltam2_bcm,formation(1),
     &                    formation(2),binstate,mergertype,'bpp')
            CALL kick(kw,mass(k),mt,mass(3-k),ecc,sep,jorb,vk,k,
     &              rad(3-k),fallback,sigmahold,kick_info,disrupt,bkick)
            sigma = sigmahold !reset sigma after possible ECSN kick dist. Remove this if u want some kick link to the intial pulsar values...

            if(mass(3-k).lt.0.d0)then
               if(kstar(3-k).lt.0.d0) mt = mt-mass(3-k)
               if(kw.eq.13.and.mt.gt.mxns) kw = 14
               CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
               kstar(k) = kw
               mass(k) = mt
               epoch(k) = tphys - age
               kstar(3-k) = 15
               mass(3-k) = 0.d0
               coel = .true.
               binstate = 1
            endif
            if(ecc.gt.1.d0)then
               kstar(k) = kw
               mass(k) = mt
               epoch(k) = tphys - age
               goto 135
            endif
            tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
            oorb = twopi/tb
         endif
         if(kw.ne.kstar(k))then
            change = .true.
            if((kw.eq.13.or.kw.eq.14).and.kstar(k).le.12)then
               snova = .true.
            endif
            mass(k) = mt
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
         endif
*
*     Determine stellar evolution timescale for nuclear burning types.
*
         if(kw.le.9)then
            CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
            dtmi(k) = MIN(dt,dtr)
*           dtmi(k) = dtr
            dtmi(k) = MAX(1.0d-07,dtmi(k))
         else
            dtmi(k) = 1.0d+10
         endif
*        dtmi(k) = MAX((tn-age),1.0d-07)
*
* Save relevent solar quantities.
*
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         radx(k) = rm
         lumin(k) = lum
         teff(k) = 1000.d0*((1130.d0*lumin(k)/
     &                    (rad(k)**2.d0))**(1.d0/4.d0))
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
 90   continue
*
      do 100 , k = 1,2
         q(k) = mass(k)/mass(3-k)
         rol(k) = rl(q(k))*sep*(1.d0-ecc)
 100  continue
      if(rad(j1).gt.rol(j1)) radx(j1) = MAX(radc(j1),rol(j1))
      do 110 , k = 1,2
         ospin(k) = jspin(k)/(k2str(k)*(mass(k)-massc(k))*radx(k)*
     &              radx(k) + k3*massc(k)*radc(k)*radc(k))
 110  continue

      if(check_dtp.eq.1)then
          CALL checkstate(dtp,dtp_original,tsave,tphys,tphysf,
     &                      iplot,isave,binstate,evolve_type,
     &                      mass(1),mass(2),kstar(1),kstar(2),sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      radc(1),radc(2),menv(1),menv(2),
     &                      renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),
     &                      tacc(1),tacc(2),epoch(1),epoch(2),
     &                      bhspin(1),bhspin(2),teff1,teff2)
      endif
*
      if((isave.and.tphys.ge.tsave).or.iplot)then
          if(B_0(1).eq.0.d0)then !PK.
              b01_bcm = 0.d0
          elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
              b01_bcm = B_0(1)
          else
              b01_bcm = B(1)
          endif
          if(B_0(2).eq.0.d0)then
              b02_bcm = 0.d0
          elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
              b02_bcm = B_0(2)
          else
              b02_bcm = B(2)
          endif
          teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
          teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
          rrl1 = rad(1)/rol(1)
          rrl2 = rad(2)/rol(2)
          dt = MAX(dtm,1.0d-12)*1.0d+06
          if(j1.eq.1)then
              deltam1_bcm = (-1.0*dm1 - dms(1))/dt
              deltam2_bcm = (dm2 - dms(2))/dt
          else
              deltam1_bcm = (dm2 - dms(1))/dt
              deltam2_bcm = (-1.0*dm1 - dms(2))/dt
          endif
* Check if PISN occurred, and if so overwrite formation
          if(pisn_track(1).ne.0) formation(1) = pisn_track(1)
          if(pisn_track(2).ne.0) formation(2) = pisn_track(2)
          CALL writetab(ip,tphys,evolve_type,
     &                  mass(1),mass(2),kstar(1),kstar(2),
     &                  sep,tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bcm')
         if(isave) tsave = tsave + dtp
         if(output) write(*,*)'bcm3:',kstar(1),kstar(2),mass(1),
     & mass(2),rad(1),rad(2),ospin(1),ospin(2),jspin(1)
*     & mass(2),rad(1),rad(2),ospin(1),ospin(2),b01_bcm,b02_bcm,jspin(1)
      endif
*
      if(tphys.ge.tphysf.and..not.inspiral) goto 140
*
      if(change)then
         change = .false.
         evolve_type = 2.0
         rrl1 = rad(1)/rol(1)
         rrl2 = rad(2)/rol(2)
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                    mass(1),mass(2),kstar(1),kstar(2),
     &                    sep,tb,ecc,rrl1,rrl2,
     &                    aj(1),aj(2),tms(1),tms(2),
     &                    massc(1),massc(2),rad(1),rad(2),
     &                    mass0(1),mass0(2),lumin(1),lumin(2),
     &                    teff1,teff2,radc(1),radc(2),
     &                    menv(1),menv(2),renv(1),renv(2),
     &                    ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                    bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                    epoch(2),bhspin(1),bhspin(2),
     &                    deltam1_bcm,deltam2_bcm,formation(1),
     &                    formation(2),binstate,mergertype,'bpp')
      endif

*
* Test whether the primary still fills its Roche lobe.
*
      if(rad(j1).gt.rol(j1).and..not.snova)then
*
* Test for a contact system
*
         if(rad(j2).gt.rol(j2)) goto 130
         iter = iter + 1
         goto 8
      else
         evolve_type = 4.0
         rrl1 = rad(1)/rol(1)
         rrl2 = rad(2)/rol(2)
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
         dtm = 0.d0
         goto 4
      endif
*
 130  continue
*
* Contact system.
*
      coel = .true.
      binstate = 1
      CALL CONCATKSTARS(kstar(j1), kstar(j2), mergertype)
*
* If *1 or *2 is giant-like this will be common-envelope evolution.
*
      m1ce = mass(j1)
      m2ce = mass(j2)
      rrl1 = MIN(999.999d0,rad(1)/rol(1))
      rrl2 = MIN(999.999d0,rad(2)/rol(2))
*
      evolve_type = 5.0
      teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
      teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
      if(B_0(1).eq.0.d0)then !PK.
         b01_bcm = 0.d0
      elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
         b01_bcm = B_0(1)
      else
         b01_bcm = B(1)
      endif
      if(B_0(2).eq.0.d0)then
         b02_bcm = 0.d0
      elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
         b02_bcm = B_0(2)
      else
         b02_bcm = B(2)
      endif
      CALL writetab(jp,tphys,evolve_type,
     &              mass(1),mass(2),kstar(1),kstar(2),sep,
     &              tb,ecc,rrl1,rrl2,
     &              aj(1),aj(2),tms(1),tms(2),
     &              massc(1),massc(2),rad(1),rad(2),
     &              mass0(1),mass0(2),lumin(1),lumin(2),
     &              teff1,teff2,radc(1),radc(2),
     &              menv(1),menv(2),renv(1),renv(2),
     &              ospin(1),ospin(2),b01_bcm,b02_bcm,
     &              bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &              epoch(2),bhspin(1),bhspin(2),
     &              deltam1_bcm,deltam2_bcm,formation(1),
     &              formation(2),binstate,mergertype,'bpp')
*
      kcomp1 = kstar(j1)
      kcomp2 = kstar(j2)
*
      if(output) write(*,*)'coal r/rl1 & r/rl2 > 0',tphys,kcomp1,kcomp2,
     & m1ce,m2ce
*
      if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
         if(j1.eq.2)then
             switchedCE = .true.
         else
             switchedCE = .false.
         endif
         evolve_type = 7.d0
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),
     &                 kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel,j1,j2,
     &               vk,kick_info,formation(j1),formation(j2),sigmahold,
     &               bhspin(j1),bhspin(j2),binstate,mergertype,
     &               jp,tphys,switchedCE,rad,tms,evolve_type,disrupt,
     &               lumin,B_0,bacc,tacc,epoch,menv,renv,bkick,
     &               deltam1_bcm,deltam2_bcm)
         if(output) write(*,*)'coal1:',tphys,kstar(j1),kstar(j2),coel,
     & mass(j1),mass(j2)
         if(j1.eq.2.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif
         if(j1.eq.1.and.kcomp2.eq.13.and.kstar(j2).eq.15.and.
     &      kstar(j1).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j1) = formation(j2)
         endif
         com = .true.
         if(com.and..not.coel.and..not.disrupt)then
* if it went through common envelope
* did not disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system is still in binary
            binstate = 0
            mergertype = -1
         elseif(com.and..not.coel.and.disrupt)then
* if it went through common envelope
* and did disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system should be marked as disrupted
            binstate = 2
            mergertype = -1
         endif
* else it merged in the common envelope
         if(binstate.eq.1.d0)then
             sep = 0.d0
             tb = 0.d0
         elseif(binstate.eq.2.d0)then
             sep = -1.d0
             tb = -1.d0
         endif
      elseif(kstar(j2).ge.2.and.kstar(j2).le.9.and.kstar(j2).ne.7)then
         if(j1.eq.1)then
             switchedCE = .true.
         else
             switchedCE = .false.
         endif
         evolve_type = 7.d0
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),
     &                 kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
         CALL comenv(mass0(j2),mass(j2),massc(j2),aj(j2),jspin(j2),
     &               kstar(j2),mass0(j1),mass(j1),massc(j1),aj(j1),
     &               jspin(j1),kstar(j1),zpars,ecc,sep,jorb,coel,j2,j1,
     &               vk,kick_info,formation(j2),formation(j1),sigmahold,
     &               bhspin(j2),bhspin(j1),binstate,mergertype,
     &               jp,tphys,switchedCE,rad,tms,evolve_type,disrupt,
     &               lumin,B_0,bacc,tacc,epoch,menv,renv,bkick,
     &               deltam1_bcm,deltam2_bcm)
         if(output) write(*,*)'coal2:',tphys,kstar(j1),kstar(j2),coel,
     & mass(j1),mass(j2)
         if(j2.eq.2.and.kcomp1.eq.13.and.kstar(j1).eq.15.and.
     &      kstar(j2).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j2) = formation(j1)
         endif
         if(j2.eq.1.and.kcomp1.eq.13.and.kstar(j1).eq.15.and.
     &      kstar(j2).eq.13)then !PK.
* In CE the NS got switched around. Do same to formation.
            formation(j2) = formation(j1)
         endif
         com = .true.
         if(com.and..not.coel.and..not.disrupt)then
* if it went through common envelope
* did not disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system is still in binary
            binstate = 0
            mergertype = -1
         elseif(com.and..not.coel.and.disrupt)then
* if it went through common envelope
* and did disrupt (from one of the objects going SN)
* and did not merge in common envelope
* then system should be marked as disrupted
            binstate = 2
            mergertype = -1
         endif
* else it merged in the common envelope
         if(binstate.eq.1.d0)then
             sep = 0.d0
             tb = 0.d0
         elseif(binstate.eq.2.d0)then
             sep = -1.d0
             tb = -1.d0
         endif
      else
         CALL mix(mass0,mass,aj,kstar,zpars,bhspin)
      endif

      if(com)then
          evolve_type = 8.0
          mass1_bpp = mass(1)
          mass2_bpp = mass(2)
          if(kstar(1).eq.15) mass1_bpp = mass0(1)
          if(kstar(2).eq.15) mass2_bpp = mass0(2)
          rrl1 = MIN(rrl1,0.99d0)
          rrl2 = MIN(rrl2,0.99d0)
          teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
          teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
          if(B_0(1).eq.0.d0)then !PK.
             b01_bcm = 0.d0
          elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
             b01_bcm = B_0(1)
          else
             b01_bcm = B(1)
          endif
          if(B_0(2).eq.0.d0)then
             b02_bcm = 0.d0
          elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
             b02_bcm = B_0(2)
          else
             b02_bcm = B(2)
          endif

          CALL writetab(jp,tphys,evolve_type,
     &                  mass1_bpp,mass2_bpp,
     &                  kstar(1),kstar(2),sep,
     &                  tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bpp')
      endif
      epoch(1) = tphys - aj(1)
      epoch(2) = tphys - aj(2)
      if(.not.coel)then
*
* Next step should be made without changing the time.
*
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
*
* Need to confirm that RLO is over
*

         evolve_type = 4.0
         rrl1 = rad(1)/rol(1)
         rrl2 = rad(2)/rol(2)
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
         if(B_0(1).eq.0.d0)then !PK.
            b01_bcm = 0.d0
         elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
            b01_bcm = B_0(1)
         else
            b01_bcm = B(1)
         endif
         if(B_0(2).eq.0.d0)then
            b02_bcm = 0.d0
         elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
            b02_bcm = B_0(2)
         else
            b02_bcm = B(2)
         endif

         CALL writetab(jp,tphys,evolve_type,
     &                 mass(1),mass(2),kstar(1),kstar(2),sep,
     &                 tb,ecc,rrl1,rrl2,
     &                 aj(1),aj(2),tms(1),tms(2),
     &                 massc(1),massc(2),rad(1),rad(2),
     &                 mass0(1),mass0(2),lumin(1),lumin(2),
     &                 teff1,teff2,radc(1),radc(2),
     &                 menv(1),menv(2),renv(1),renv(2),
     &                 ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                 bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                 epoch(2),bhspin(1),bhspin(2),
     &                 deltam1_bcm,deltam2_bcm,formation(1),
     &                 formation(2),binstate,mergertype,'bpp')
         dtm = 0.d0
*
* Reset orbital parameters as separation may have changed.
*
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb

         goto 4
      endif
*
 135  continue
*
      sgl = .true.
      if(kstar(1).eq.13.and.mergemsp.eq.1.and.
     &   notamerger.eq.0)then
         s = (twopi*yearsc)/ospin(1)
         if(s.lt.0.03d0.and.B(1).gt.0.d0)then
            merge_mem = 1
         endif
      endif
      if(kstar(2).eq.13.and.mergemsp.eq.1.and.
     &   notamerger.eq.0)then
         s = (twopi*yearsc)/ospin(2)
         if(s.lt.0.03d0.and.B(2).gt.0.d0)then
            merge_mem = 1
         endif
      endif
      if(kstar(1).ne.15.or.kstar(2).ne.15)then
         if(com)then
            com = .false.
         else
            mass1_bpp = mass(1)
            mass2_bpp = mass(2)
*
* KB remove these 20 jul 23: we want the stars to not have mass
* if they are kstar=15
*
*            if(kstar(1).eq.15) mass1_bpp = mass0(1)
*            if(kstar(2).eq.15) mass2_bpp = mass0(2)
            if(coel)then
                evolve_type = 6.0
                teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
                teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
                if(B_0(1).eq.0.d0)then !PK.
                   b01_bcm = 0.d0
                elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                   b01_bcm = B_0(1)
                else
                   b01_bcm = B(1)
                endif
                if(B_0(2).eq.0.d0)then
                   b02_bcm = 0.d0
                elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                   b02_bcm = B_0(2)
                else
                   b02_bcm = B(2)
                endif
                CALL writetab(jp,tphys,evolve_type,
     &                        mass1_bpp,mass2_bpp,
     &                        kstar(1),kstar(2),0.d0,
     &                        0.d0,-1.d0,0.d0,ngtv,
     &                        aj(1),aj(2),tms(1),tms(2),
     &                        massc(1),massc(2),rad(1),rad(2),
     &                        mass0(1),mass0(2),lumin(1),lumin(2),
     &                        teff1,teff2,radc(1),radc(2),
     &                        menv(1),menv(2),renv(1),renv(2),
     &                        ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                        bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                        epoch(2),bhspin(1),bhspin(2),
     &                        deltam1_bcm,deltam2_bcm,formation(1),
     &                        formation(2),binstate,mergertype,'bpp')
            elseif(ecc.gt.1.d0)then
*
* Binary dissolved by a supernova or tides.
*
                evolve_type = 11.0
                binstate = 2
                mergertype = -1
                tb = -1.d0
                sep = -1.d0
                ecc = -1.d0
                teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
                teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
                if(B_0(1).eq.0.d0)then !PK.
                   b01_bcm = 0.d0
                elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                   b01_bcm = B_0(1)
                else
                   b01_bcm = B(1)
                endif
                if(B_0(2).eq.0.d0)then
                   b02_bcm = 0.d0
                elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                   b02_bcm = B_0(2)
                else
                   b02_bcm = B(2)
                endif
                CALL writetab(jp,tphys,evolve_type,
     &                        mass1_bpp,mass2_bpp,
     &                        kstar(1),kstar(2),sep,
     &                        tb,ecc,0.d0,ngtv2,
     &                        aj(1),aj(2),tms(1),tms(2),
     &                        massc(1),massc(2),rad(1),rad(2),
     &                        mass0(1),mass0(2),lumin(1),lumin(2),
     &                        teff1,teff2,radc(1),radc(2),
     &                        menv(1),menv(2),renv(1),renv(2),
     &                        ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                        bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                        epoch(2),bhspin(1),bhspin(2),
     &                        deltam1_bcm,deltam2_bcm,formation(1),
     &                        formation(2),binstate,mergertype,'bpp')
            else
                evolve_type = 9.0
                teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
                teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
                if(B_0(1).eq.0.d0)then !PK.
                   b01_bcm = 0.d0
                elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                   b01_bcm = B_0(1)
                else
                   b01_bcm = B(1)
                endif
                if(B_0(2).eq.0.d0)then
                   b02_bcm = 0.d0
                elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                   b02_bcm = B_0(2)
                else
                   b02_bcm = B(2)
                endif
                CALL writetab(jp,tphys,evolve_type,
     &                        mass1_bpp,mass2_bpp,
     &                        kstar(1),kstar(2),0.d0,
     &                        0.d0,0.d0,0.d0,ngtv,
     &                        aj(1),aj(2),tms(1),tms(2),
     &                        massc(1),massc(2),rad(1),rad(2),
     &                        mass0(1),mass0(2),lumin(1),lumin(2),
     &                        teff1,teff2,radc(1),radc(2),
     &                        menv(1),menv(2),renv(1),renv(2),
     &                        ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                        bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                        epoch(2),bhspin(1),bhspin(2),
     &                        deltam1_bcm,deltam2_bcm,formation(1),
     &                        formation(2),binstate,mergertype,'bpp')
            endif
         endif
         if(kstar(2).eq.15)then
            kmax = 1
            rol(2) = -1.d0*rad(2)
            dtmi(2) = tphysf
         elseif(kstar(1).eq.15)then
            kmin = 2
            rol(1) = -1.d0*rad(1)
            dtmi(1) = tphysf
         endif
* Makes sure coalesced NSs are reset. PK.
         if(kstar(1).eq.13.and.ecc.le.1.d0.and.pulsar.gt.0.and.
     &      notamerger.eq.0)then
            age = 0.d0
            epoch(1) = tphys
         endif
         if(kstar(2).eq.13.and.ecc.le.1.d0.and.pulsar.gt.0.and.
     &      notamerger.eq.0)then
            age = 0.d0
            epoch(2) = tphys
         endif
         ecc = -1.d0
         if(binstate.eq.2)then
*            Check if disrupted then we want sep=ecc=porb=-1
             sep = -1.d0
         elseif(binstate.eq.1)then
*            check if merge then sep=0
             tb = 0.d0
             sep = 0.d0
         endif
         dtm = 0.d0
         coel = .false.
         goto 4
      endif
*
 140  continue
*
      if(com)then
         com = .false.
      else
          mass1_bpp = mass(1)
          mass2_bpp = mass(2)

          !added by PA for systems that manage to reach here
          !without encountering write bpp at all
          if (jp<1) then
              evolve_type = 1.d0
              rrl1 = rad(1)/rol(1)
              rrl2 = rad(2)/rol(2)
              teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
              teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
              if(B_0(1).eq.0.d0)then !PK.
                 b01_bcm = 0.d0
              elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                 b01_bcm = B_0(1)
              else
                 b01_bcm = B(1)
              endif
              if(B_0(2).eq.0.d0)then
                 b02_bcm = 0.d0
              elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                 b02_bcm = B_0(2)
              else
                 b02_bcm = B(2)
              endif

              CALL writetab(jp,tphys,evolve_type,
     &                  mass(1),mass(2),kstar(1),kstar(2),sep,
     &                  tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bpp')
          endif
          
*          if(kstar(1).eq.15.and.bpp(jp,4).lt.15.0)then
*              mass1_bpp = mass0(1)
*          endif

*          if(kstar(2).eq.15.and.bpp(jp,5).lt.15.0)then
*              mass2_bpp = mass0(2)
*          endif
          
          if(coel)then
              evolve_type = 6.0
              teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
              teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
              if(B_0(1).eq.0.d0)then !PK.
                 b01_bcm = 0.d0
              elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                 b01_bcm = B_0(1)
              else
                 b01_bcm = B(1)
              endif
              if(B_0(2).eq.0.d0)then
                 b02_bcm = 0.d0
              elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                 b02_bcm = B_0(2)
              else
                 b02_bcm = B(2)
              endif
              CALL writetab(jp,tphys,evolve_type,
     &                      mass1_bpp,mass2_bpp,
     &                      kstar(1),kstar(2),0.d0,
     &                      0.d0,-1.d0,0.d0,ngtv,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      teff1,teff2,radc(1),radc(2),
     &                      menv(1),menv(2),renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin(1),bhspin(2),
     &                      deltam1_bcm,deltam2_bcm,formation(1),
     &                      formation(2),binstate,mergertype,'bpp')
          elseif(kstar(1).eq.15.and.kstar(2).eq.15)then
*
* Cases of accretion induced supernova or single star supernova.
* No remnant is left in either case.
*
              evolve_type = 9.0
              teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
              teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
              if(B_0(1).eq.0.d0)then !PK.
                 b01_bcm = 0.d0
              elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                 b01_bcm = B_0(1)
              else
                 b01_bcm = B(1)
              endif
              if(B_0(2).eq.0.d0)then
                 b02_bcm = 0.d0
              elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                 b02_bcm = B_0(2)
              else
                 b02_bcm = B(2)
              endif
              CALL writetab(jp,tphys,evolve_type,
     &                      mass1_bpp,mass2_bpp,
     &                      kstar(1),kstar(2),0.d0,
     &                      0.d0,0.d0,0.d0,ngtv2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      teff1,teff2,radc(1),radc(2),
     &                      menv(1),menv(2),renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin(1),bhspin(2),
     &                      deltam1_bcm,deltam2_bcm,formation(1),
     &                      formation(2),binstate,mergertype,'bpp')
          else
              evolve_type = 10.0
              !added by PA for systems that stop evolving halfway
              if(iter.ge.loop) evolve_type = 100.0
              rrl1 = rad(1)/rol(1)
              rrl2 = rad(2)/rol(2)
              teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
              teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
              if(B_0(1).eq.0.d0)then !PK.
                 b01_bcm = 0.d0
              elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
                 b01_bcm = B_0(1)
              else
                 b01_bcm = B(1)
              endif
              if(B_0(2).eq.0.d0)then
                 b02_bcm = 0.d0
              elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
                 b02_bcm = B_0(2)
              else
                 b02_bcm = B(2)
              endif
              CALL writetab(jp,tphys,evolve_type,
     &                      mass1_bpp,mass2_bpp,
     &                      kstar(1),kstar(2),sep,
     &                      tb,ecc,rrl1,rrl2,
     &                      aj(1),aj(2),tms(1),tms(2),
     &                      massc(1),massc(2),rad(1),rad(2),
     &                      mass0(1),mass0(2),lumin(1),lumin(2),
     &                      teff1,teff2,radc(1),radc(2),
     &                      menv(1),menv(2),renv(1),renv(2),
     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                      bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                      epoch(2),bhspin(1),bhspin(2),
     &                      deltam1_bcm,deltam2_bcm,formation(1),
     &                      formation(2),binstate,mergertype,'bpp')
          endif
      endif
*
*      CALL checkstate(dtp,dtp_original,tsave,tphys,tphysf,
*     &                      iplot,isave,binstate,evolve_type,
*     &                      mass(1),mass(2),kstar(1),kstar(2),sep,
*     &                      tb,ecc,rrl1,rrl2,
*     &                      aj(1),aj(2),tms(1),tms(2),
*     &                      massc(1),massc(2),rad(1),rad(2),
*     &                      mass0(1),mass0(2),lumin(1),lumin(2),
*     &                      radc(1),radc(2),menv(1),menv(2),
*     &                      renv(1),renv(2),
*     &                      ospin(1),ospin(2),b01_bcm,b02_bcm,
*     &                      bacc(1),bacc(2),
*     &                      tacc(1),tacc(2),epoch(1),epoch(2),
*     &                      bhspin(1),bhspin(2))
      if((isave.and.tphys.ge.tsave).or.iplot)then
          if(B_0(1).eq.0.d0)then !PK.
              b01_bcm = 0.d0
          elseif(B_0(1).gt.0.d0.and.B(1).eq.0.d0)then
              b01_bcm = B_0(1)
          else
              b01_bcm = B(1)
          endif
          if(B_0(2).eq.0.d0)then
              b02_bcm = 0.d0
          elseif(B_0(2).gt.0.d0.and.B(2).eq.0.d0)then
              b02_bcm = B_0(2)
          else
              b02_bcm = B(2)
          endif
          teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
          teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
          rrl1 = rad(1)/rol(1)
          rrl2 = rad(2)/rol(2)
          dt = MAX(dtm,1.0d-12)*1.0d+06
          if(j1.eq.1)then
              deltam1_bcm = (-1.0*dm1 - dms(1))/dt
              deltam2_bcm = (dm2 - dms(2))/dt
          else
              deltam1_bcm = (dm2 - dms(1))/dt
              deltam2_bcm = (-1.0*dm1 - dms(2))/dt
          endif
* Check if PISN occurred, and if so overwrite formation
          if(pisn_track(1).ne.0) formation(1) = pisn_track(1)
          if(pisn_track(2).ne.0) formation(2) = pisn_track(2)
          CALL writetab(ip,tphys,evolve_type,
     &                  mass(1),mass(2),kstar(1),kstar(2),
     &                  sep,tb,ecc,rrl1,rrl2,
     &                  aj(1),aj(2),tms(1),tms(2),
     &                  massc(1),massc(2),rad(1),rad(2),
     &                  mass0(1),mass0(2),lumin(1),lumin(2),
     &                  teff1,teff2,radc(1),radc(2),
     &                  menv(1),menv(2),renv(1),renv(2),
     &                  ospin(1),ospin(2),b01_bcm,b02_bcm,
     &                  bacc(1),bacc(2),tacc(1),tacc(2),epoch(1),
     &                  epoch(2),bhspin(1),bhspin(2),
     &                  deltam1_bcm,deltam2_bcm,formation(1),
     &                  formation(2),binstate,mergertype,'bcm')
         if(output) write(*,*)'bcm4:',kstar(1),kstar(2),mass(1),
     & mass(2),rad(1),rad(2),ospin(1),ospin(2),jspin(1),
     & tphys,tphysf
*     & mass(2),rad(1),rad(2),ospin(1),ospin(2),b01_bcm,b02_bcm,jspin(1),
         if(isave) tsave = tsave + dtp
         if(tphysf.le.0.d0)then
            ip = ip + 1
            do 145 , k = 1,38
               bcm(ip,k) = bcm(ip-1,k)
 145        continue
         endif
*
      elseif((kstar(1).eq.15.and.kstar(2).eq.15))then
         tphys = tphysf
         evolve_type = 10.0
         goto 135
      endif
      tphysfhold = tphysf
      tphysf = tphys
      if(sgl)then
         if(ecc.ge.0.d0.and.ecc.le.1.d0) ecc = -1.d0
         tb = -1.d0
      endif
      tb = tb*yeardy

      if(jp.ge.1000)then
         WRITE(*,*)' STOP: EVOLV2 ARRAY ERROR '
*         CALL exit(0)
*         STOP
      elseif(jp.ge.40)then
         WRITE(99,*)' EVOLV2 ARRAY WARNING ',mass1i,mass2i,tbi,ecci,jp
      endif
      if(iter.ge.loop)then
         WRITE(99,*)'ITER>=LOOP:',jp,tphys,tphysfhold,dtp,kstar,age,kst,
     & id1_pass,id2_pass,mass(1),mass(2),iter,loop
*         CALL exit(0)
*         STOP
      endif
      bcm(ip+1,1) = -1.0
      bpp(jp+1,1) = -1.0

      if(using_cmc.eq.0)then
          bcm_index_out = ip
          bpp_index_out = jp
          kick_info_out = kick_info
      endif
*

      END SUBROUTINE evolv2
***
