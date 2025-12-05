      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,sn,
     &                r2,fallback,sigmahold,kick_info,disrupt,bkick)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* Variables
* ---------
* kw: integer
*    Stellar type of the exploding star
* m1: real*8
*    Mass of the exploding star
* m1n: real*8
*    Mass of the compact remnant post-SN
* m2: real*8
*    Mass of the companion star
* ecc: real*8
*    Eccentricity of the binary pre-SN
* sep: real*8
*    Semi-major axis of the binary pre-SN
* sn: integer
*    Which star is going supernova (1 or 2)
* r2: real*8
*    Radius of the companion star
* fallback: real*8
*    Fallback mass fraction
* sigmahold: real*8
*    Original sigma value for the kick
* kick_info: real*8
*    Array with information about the supernova kicks (details below)
* bkick: real*8
*    Array with information about the kicks for CMC (details below)
* jorb: real*8, output
*    Total orbital angular momentum of the binary
* vk: real*8, output
*    Magnitude of the natal kick
* disrupt: logical, output
*    Whether the system is disrupted by the supernova

      integer kw,sn
      real*8 m1,m2,m1n,ecc,sep,jorb,vk,r2,fallback,sigmahold
      real*8 kick_info(2,18), bkick(20)
      logical disrupt

* Use one of the two kick prescriptions based on the kickflag
      if(kickflag.lt.0)then
* Original Kiel & Hurley 2009 prescription
         call kick_kiel(kw,m1,m1n,m2,ecc,sep,jorb,vk,sn,
     &                  r2,fallback,sigmahold,kick_info,disrupt,bkick)
      else
* New Pfahl et al. 2002 prescription
         call kick_pfahl(kw,m1,m1n,m2,ecc,sep,jorb,vk,sn,
     &                   r2,fallback,sigmahold,kick_info,disrupt,bkick)
      end if
      RETURN
      END


      SUBROUTINE kick_pfahl(kw,m1,m1n,m2,ecc,sep,jorb,vk,sn,r2,
     &                      fallback,sigmahold,kick_info,disrupt,bkick)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* This function is implements an entirely new kick prescription based
* on Appendix B of Pfahl et al. 2002
* https://ui.adsabs.harvard.edu/abs/2002ApJ...573..283P/abstract
* instead of Kiel & Hurley 2009.
*
* This prescription better accounts for how secondary stars are ejected
* from disrupted binaries and also corrects a few minor bugs in the
* implementation of K&H09.
*
* Specific kick magnitudes, angles, eccentric anomaly, and random seeds
* can be supplied in the initialization file with natal_kick_array, a
* (2,5) array with the first row being for sn=1 and second for sn=2
*
* kick_info[i,1]: sn of exploding star
* kick_info[i,2]: disrupted (0=no, 1=yes)
* kick_info[i,3]: magnitude of the natal kick
* kick_info[i,4-5]: phi and theta (in the frame of the exploding star)
* kick_info[i,6]: eccentric anamoly
* kick_info[i,7-9]: change in 3D systemic velocity of the binary, or the
*       change in 3D velocity of sn=1 if the system is disrupted
* kick_info[i,10]: magnitude of systemic velocity of the binary if bound
*       or magnitude of total velocity of sn=1 if disrupted,
*       accounting for both SNe
* kick_info[i,11-13]: change in 3D velocity of the sn=2 if system
*       is disrupted
* kick_info[i,14]: magnitude of velocity of sn=2 if disrupted,
*       accounting for both SNe
* kick_info[i,15]: First Euler angle of rotation of orbital plane after each SN
* kick_info[i,16]: Second Euler angle of rotation of orbital plane after each SN
* kick_info[i,17]: Third Euler angle of rotation of orbital plane after each SN
* kick_info[i,18]: random seed at the start of call to kick.f
*
* For cmc kick_info array is zero, not negative.
      integer kw,k,sn,safety,abskickflag

      real*8 m1,m2,m1n
      real*8 ecc,ecc_2,sep
      real*8 pi,twopi,yearsc,rsunkm,G_const
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mean_anom,ecc_anom,dif,der,del
      real*8 u1,u2,vk,vk2,v(4),s,sigmah
      real*8 theta,phi,sin_phi,cos_phi,sin_theta,cos_theta
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
      real*8 bkick(20),r2,jorb
      real*8 ecc_prev,a_prev,mtot,mtot_prev
      real*8 natal_kick(3), sep_vec(3), v_rel(3), v_rel_prev(3)
      real*8 a_prev_2, a_prev_3, cos_ecc_anom, sin_ecc_anom
      real*8 sqrt_1m_ecc_prev_2, sep_prev, prefactor, omega
      real*8 h_prev(3),h(3), h_hat(3), h_mag
      real*8 LRL_prev(3), LRL(3), e_hat(3)
      real*8 v_cm(3), v_sn(3), v_comp(3), v_inf_vec(3), v_inf
      real*8 v_sn_rot(3), v_comp_rot(3), v_cm_rot(3)
      real*8 h_cross_e_hat(3)
      real*8 thetaE, phiE, psiE
      real*8 psiplusphi, orbital_pivot_axis(3), unsigned_phi
      real*8 LRL_prev_dot_h, LRL_dot_h_prev, unsigned_psi
      real*8 disberg_mean
      integer i
      logical ECSN_or_USSN
* Output
      logical output,disrupt,collide
*
      real*8 kick_info(2,18)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      collide = .false.
      ECSN_or_USSN = .false.
      safety = 0
      abskickflag = ABS(kickflag)

* ----------------------------------------------------------------------
* -------------- Initialise variables and constants --------------------
* ----------------------------------------------------------------------

* Set up empty arrays and constants
      u1 = 0.d0
      u2 = 0.d0
      vk = 0.d0
      disberg_mean = 5.60d0
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Gravitational constant in units of km^3 / (Msun * s^2)
      G_const = 1.3271244d+11

* Set values for mean NS mass and mean ejecta as in Giacobbo & Mapelli 2020
      mean_mns = 1.2d0
      mean_mej = 9.0d0
* Set values for alpha and beta as in Bray & Eldridge 2016
      alphakick = 70.0d0
      betakick = 120.0d0

      if(using_cmc.eq.0)then
* check if we have supplied a randomseed for this SN from kick_info
* already
          if(natal_kick_array(sn,5).gt.0.d0)then
* if we have we need to run ran3 enough times until
* we are at the same state of the random number generator
* as we were before
              do while (natal_kick_array(sn,5).ne.idum1
     &                          .and.safety.le.20)
                  xx = RAN3(idum1)
                  safety = safety + 1
              end do
          endif
      endif
* save the current idum1
      natal_kick_array(sn,5) = idum1
      kick_info(sn,18) = idum1

* set the sn of the exploding object in the kick_info array
      kick_info(sn,1) = sn

* if the system was disrupted from the first supernova, marked as disrupted
      if(kick_info(1,2).eq.1) kick_info(2,2)=1

* sigma is negative for ECSN
      if((sigma.lt.0.d0).and.(abskickflag.eq.1.or.abskickflag.eq.5))then
         sigma = -1.d0*sigma
         ECSN_or_USSN = .true.
* for kick prescriptions other than default, revert to original sigma
      elseif((sigma.lt.0.d0).and.(abskickflag.gt.1))then
         sigma = sigmahold
      endif
      sigmah = sigma

* scale down BH kicks if bhsigmafrac is specified
      if(abskickflag.eq.1.or.abskickflag.eq.5)then
         if(kw.eq.14.or.(kw.eq.13.and.(m1n.ge.mxns)))then
            sigma = sigmah*bhsigmafrac
            disberg_mean = disberg_mean * bhsigmafrac
         endif
      endif


* ----------------------------------------------------------------------
* ----- Draw and scale a natal kick magnitude based on input dists -----
* ----------------------------------------------------------------------

* Before we draw the kick from the maxwellian and then scale it
* as desired, let us see if a pre-supplied natal kick magnitude
* was passed.
      if(natal_kick_array(sn,1).ge.0.d0)then
          vk = natal_kick_array(sn,1)
          vk2 = vk*vk
* per supplied kick value we mimic a call to random number generator
          xx = RAN3(idum1)
          xx = RAN3(idum1)
          xx = RAN3(idum1)
          xx = RAN3(idum1)
      else
* If no pre-supplied kick magnitude, we draw a kick from a distribution
* If the kickflag is 5 then use the log-normal distribution described
* by Disberg & Mandel 2025
          if(abskickflag.eq.5.and..not.ECSN_or_USSN)then
             call RandomLogNormal(disberg_mean,0.69d0,vk,idum1,twopi)
             vk2 = vk*vk
          else
* Otherwise use the Hobbs et al. 2005 Maxwellian distribution
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
            do 25 k = 1,2
                u1 = RAN3(idum1)
                u2 = RAN3(idum1)
                if(u1.gt.0.9999d0) u1 = 0.9999d0
                if(u2.gt.1.d0) u2 = 1.d0
* Generate two velocities from polar coordinates S & THETA.
                s = -2.d0*LOG(1.d0 - u1)
                s = sigma*SQRT(s)
                theta = twopi*u2
                v(2*k-1) = s*COS(theta)
                v(2*k) = s*SIN(theta)
25          continue
            vk2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
            vk = SQRT(vk2)
          endif

          if(abskickflag.eq.1.or.abskickflag.eq.5)then
* Limit BH kick with fallback mass fraction.
             if(kw.eq.14.and.bhflag.eq.0)then
                vk2 = 0.d0
                vk = 0.d0
             elseif(kw.eq.14.and.bhflag.eq.1)then
                 fallback = MIN(fallback,1.d0)
                 vk = MAX((1.d0-fallback)*vk,0.d0)
                 vk2 = vk*vk
             elseif(kw.eq.14.and.bhflag.eq.2)then
                vk = vk * mxns / m1n
                vk2 = vk*vk
             endif
          elseif(abskickflag.eq.2)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 1
             vk = vk * ((m1-m1n)/mean_mej) * (mean_mns/m1n)
             vk2 = vk*vk
          elseif(abskickflag.eq.3)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 2
             vk = vk * ((m1-m1n)/mean_mej)
             vk2 = vk*vk
          elseif(abskickflag.eq.4)then
* Use kick scaling from Bray & Eldridge 2016, Eq. 1
             vk = alphakick * ((m1-m1n)/m1n) + betakick
             vk2 = vk*vk
          endif

* If a massless remnant is produced then artificially set the kick to
* a large value that will cause a disruption. This avoids infinite kicks
* for prescriptions that scale with the mass of the remnant.
* See Issue #687
          if(m1n.le.0.d0)then
             vk = 10000.d0
             vk2 = vk * vk
          endif

      endif
      sigma = sigmah

* save natal kick velocity in the kick_info array and natal_kick_array
      kick_info(sn,3) = vk
      if(using_cmc.eq.0)then
          natal_kick_array(sn,1) = vk
      endif

* ----------------------------------------------------------------------
* --------- Now input or draw supernova natal kick angles --------------
* ----------------------------------------------------------------------

* Before we randomly draw a phi and theta for the natal kick,
* see if a pre-supplied set of phi/theta is passed
      if((natal_kick_array(sn,2).ge.(-90.d0)).and.
     &       (natal_kick_array(sn,2).le.(90.d0)))then
          phi = natal_kick_array(sn,2)*pi/180.d0
          sin_phi = SIN(phi)
* per supplied kick value we mimic a call to random number generator
          xx = RAN3(idum1)
          xx = RAN3(idum1)
      else
* CLR - Allow for a restricted opening angle for SN kicks
*       Only relevant for binaries, obviously
*       Default value for polar_kick_angle = 90.0
          bound = SIN((90.d0 - polar_kick_angle)*pi/180.d0)
          sin_phi = (1.d0-bound)*ran3(idum1) + bound
          phi = ASIN(sin_phi)
* MJZ - The constrained kick will hit at either the north
*       or south pole, so randomly choose the hemisphere
          if(RAN3(idum1).ge.0.5)then
            phi = -phi
            sin_phi = SIN(phi)
          endif
      endif
      cos_phi = COS(phi)

      if((natal_kick_array(sn,3).ge.(0.d0)).and.
     &       (natal_kick_array(sn,3).le.(360.d0)))then
          theta = natal_kick_array(sn,3)*pi/180.d0
* per supplied kick value we mimic a call to random number generator
           xx = RAN3(idum1)
      else
          theta = twopi*ran3(idum1)
      endif
      sin_theta = SIN(theta)
      cos_theta = COS(theta)

*     save theta and phi in the kick_info and
*     natal_kick_array
      kick_info(sn,4) = phi*180/pi
      kick_info(sn,5) = theta*180/pi
      if(using_cmc.eq.0)then
          natal_kick_array(sn,2) = phi*180/pi
          natal_kick_array(sn,3) = theta*180/pi
      endif

* create a vector for the natal kick
      natal_kick(1) = vk * cos_phi * cos_theta
      natal_kick(2) = vk * cos_phi * sin_theta
      natal_kick(3) = vk * sin_phi

* ----------------------------------------------------------------------
* ----- Natal kick all done, check for pre-disruption as quick exit ----
* ----------------------------------------------------------------------

* Check if the system is already not a bound binary
      if((sn.eq.2.and.kick_info(1,2).eq.1)
     &   .or.sep.le.0.or.ecc.lt.0
     &   .or.sep.ne.sep.or.ecc.ne.ecc)then
* if so, only apply kick to the current star
         disrupt = .true.
         kick_info(sn,2) = 1
         if(sn.eq.1)then
            kick_info(sn,7) = natal_kick(1)
            kick_info(sn,8) = natal_kick(2)
            kick_info(sn,9) = natal_kick(3)
         elseif(sn.eq.2)then
            kick_info(sn,11) = natal_kick(1)
            kick_info(sn,12) = natal_kick(2)
            kick_info(sn,13) = natal_kick(3)
         endif
         goto 78
      endif

* ----------------------------------------------------------------------
* ------ Draw or input mean anomaly, solve for eccentric anomaly -------
* ----------------------------------------------------------------------

* Find the initial separation by randomly choosing a mean anomaly.
* check is user supplied mean anomaly
      xx = RAN3(idum1)
      if((natal_kick_array(sn,4).ge.(0.d0)).and.
     &   (natal_kick_array(sn,4).le.(360.d0)))then

         mean_anom = natal_kick_array(sn,4) * pi / 180.d0
      else
         mean_anom = xx * twopi
      endif
* Solve Kepler's equation for the eccentric anomaly from mean anomaly
* https://en.wikipedia.org/wiki/Eccentric_anomaly
      ecc_anom = mean_anom

      if(mean_anom.eq.0.d0) goto 8

  9   dif = ecc_anom - ecc * SIN(ecc_anom) - mean_anom
      if(ABS(dif / mean_anom).le.1.0d-04) goto 8
      der = 1.d0 - ecc * COS(ecc_anom)
      del = dif/der
      ecc_anom = ecc_anom - del
      goto 9

 8    continue

* ----------------------------------------------------------------------
* ------ Calculate whether system disrupts and CM velocity change ------
* ----------------------------------------------------------------------

* Some helper variables for calculations below
      mtot_prev = m1 + m2
      mtot = m1n + m2
      ecc_prev = ecc
* Convert the separation to km
      a_prev = sep * rsunkm
      a_prev_2 = a_prev * a_prev
      a_prev_3 = a_prev_2 * a_prev
      cos_ecc_anom = COS(ecc_anom)
      sin_ecc_anom = SIN(ecc_anom)
      sqrt_1m_ecc_prev_2 = SQRT(1.d0 - ecc_prev * ecc_prev)

* Orbital frequency pre-SN (in 1/s)
      omega = SQRT(G_const * mtot_prev / a_prev_3)

* Separation vector before the supernova (in km)
      sep_vec(1) = a_prev * (cos_ecc_anom - ecc_prev)
      sep_vec(2) = a_prev * sqrt_1m_ecc_prev_2 * sin_ecc_anom
      sep_vec(3) = 0.d0
      call VectorMagnitude(sep_vec, sep_prev)

* Relative velocity vector before the supernova (in km/s)
      prefactor = omega * a_prev_2 / sep_prev
      v_rel_prev(1) = -prefactor * sin_ecc_anom
      v_rel_prev(2) = prefactor * sqrt_1m_ecc_prev_2 * cos_ecc_anom
      v_rel_prev(3) = 0.d0

* Specific angular momentum vector pre-SN (in km^2/s)
      call CrossProduct(sep_vec, v_rel_prev, h_prev)

* Laplace-Runge-Lenz vector pre-SN (unitless)
      call CrossProduct(v_rel_prev, h_prev, LRL_prev)
      do i = 1, 3
         LRL_prev(i) = LRL_prev(i) / (G_const * mtot_prev)
     &               - sep_vec(i) / sep_prev
      end do

* Calculate the new systemic velocity of the center of mass (in km/s)
      do i = 1, 3
         v_cm(i) = (-m2 * (m1 - m1n) / mtot_prev / mtot)
     &           * v_rel_prev(i)
     &           + (m1n / mtot * natal_kick(i))
      end do

* New velocity vectors after SN (in km/s)
      v_rel(1) = v_rel_prev(1) + natal_kick(1)
      v_rel(2) = v_rel_prev(2) + natal_kick(2)
      v_rel(3) = v_rel_prev(3) + natal_kick(3)

* Updated specific orbital angular momentum vector (in km^2/s)
      call CrossProduct(sep_vec, v_rel, h)

* Updated Laplace-Runge-Lenz vector (unitless)
      call CrossProduct(v_rel, h, LRL)
      DO i = 1, 3
         LRL(i) = LRL(i) / (G_const * mtot)
     &          - sep_vec(i) / sep_prev
      END DO

* Get the Euler angles from previous kick for the rotation matrix
      thetaE = kick_info(1,15) * pi / 180.d0
      phiE = kick_info(1,16) * pi / 180.d0
      psiE = kick_info(1,17) * pi / 180.d0

* Get the new eccentricity
      call VectorMagnitude(LRL, ecc)
      ecc_2 = ecc * ecc

* Set the new semi-major axis (back in Rsun now)
      call VectorMagnitude(h, h_mag)
      sep = h_mag * h_mag / (G_const * mtot * (1 - ecc_2)) / rsunkm

* ----------------------------------------------------------------------
* -------- Split based on whether this kick disrupts the system --------
* ----------------------------------------------------------------------

      if(ecc.gt.1.d0)then
* System is now disrupted
         disrupt = .true.
* Set that it is disrupted in the kick_info array
         kick_info(sn,2) = 1
         call VectorHat(LRL, ecc, e_hat)
         call VectorHat(h, h_mag, h_hat)
         call CrossProduct(h_hat, e_hat, h_cross_e_hat)

* Velocity at infinity (in km/s)
         v_inf = G_const * mtot / h_mag * sqrt(ecc_2 - 1.d0)
         do i = 1, 3
            v_inf_vec(i) = v_inf * ((-1.d0 * e_hat(i) / ecc)
     &          + SQRT(1 - 1.d0 / ecc_2) * h_cross_e_hat(i))
         end do

* Velocity of the star going supernova post-SN (in km/s)
         do i = 1, 3
            v_sn(i) = (m2 / mtot) * v_inf_vec(i) + v_cm(i)
         end do

* Velocity of the companion star post-SN (in km/s)
         do i = 1, 3
            v_comp(i) = -(m1n / mtot) * v_inf_vec(i) + v_cm(i)
         end do

* if second supernova, need to change basis to the original orbital plane
         if(sn.eq.2)then
            call ChangeBasis(v_sn, thetaE, phiE, psiE, v_sn_rot)
            call ChangeBasis(v_comp, thetaE, phiE, psiE, v_comp_rot)
         else
            v_sn_rot = v_sn
            v_comp_rot = v_comp
         endif

* save the velocities to the kick_info table
         if(sn.eq.1)then
            kick_info(sn,7) = v_sn_rot(1)
            kick_info(sn,8) = v_sn_rot(2)
            kick_info(sn,9) = v_sn_rot(3)
            kick_info(sn,11) = v_comp_rot(1)
            kick_info(sn,12) = v_comp_rot(2)
            kick_info(sn,13) = v_comp_rot(3)

            bkick(1) = float(sn)
            bkick(2) = kick_info(sn,7)
            bkick(3) = kick_info(sn,8)
            bkick(4) = kick_info(sn,9)
            bkick(5) = float(sn)
            bkick(6) = kick_info(sn,11)
            bkick(7) = kick_info(sn,12)
            bkick(8) = kick_info(sn,13)
*
         elseif(sn.eq.2)then
            kick_info(sn,11) = v_sn_rot(1)
            kick_info(sn,12) = v_sn_rot(2)
            kick_info(sn,13) = v_sn_rot(3)
            kick_info(sn,7) = v_comp_rot(1)
            kick_info(sn,8) = v_comp_rot(2)
            kick_info(sn,9) = v_comp_rot(3)

            bkick(5) = float(sn)
            bkick(6) = kick_info(sn,11)
            bkick(7) = kick_info(sn,12)
            bkick(8) = kick_info(sn,13)
            bkick(9) = float(sn)
            bkick(10) = kick_info(sn,7)
            bkick(11) = kick_info(sn,8)
            bkick(12) = kick_info(sn,9)
         endif

* lastly, check if this supernova results in a collision between stars
         call CollisionCheck(sep_vec, v_sn_rot, v_comp_rot, r2, collide)
* if it does, assume the supernova star plows right through the
* companion and obliterates it
         if(collide)then
            kick_info(sn,7) = v_sn_rot(1)
            kick_info(sn,8) = v_sn_rot(2)
            kick_info(sn,9) = v_sn_rot(3)
            kick_info(sn,11) = 0.d0
            kick_info(sn,12) = 0.d0
            kick_info(sn,13) = 0.d0
            bkick(6) = v_sn_rot(1)
            bkick(7) = v_sn_rot(2)
            bkick(8) = v_sn_rot(3)
            bkick(10) = 0.d0
            bkick(11) = 0.d0
            bkick(12) = 0.d0
            m2 = -1.d0*m2
         endif

         call AngleBetweenVectors(h, h_prev, thetaE)
         phiE = ran3(idum1) * twopi
         psiE = ran3(idum1) * twopi

* ----------------------------------------------------------------------
* The system is still bound
      else
* Record the mean anomaly in the arrays
         kick_info(sn,6) = mean_anom * 180 / pi
         if (using_cmc.eq.0) then
            natal_kick_array(sn,4) = mean_anom * 180 / pi
         endif

* Update the total orbital angular momentum (in Msun Rsun^2/yr)
         jorb = m1n * m2 / mtot * h_mag / rsunkm / rsunkm * yearsc

         if (sn.eq.2) then
            call ChangeBasis(v_cm, thetaE, phiE, psiE, v_cm_rot)
         else
            v_cm_rot = v_cm
         endif

*
* If system survives the SN, save the components of the change in
* centre-of-mass velocity
         kick_info(sn,7) = v_cm_rot(1)
         kick_info(sn,8) = v_cm_rot(2)
         kick_info(sn,9) = v_cm_rot(3)
         kick_info(sn,11) = 0
         kick_info(sn,12) = 0
         kick_info(sn,13) = 0

* 1st time with kick.
         if(bkick(1).le.0.d0)then
            bkick(1) = float(sn)
            bkick(2) = v_cm_rot(1)
            bkick(3) = v_cm_rot(2)
            bkick(4) = v_cm_rot(3)
* 2nd time with kick.
         elseif(bkick(5).le.0.d0)then
            bkick(5) = float(sn)
            bkick(6) = v_cm_rot(1)
            bkick(7) = v_cm_rot(2)
            bkick(8) = v_cm_rot(3)
* 2nd time with kick if already disrupted.
* MJZ - would this if statement ever be hit?
         elseif(bkick(5).gt.0.d0)then
            bkick(9) = float(sn)
            bkick(10) = v_cm_rot(1)
            bkick(11) = v_cm_rot(2)
            bkick(12) = v_cm_rot(3)
         endif
* In the impossible chance that the system is exactly parabolic...
         if(ecc.eq.1.d0.and.sn.eq.1)then
            kick_info(sn,7) = v_cm_rot(1)
            kick_info(sn,8) = v_cm_rot(2)
            kick_info(sn,9) = v_cm_rot(3)
            kick_info(sn,11) = -v_cm_rot(1)
            kick_info(sn,12) = -v_cm_rot(2)
            kick_info(sn,13) = -v_cm_rot(3)
            bkick(1) = float(sn)
            bkick(2) = v_cm_rot(1)
            bkick(3) = v_cm_rot(2)
            bkick(4) = v_cm_rot(3)
            bkick(5) = float(sn)
            bkick(6) = -v_cm_rot(1)
            bkick(7) = -v_cm_rot(2)
            bkick(8) = -v_cm_rot(3)
         elseif(ecc.eq.1.d0.and.sn.eq.2)then
            kick_info(sn,7) = -v_cm_rot(1)
            kick_info(sn,8) = -v_cm_rot(2)
            kick_info(sn,9) = -v_cm_rot(3)
            kick_info(sn,11) = v_cm_rot(1)
            kick_info(sn,12) = v_cm_rot(2)
            kick_info(sn,13) = v_cm_rot(3)
            bkick(5) = float(sn)
            bkick(6) = v_cm_rot(1)
            bkick(7) = v_cm_rot(2)
            bkick(8) = v_cm_rot(3)
            bkick(9) = float(sn)
            bkick(10) = -v_cm_rot(1)
            bkick(11) = -v_cm_rot(2)
            bkick(12) = -v_cm_rot(3)
         endif

* Update the Euler angles for the orbital plane rotation
         call AngleBetweenVectors(h, h_prev, thetaE)

* first two special cases for orbital A.M. remaining unchanged in angle
* since the cross product is not well defined in this case
         if(thetaE.eq.0.d0)then
            xx = ran3(idum1)
            phiE = twopi * xx
            xx = ran3(idum1)
* we can only calculate the angle if the eccentricity is nonzero
            if (ecc_prev.gt.0.d0.and.ecc.gt.0.d0)then
               call AngleBetweenVectors(LRL, LRL_prev, psiPlusPhi)
               psiE = psiPlusPhi - phiE
            else
               psiE = twopi * xx
            end if            
         else if(thetaE.eq.pi)then
            xx = ran3(idum1)
            phiE = twopi * xx
            xx = ran3(idum1)
* we can only calculate the angle if the eccentricity is nonzero
            if (ecc_prev.gt.0.d0.and.ecc.gt.0.d0)then
               call AngleBetweenVectors(LRL, LRL_prev, psiPlusPhi)
               psiE = phiE + psiPlusPhi
            else
               psiE = twopi * xx
            end if
* now we can actually use the cross product to get the pivot axis
         else
            call CrossProduct(h_prev, h, orbital_pivot_axis)

* first handled phiE, need to check ecc_prev is nonzero
* since otherwise LRL_prev is not well-defined
            xx = ran3(idum1)
            if(ecc_prev.eq.0.d0)then
               phiE = twopi * xx
            else
               call DotProduct(LRL_prev, h, LRL_prev_dot_h)
               call AngleBetweenVectors(LRL_prev, orbital_pivot_axis,
     &                                  unsigned_phi)
               if (LRL_prev_dot_h.ge.0.d0) then
                  phiE = unsigned_phi
               else
                  phiE = -unsigned_phi
               endif
            endif

* repeat for psi, now focusing on ecc instead of ecc_prev
            xx = ran3(idum1)
            if(ecc.eq.0.d0)then
               psiE = twopi * xx
            else
               call DotProduct(LRL, h_prev, LRL_dot_h_prev)
               call AngleBetweenVectors(LRL, orbital_pivot_axis,
     &                                  unsigned_psi)
               if (LRL_dot_h_prev.ge.0.d0) then
                  psiE = unsigned_psi
               else
                  psiE = -unsigned_psi
               endif
            endif
         endif
      endif

* TODO: Does Katie want to randomise the Psi angle same as COMPAS?

* save Euler angles in the kick_info array
      kick_info(sn,15) = thetaE * 180 / pi
      kick_info(sn,16) = phiE * 180 / pi
      kick_info(sn,17) = psiE * 180 / pi

* For systems that were distrupted in the first SN, skip to here
 78   continue
*
* Set systemic velocity magnitudes in the kick_info array
* For first SN, this should be identical to the magnitude
* of the three component vectors. For the second SN, this
* will be the systemic velocity relative to the initial frame.
      if(sn.eq.1)then
         kick_info(sn,10) = SQRT(kick_info(sn,7)*kick_info(sn,7) +
     &              kick_info(sn,8)*kick_info(sn,8) +
     &              kick_info(sn,9)*kick_info(sn,9))
         kick_info(sn,14) = SQRT(kick_info(sn,11)*kick_info(sn,11) +
     &              kick_info(sn,12)*kick_info(sn,12) +
     &              kick_info(sn,13)*kick_info(sn,13))
      elseif(sn.eq.2)then
         kick_info(sn,10) = SQRT(
     &      (kick_info(1,7)+kick_info(2,7))*
     &      (kick_info(1,7)+kick_info(2,7)) +
     &      (kick_info(1,8)+kick_info(2,8))*
     &      (kick_info(1,8)+kick_info(2,8)) +
     &      (kick_info(1,9)+kick_info(2,9))*
     &      (kick_info(1,9)+kick_info(2,9)))
         kick_info(sn,14) = SQRT(
     &      (kick_info(1,11)+kick_info(2,11))*
     &      (kick_info(1,11)+kick_info(2,11)) +
     &      (kick_info(1,12)+kick_info(2,12))*
     &      (kick_info(1,12)+kick_info(2,12)) +
     &      (kick_info(1,13)+kick_info(2,13))*
     &      (kick_info(1,13)+kick_info(2,13)))
      endif
      
      RETURN
      END

* ======================================================================
* ================== Old kick routine follows ==========================
* ======================================================================
***
      SUBROUTINE kick_kiel(kw,m1,m1n,m2,ecc,sep,jorb,vk,snstar,r2,
     &                     fallback,sigmahold,kick_info,disrupt,bkick)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* WARNINGS (from Tom Wagg)
* ------------------------
* Here there be dragons...this prescription seems to have some issues.
*     1. Natal kick strongly affects ejection velocity of secondaries
*     2. Some coordinate transformations seem to be incorrect
*     3. Criteria for collisions are a little peculiar
* I recommend that the kick_pfahl routine be used instead
* ----------------------------------------------------------------------
*
* Updated JRH kick routine by PDK (see Kiel & Hurley 2009).
*
* Here theta is the \omega angle within the HTP02 paper (thus phi is phi).
* Also, mu is the \nu angle in the HTP02 paper and omega its azimuth
*
* Produces a random SN kick.
* Evolution was added for binaries in which the kick creates
* an eccentricity of greater than unity (i.e. hyperbolic orbit).
*
* Specific kick magnitudes, angles, eccentric anomaly, and random seeds
* can be supplied in the initialization file with natal_kick_array, a
* (2,5) array with the first row being for snstar=1 and second for snstar=2
*
* SBC (September 2020)
* put bkick array back in for compatibility with CMC. bkick array is not used
* for COSMIC BSE in anyway
*
* MJZ/SBC (April 2020)
* kick_info is a (2,17) array that tracks information about the supernova
* kicks. This allows us to track the total change to the systemic
* velocity and the total change in the orbital plane tilt after both
* supernovae, as well as reproduce systems.
* The first row contains information about the first supernova that
* occurs, the second row the second supernova.
* Note that some values the second row will take into account the
* effect of the first SN (e.g., kick_info[2,10] is the total systemic
* velocity after both supernovae).
*
* kick_info[i,1]: snstar of exploding star
* kick_info[i,2]: disrupted (0=no, 1=yes)
* kick_info[i,3]: magnitude of the natal kick
* kick_info[i,4-5]: phi and theta (in the frame of the exploding star)
* kick_info[i,6]: eccentric anamoly
* kick_info[i,7-9]: change in 3D systemic velocity of the binary, or the
*       change in 3D velocity of snstar=1 if the system is disrupted
* kick_info[i,10]: magnitude of systemic velocity of the binary if bound
*       or magnitude of total velocity of snstar=1 if disrupted,
*       accounting for both SNe
* kick_info[i,11-13]: change in 3D velocity of the snstar=2 if system
*       is disrupted
* kick_info[i,14]: magnitude of velocity of snstar=2 if disrupted,
*       accounting for both SNe
* kick_info[i,15]: (total) tilt of the orbital plane after each SN
*       w.r.t. the original angular momentum axis after each SN
* kick_info[i,16]: azimuthal angle of the orbital plane w.r.t. spins
* kick_info[i,18]: random seed at the start of call to kick.f
*
* For cmc kick_info array is zero, not negative.
      integer kw,k,snstar,sn,safety,abskickflag

      real*8 m1,m2,m1n,mbi,mbf,mdif
      real*8 ecc,sep,sepn,jorb,ecc2
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 x_tilt,y_tilt,z_tilt
      real*8 mu,cmu,smu,omega,comega,somega
      real*8 cmu1,smu1,comega1,somega1
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 vs(3),v1,v2,v3
      real*8 mx1,mx2,r2
      real*8 sigmah,RotInvX
      real*8 signs,sigc,psins,psic,cpsins,spsins,cpsic,spsic
      real*8 csigns
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
      real*8 bkick(20)
* Output
      logical output,disrupt
*
      real*8 kick_info(2,18)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      safety = 0
      abskickflag = ABS(kickflag)

* Set up empty arrays and constants
      do k = 1,3
         vs(k) = 0.d0
      enddo
      u1 = 0.d0
      u2 = 0.d0
      vk = 0.d0
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05

* Set values for mean NS mass and mean ejecta as in Giacobbo & Mapelli 2020
      mean_mns = 1.2d0
      mean_mej = 9.0d0
* Set values for alpha and beta as in Bray & Eldridge 2016
      alphakick = 70.0d0
      betakick = 120.0d0

* find whether this is the first or second supernova
      if(using_cmc.eq.0)then
          if(kick_info(1,1).eq.0) sn=1
          if(kick_info(1,1).gt.0) sn=2
      else
          if(bkick(1).eq.0) sn=1
          if(bkick(1).gt.0) sn=2
      endif

      if(using_cmc.eq.0)then
* check if we have supplied a randomseed for this SN from kick_info
* already
          if(natal_kick_array(snstar,5).gt.0.d0)then
* if we have we need to run ran3 enough times until
* we are at the same state of the random number generator
* as we were before
              do while (natal_kick_array(snstar,5).ne.idum1
     &                          .and.safety.le.20)
                  xx = RAN3(idum1)
                  safety = safety + 1
              end do
          endif
      endif
* save the current idum1
      natal_kick_array(snstar,5) = idum1
      kick_info(sn,18) = idum1
* add a blank column (not used in this prescription)
      kick_info(sn,17) = 0.d0

* set the SNstar of the exploding object in the kick_info array
      kick_info(sn,1) = snstar

* if the system was disrupted from the first supernova, marked as disrupted
      if(kick_info(1,2).eq.1) kick_info(2,2)=1

* sigma is negative for ECSN
      if((sigma.lt.0.d0).and.(abskickflag.eq.1))then
         sigma = -1.d0*sigma
* for kick prescriptions other than default, revert to original sigma
      elseif((sigma.lt.0.d0).and.(abskickflag.gt.1))then
         sigma = sigmahold
      endif
      sigmah = sigma

* scale down BH kicks if bhsigmafrac is specified
      if(abskickflag.eq.1)then
         if(kw.eq.14.or.(kw.eq.13.and.(m1n.ge.mxns)))then
            sigma = sigmah*bhsigmafrac
         endif
      endif

*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then

* check is user supplied mean anomaly
         if((natal_kick_array(snstar,4).ge.(0.d0)).and.
     &       (natal_kick_array(snstar,4).le.(360.d0)))then

             mm = natal_kick_array(snstar,4)*pi/180.d0
* per supplied kick value we mimic a call to random number generator
             xx = RAN3(idum1)

* Solve for the eccentric anomaly from the mean anomaly
* https://en.wikipedia.org/wiki/Eccentric_anomaly
             em = mm

             if(mm.eq.0.d0) goto 3

  4          dif = em - ecc*SIN(em) - mm
             if(ABS(dif/mm).le.1.0d-04) goto 3
             der = 1.d0 - ecc*COS(em)
             del = dif/der
             em = em - del
             goto 4

         endif

         xx = RAN3(idum1)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
*
* Find the initial relative velocity vector.
* With a randomly selected quadrant of the orbit.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
* Before we draw the kick from the maxwellian and then scale it
* as desired, let us see if a pre-supplied natal kick maganitude
* was passed.
      if(natal_kick_array(snstar,1).ge.0.d0)then
          vk = natal_kick_array(snstar,1)
          vk2 = vk*vk
* per supplied kick value we mimic a call to random number generator
          xx = RAN3(idum1)
          xx = RAN3(idum1)
          xx = RAN3(idum1)
          xx = RAN3(idum1)
      else
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
          do 20 k = 1,2
             u1 = RAN3(idum1)
             u2 = RAN3(idum1)
             if(u1.gt.0.9999d0) u1 = 0.9999d0
             if(u2.gt.1.d0) u2 = 1.d0
* Generate two velocities from polar coordinates S & THETA.
             s = -2.d0*LOG(1.d0 - u1)
             s = sigma*SQRT(s)
             theta = twopi*u2
             v(2*k-1) = s*COS(theta)
             v(2*k) = s*SIN(theta)
 20          continue
          vk2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
          vk = SQRT(vk2)

          if(abskickflag.eq.1)then
* Limit BH kick with fallback mass fraction.
             if(kw.eq.14.and.bhflag.eq.0)then
                vk2 = 0.d0
                vk = 0.d0
             elseif(kw.eq.14.and.bhflag.eq.1)then
                 fallback = MIN(fallback,1.d0)
                 vk = MAX((1.d0-fallback)*vk,0.d0)
                 vk2 = vk*vk
             elseif(kw.eq.14.and.bhflag.eq.2)then
                vk = vk * mxns / m1n
                vk2 = vk*vk
             endif
          elseif(abskickflag.eq.2)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 1
             vk = vk * ((m1-m1n)/mean_mej) * (mean_mns/m1n)
             vk2 = vk*vk
          elseif(abskickflag.eq.3)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 2
             vk = vk * ((m1-m1n)/mean_mej)
             vk2 = vk*vk
          elseif(abskickflag.eq.4)then
* Use kick scaling from Bray & Eldridge 2016, Eq. 1
             vk = alphakick * ((m1-m1n)/m1n) + betakick
             vk2 = vk*vk
          endif

      endif
      sigma = sigmah

* save natal kick velocity in the kick_info array and natal_kick_array
      kick_info(sn,3) = vk
      if(using_cmc.eq.0)then
          natal_kick_array(snstar,1) = vk
      endif
* Before we randomly draw a phi and theta for the natal kick,
* see if a pre-supplied set of phi/theta is passed
      if((natal_kick_array(snstar,2).ge.(-90.d0)).and.
     &       (natal_kick_array(snstar,2).le.(90.d0)))then
          phi = natal_kick_array(snstar,2)*pi/180.d0
          sphi = SIN(phi)
* per supplied kick value we mimic a call to random number generator
          xx = RAN3(idum1)
          xx = RAN3(idum1)
      else
* CLR - Allow for a restricted opening angle for SN kicks
*       Only relevant for binaries, obviously
*       Default value for polar_kick_angle = 90.0
          bound = SIN((90.d0 - polar_kick_angle)*pi/180.d0)
          sphi = (1.d0-bound)*ran3(idum1) + bound
          phi = ASIN(sphi)
* MJZ - The constrained kick will hit at either the north
*       or south pole, so randomly choose the hemisphere
          if(RAN3(idum1).ge.0.5)then
            phi = -phi
            sphi = SIN(phi)
          endif
      endif
      cphi = COS(phi)

      if((natal_kick_array(snstar,3).ge.(0.d0)).and.
     &       (natal_kick_array(snstar,3).le.(360.d0)))then
          theta = natal_kick_array(snstar,3)*pi/180.d0
* per supplied kick value we mimic a call to random number generator
           xx = RAN3(idum1)
      else
          theta = twopi*ran3(idum1)
      endif
      stheta = SIN(theta)
      ctheta = COS(theta)

*     save theta and phi (exploding star frame) in the kick_info and
*     natal_kick_array
      kick_info(sn,4) = phi*180/pi
      kick_info(sn,5) = theta*180/pi
      if(using_cmc.eq.0)then
          natal_kick_array(snstar,2) = phi*180/pi
          natal_kick_array(snstar,3) = theta*180/pi
      endif

* If the system is already disrupted, apply this kick only to the
* exploding star, and skip ahead.
      if(sn.eq.2.and.kick_info(1,2).eq.1)then
         if(snstar.eq.1)then
            kick_info(sn,7) = vk*COS(theta)*SIN(phi)
            kick_info(sn,8) = vk*SIN(theta)*SIN(phi)
            kick_info(sn,9) = vk*COS(phi)
         elseif(snstar.eq.2)then
            kick_info(sn,11) = vk*COS(theta)*SIN(phi)
            kick_info(sn,12) = vk*SIN(theta)*SIN(phi)
            kick_info(sn,13) = vk*COS(phi)
         endif
         GOTO 73
      endif

* CLR - if the orbit has already been kicked, then any polar kick
*       needs to be tilted as well (since L_hat and S_hat are no longer
*       aligned).
* MJZ - to track the total angular change in the binary's orbital plane,
*       we use both the value of \mu and \omega from the first SN.
*       We first rotate by \mu about the x-axis, then by \omega
*       about the z-axis.
*       We do this when the system has already had one SN (sn=2)
*  NOTE: This prescription does not account for realignment between SNe!
      if(sn.eq.2)then
        cmu = COS(kick_info(1,15)*pi/180)
        smu = SIN(kick_info(1,15)*pi/180)
        comega = COS(kick_info(1,16)*pi/180)
        somega = SIN(kick_info(1,16)*pi/180)

        x_tilt = ctheta*cphi*comega + smu*sphi*somega -
     &                   cmu*cphi*stheta*somega
        y_tilt = cmu*cphi*comega*stheta + ctheta*cphi*somega -
     &                   comega*smu*sphi
        z_tilt = cmu*sphi + cphi*smu*stheta

        phi = ASIN(z_tilt)
        sphi = z_tilt
        cphi = COS(phi)
        theta = ATAN(y_tilt/x_tilt)
        stheta = SIN(theta)
        ctheta = COS(theta)
      endif

      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90

*
* Determine the magnitude of the new relative velocity.
* CLR - fixed a spurious minus sign in the parenthesis here; only
*       relevant for eccentric orbits
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha+stheta*cphi*calpha)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      v3 = vk2*stheta*stheta*cphi*cphi
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(0.d0,ecc2)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the polar angle (mu) and azimuthal angle (omega)
* between the new and old orbital angular momentum vectors
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
      smu = SIN(mu)

      comega = (vr*salpha-vk*ctheta*cphi)/SQRT(v3 + v2)
      omega = ACOS(comega)
      somega = SIN(omega)

* Write angle between initial and current orbital angular momentum vectors
* and mean anomaly if system is still bound
      if(sn.eq.1.and.ecc.le.1)then
        kick_info(sn,15) = mu*180/pi
        kick_info(sn,16) = omega*180/pi
        kick_info(sn,6) = mm*180/pi
        if(using_cmc.eq.0)then
            natal_kick_array(snstar,4) = mm*180/pi
        endif
      elseif(sn.eq.2)then
* MJZ - Here we calculate the total change in the orbital plane
*       from both SN. Note that these angles mu and omega are in
*       typical spherical coordinates rather than colateral coordinates,
*       so the rotations are slightly different than above.
*       We rotate about z-axis by omega1 then y-axis by mu1
        cmu1 = COS(kick_info(1,15)*pi/180)
        smu1 = SIN(kick_info(1,15)*pi/180)
        comega1 = COS(kick_info(1,16)*pi/180)
        somega1 = SIN(kick_info(1,16)*pi/180)

        x_tilt = cmu1*comega*comega1*smu + cmu*smu1
     &                - cmu1*smu*somega*somega1
        y_tilt = comega1*smu*somega + comega*smu*somega1
        z_tilt = cmu*cmu1 + smu*smu1*somega*somega1
     &               - comega*comega1*smu*smu1

        kick_info(sn,15) = ACOS(z_tilt)*180/pi
* If both kicks were 0 then x_tilt for SN2 will be 0 since smu=0
        if(x_tilt.eq.0)then
          kick_info(sn,16) = 0
        else
          kick_info(sn,16) = ATAN(y_tilt/x_tilt)*180/pi
        endif
        kick_info(sn,6) = mm*180/pi
        if(using_cmc.eq.0)then
            natal_kick_array(snstar,4) = mm*180/pi
        endif

      endif

* Determine if orbit becomes hyperbolic.
 90   continue

* Calculate systemic velocity
      mx1 = vk*m1n/(m1n+m2)
      mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
      vs(1) = mx1*ctheta*cphi + mx2*salpha
      vs(2) = mx1*stheta*cphi + mx2*calpha
      vs(3) = mx1*sphi
*
* See if system was disrupted
      if(ecc.gt.1.d0)then
         disrupt = .true.
* Set that it is disrupted in the kick_info array
         kick_info(sn,2) = 1
*
*************************
* Kiel & Hurley method: * Similar to Belczynski et al. (2008) method.
*************************
* Find the semimajor axis and semiminor axis.
*         sepn = hn2/(gmrkm*(m1n+m2)*(ecc*ecc - 1.d0))
         sepn = -sep
*         bb = sqrt(ecc*ecc-1.d0)*sepn
* Direction for hyperbolic velocity in x (csig) and y (ssig).
         csigns = 1.d0/ecc
* Calculate the velocity magnitude at infinity for hyperbolic orbit.
         v1 = SQRT((gmrkm*(m1n+m2))/sepn)
         signs = ACOS(MIN(1.d0,csigns))
         sigc = signs
* Calculating position of NS compared to companion for
* rotation calculation around the z-axis.
         semilatrec = gmrkm*(m1n+m2)
         semilatrec = hn2/semilatrec
         cangleofdeath = (1.d0/ecc)*((semilatrec/r) - 1.d0)
         angleofdeath = ACOS(MIN(1.d0,cangleofdeath))
* Find which hemisphere (in x) the NS/companion originated in...
         psins = angleofdeath
         if((stheta*cphi - calpha).gt.0.d0) psins = -psins
         psic = psins
* Accounting for rotation in x-y plane due to SN event.
         cpsins = COS(signs-psins)
         spsins = SIN(signs-psins)
         cpsic = COS(sigc-psic)
         spsic = SIN(sigc-psic)
* Inverse rotation matrix accounting for our assumed alignment of
* the post-SN orbital angular mometum with the z-axis.
*         RotInvX = vr*salpha + vk*(sphi - ctheta*cphi)
*         mbi = RotInvX
*         RotInvZ = vr*salpha - vk*(ctheta*cphi+sphi)
*         RotInvX = RotInvX/SQRT(RotInvX*RotInvX + RotInvZ*RotInvZ)
*         RotInvZ = RotInvZ/SQRT(mbi*mbi + RotInvZ*RotInvZ)
         RotInvX = cmu
*
         mbi = m1+m2
         mbf = m1n+m2
         mdif = m1 - m1n
* Energy calculation, for interest - should be positive for unbound system.
         energy = vn2/2.d0 - gmrkm*(m1n+m2)/r
*
* Set values in the kick_info array for this disrupted system
* First specify that the system was disrupted from this SN
         kick_info(sn,2) = 1
* Now save the velocity of each star. Note that kick_info[4-6] is for snstar=1.
         if(snstar.eq.1)then
            kick_info(sn,7) = ((m1n/mbf)*vk*ctheta*cphi +
     &           (mdif*m2)/(mbi*mbf)*vr*salpha) -
     &           v1*(spsins)*RotInvX*(m2/mbf)
*
            kick_info(sn,8) = (m1n/mbf)*vk*stheta*cphi +
     &           (mdif*m2)/(mbi*mbf)*vr*calpha -
     &           v1*(cpsins)*(m2/mbf)
*
            kick_info(sn,9) = (m1n/mbf)*vk*sphi
*
            kick_info(sn,11) = ((mdif*m2)/(mbi*mbf)*vr*salpha +
     &           (m1n/mbf)*vk*ctheta*cphi) +
     &           v1*(spsic)*(m1n/mbf)*RotInvX
*
            kick_info(sn,12) = (mdif*m2)/(mbi*mbf)*vr*calpha +
     &           (m1n/mbf)*vk*stheta*cphi +
     &           v1*(cpsic)*(m1n/mbf)
*
            kick_info(sn,13) = (m1n/mbf)*vk*sphi

            bkick(1) = float(snstar)
            bkick(2) = kick_info(sn,7)
            bkick(3) = kick_info(sn,8)
            bkick(4) = kick_info(sn,9)
            bkick(5) = float(snstar)
            bkick(6) = kick_info(sn,11)
            bkick(7) = kick_info(sn,12)
            bkick(8) = kick_info(sn,13)

            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  kick_info(sn,7) = vs(1)
                  kick_info(sn,8) = vs(2)
                  kick_info(sn,9) = vs(3)
                  kick_info(sn,11) = 0.d0
                  kick_info(sn,12) = 0.d0
                  kick_info(sn,13) = 0.d0
                  bkick(2) = vs(1)
                  bkick(3) = vs(2)
                  bkick(4) = vs(3)
                  bkick(6) = 0.d0
                  bkick(7) = 0.d0
                  bkick(8) = 0.d0
                  m2 = -1.d0*m2
               endif
            endif
*
         elseif(snstar.eq.2)then
            kick_info(sn,11) = ((m1n/mbf)*vk*ctheta*cphi +
     &           (mdif*m2)/(mbi*mbf)*vr*salpha) -
     &           v1*(spsins)*RotInvX*(m2/mbf)
*
            kick_info(sn,12) = (m1n/mbf)*vk*stheta*cphi +
     &           (mdif*m2)/(mbi*mbf)*vr*calpha -
     &           v1*(cpsins)*(m2/mbf)
*
            kick_info(sn,13) = (m1n/mbf)*vk*sphi
*
            kick_info(sn,7) = ((mdif*m2)/(mbi*mbf)*vr*salpha +
     &           (m1n/mbf)*vk*ctheta*cphi) +
     &           v1*(spsic)*(m1n/mbf)*RotInvX
*
            kick_info(sn,8) = (mdif*m2)/(mbi*mbf)*vr*calpha +
     &           (m1n/mbf)*vk*stheta*cphi +
     &           v1*(cpsic)*(m1n/mbf)
*
            kick_info(sn,9) = (m1n/mbf)*vk*sphi

            bkick(5) = float(snstar)
            bkick(6) = kick_info(sn,11)
            bkick(7) = kick_info(sn,12)
            bkick(8) = kick_info(sn,13)
            bkick(9) = float(snstar)
            bkick(10) = kick_info(sn,7)
            bkick(11) = kick_info(sn,8)
            bkick(12) = kick_info(sn,9)

            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  kick_info(sn,7) = vs(1)
                  kick_info(sn,8) = vs(2)
                  kick_info(sn,9) = vs(3)
                  kick_info(sn,11) = 0.d0
                  kick_info(sn,12) = 0.d0
                  kick_info(sn,13) = 0.d0
                  bkick(6) = vs(1)
                  bkick(7) = vs(2)
                  bkick(8) = vs(3)
                  bkick(10) = 0.d0
                  bkick(11) = 0.d0
                  bkick(12) = 0.d0
                  m2 = -1.d0*m2
               endif
            endif
*
         endif
         ecc = MIN(ecc,99.99d0)
      endif
*
* If system survives the SN, save the components of the change in
* centre-of-mass velocity
      if(ecc.lt.1.d0)then
         kick_info(sn,7) = vs(1)
         kick_info(sn,8) = vs(2)
         kick_info(sn,9) = vs(3)
         kick_info(sn,11) = 0
         kick_info(sn,12) = 0
         kick_info(sn,13) = 0
      endif

      if(ecc.lt.1.d0)then
*         if(ecc.eq.1.d0.or.ecc.lt.0.d0) m2 = -1.d0 * m2
* 1st time with kick.
         if(bkick(1).le.0.d0)then
            bkick(1) = float(snstar)
            bkick(2) = vs(1)
            bkick(3) = vs(2)
            bkick(4) = vs(3)
* 2nd time with kick.
         elseif(bkick(5).le.0.d0)then
            bkick(5) = float(snstar)
            bkick(6) = vs(1)
            bkick(7) = vs(2)
            bkick(8) = vs(3)
* 2nd time with kick if already disrupted.
* MJZ - would this if statement ever be hit?
         elseif(bkick(5).gt.0.d0)then
            bkick(9) = float(snstar)
            bkick(10) = vs(1)
            bkick(11) = vs(2)
            bkick(12) = vs(3)
         endif
      endif

* In the impossible chance that the system is exactly parabolic...
      if(ecc.eq.1.d0.and.snstar.eq.1)then
         kick_info(sn,7) = vs(1)
         kick_info(sn,8) = vs(2)
         kick_info(sn,9) = vs(3)
         kick_info(sn,11) = -vs(1)
         kick_info(sn,12) = -vs(2)
         kick_info(sn,13) = -vs(3)
         bkick(1) = float(snstar)
         bkick(2) = vs(1)
         bkick(3) = vs(2)
         bkick(4) = vs(3)
         bkick(5) = float(snstar)
         bkick(6) = -vs(1)
         bkick(7) = -vs(2)
         bkick(8) = -vs(3)
      elseif(ecc.eq.1.d0.and.snstar.eq.2)then
         kick_info(sn,7) = -vs(1)
         kick_info(sn,8) = -vs(2)
         kick_info(sn,9) = -vs(3)
         kick_info(sn,11) = vs(1)
         kick_info(sn,12) = vs(2)
         kick_info(sn,13) = vs(3)
         bkick(5) = float(snstar)
         bkick(6) = vs(1)
         bkick(7) = vs(2)
         bkick(8) = vs(3)
         bkick(9) = float(snstar)
         bkick(10) = -vs(1)
         bkick(11) = -vs(2)
         bkick(12) = -vs(3)
      endif

* Uncomment to randomly rotate system velocities
      CALL randomness3(idum1,
     &            kick_info(sn,7),kick_info(sn,8),kick_info(sn,9),
     &            kick_info(sn,11),kick_info(sn,12),kick_info(sn,13))
*
* Put a cap on the eccentricity
      if(ecc.gt.99.9d0) ecc = 99.9d0

* For systems that were distrupted in the first SN, skip to here
 73   continue
*
* Set systemic velocity magnitudes in the kick_info array
* For first SN, this should be identical to the magnitude
* of the three component vectors. For the second SN, this
* will be the systemic velocity relative to the initial frame.
      if(sn.eq.1)then
         kick_info(sn,10) = SQRT(kick_info(sn,7)*kick_info(sn,7) +
     &              kick_info(sn,8)*kick_info(sn,8) +
     &              kick_info(sn,9)*kick_info(sn,9))
         kick_info(sn,14) = SQRT(kick_info(sn,11)*kick_info(sn,11) +
     &              kick_info(sn,12)*kick_info(sn,12) +
     &              kick_info(sn,13)*kick_info(sn,13))
      elseif(sn.eq.2)then
         kick_info(sn,10) = SQRT(
     &      (kick_info(1,7)+kick_info(2,7))*
     &      (kick_info(1,7)+kick_info(2,7)) +
     &      (kick_info(1,8)+kick_info(2,8))*
     &      (kick_info(1,8)+kick_info(2,8)) +
     &      (kick_info(1,9)+kick_info(2,9))*
     &      (kick_info(1,9)+kick_info(2,9)))
         kick_info(sn,14) = SQRT(
     &      (kick_info(1,11)+kick_info(2,11))*
     &      (kick_info(1,11)+kick_info(2,11)) +
     &      (kick_info(1,12)+kick_info(2,12))*
     &      (kick_info(1,12)+kick_info(2,12)) +
     &      (kick_info(1,13)+kick_info(2,13))*
     &      (kick_info(1,13)+kick_info(2,13)))
      endif

* Write some output information
*
      RETURN
      END
***
      SUBROUTINE randomness3(idum,vx1,vy1,vz1,vx2,vy2,vz2)
*
      implicit none
*
      INTEGER idum
      real ran3
      external ran3
      REAL*8 vx1,vy1,vz1,alpha,gamma,beta,pi,twopi,vx2,vy2,vz2
      REAL*8 cg,sg,ca,sa,cb,sb,vx1s,vy1s,vz1s,vx2s,vy2s,vz2s
*
* Introduce random orientation of binary system to the Galaxy/GC... u bastard.
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
      alpha = twopi*ran3(idum) - pi
      gamma = twopi*ran3(idum) - pi
      beta = pi*ran3(idum) - (pi/2.d0)
*      write(80,*)alpha,gamma,beta
      cg = COS(gamma)
      sg = SIN(gamma)
      ca = COS(alpha)
      sa = SIN(alpha)
      cb = COS(beta)
      sb = SIN(beta)
*
* Randomized orientation
*
      vx1s = vx1
      vy1s = vy1
      vz1s = vz1
      vx2s = vx2
      vy2s = vy2
      vz2s = vz2
* theta = b, phi = a, g = psi; pitch-roll-yaw; z-y-x axis rotations
      vx1 = vx1s*(cb*ca)
      vx1 = vx1 + vy1s*(cb*sa)
      vx1 = vx1 - vz1s*sb
*
      vy1 = vx1s*(sg*sb*ca - cg*sa)
      vy1 = vy1 + vy1s*(sg*sb*sa + cg*ca)
      vy1 = vy1 + vz1s*(cb*sg)
*
      vz1 = vx1s*(cg*sb*ca + sg*sa)
      vz1 = vz1 + vy1s*(cg*sb*sa - sg*ca)
      vz1 = vz1 + vz1s*(cb*cg)
**
**
      vx2 = vx2s*(cb*ca)
      vx2 = vx2 + vy2s*(cb*sa)
      vx2 = vx2 - vz2s*sb
*
      vy2 = vx2s*(sg*sb*ca - cg*sa)
      vy2 = vy2 + vy2s*(sg*sb*sa + cg*ca)
      vy2 = vy2 + vz2s*(cb*sg)
*
      vz2 = vx2s*(cg*sb*ca + sg*sa)
      vz2 = vz2 + vy2s*(cg*sb*sa - sg*ca)
      vz2 = vz2 + vz2s*(cb*cg)
**
*
      RETURN
      END
*



* ======================================================================
* ================== Vector helper functions follow ====================
* ======================================================================

      SUBROUTINE ChangeBasis(Vector, ThetaE, PhiE, PsiE, Result)
* Redefine a vector from one coordinate basis to another using Euler Angles
* Vector is the input vector in the new basis (X', Y', Z')
* Result is the transformed vector in the original basis (X, Y, Z)
* ThetaE, PhiE, PsiE are the Euler angles
  
      real*8 Vector(3), Result(3)
      real*8 ThetaE, PhiE, PsiE
      real*8 cTheta, cPhi, cPsi, sTheta, sPhi, sPsi
      real*8 rotationMatrix(3,3)
      integer i, j
    
* define trigonometric values
      cTheta = COS(ThetaE)
      sTheta = SIN(ThetaE)
      cPhi   = COS(PhiE)
      sPhi   = SIN(PhiE)
      cPsi   = COS(PsiE)
      sPsi   = SIN(PsiE)
    
* define the Rotation Matrix
      rotationMatrix(1,1) = cPhi * cPsi - sPhi * cTheta * sPsi
      rotationMatrix(1,2) = -cPhi * sPsi - sPhi * cTheta * cPsi
      rotationMatrix(1,3) = sTheta * sPhi
      rotationMatrix(2,1) = sPhi * cPsi + cPhi * cTheta * sPsi
      rotationMatrix(2,2) = -sPhi * sPsi + cPhi * cTheta * cPsi
      rotationMatrix(2,3) = -sTheta * cPhi
      rotationMatrix(3,1) = sTheta * sPsi
      rotationMatrix(3,2) = sTheta * cPsi
      rotationMatrix(3,3) = cTheta
    
* initialize the result to zero
      DO i = 1, 3
         Result(i) = 0.0D0
      END DO
    
* apply rotation to the vector
      DO i = 1, 3
         DO j = 1, 3
            Result(i) = Result(i) + Vector(j) * rotationMatrix(i, j)
         END DO
      END DO
    
      RETURN
      END


      SUBROUTINE CrossProduct(A, B, C)
* This function computes the cross product of two vectors A and B
* A, B are input vectors of dimension 3
* C is the resulting vector, also of dimension 3

      real*8 A(3), B(3), C(3)

* Calculate each component of the cross product
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)

      RETURN
      END



      SUBROUTINE VectorMagnitude(A, magnitude)
* This function computes the magnitude of a vector A
* A is the input vector of dimension 3

      real*8 A(3)
      real*8 magnitude

* Calculate the magnitude of the vector
      magnitude = SQRT(A(1) * A(1) + A(2) * A(2) + A(3) * A(3))

      RETURN
      END

      SUBROUTINE VectorHat(A, A_mag, A_hat)
* This function computes the unit vector of a vector A
* A is the input vector of dimension 3
* A_mag is the magnitude of the vector A
* A_hat is the resulting unit vector, also of dimension 3

      real*8 A(3), A_hat(3), A_mag

* Calculate the unit vector
      A_hat(1) = A(1) / A_mag
      A_hat(2) = A(2) / A_mag
      A_hat(3) = A(3) / A_mag
    
      RETURN
      END

      SUBROUTINE DotProduct(A, B, dot)
* This function computes the dot product of two vectors A and B
* A, B are input vectors of dimension 3
* dot is the resulting scalar

      real*8 A(3), B(3), dot

* Calculate the dot product
      dot = A(1) * B(1) + A(2) * B(2) + A(3) * B(3)

      RETURN
      END

      SUBROUTINE AngleBetweenVectors(A, B, angle)
* This function computes the angle between two vectors A and B
* A, B are input vectors of dimension 3
* angle is the resulting angle in radians

      real*8 A(3), B(3), angle
      real*8 dot, magA, magB

* Calculate the dot product of the two vectors
      call DotProduct(A, B, dot)

* Calculate the magnitudes of the two vectors
      call VectorMagnitude(A, magA)
      call VectorMagnitude(B, magB)

* Calculate the angle between the two vectors
      angle = ACOS(dot / (magA * magB))

      RETURN
      END


      SUBROUTINE CollisionCheck(sep_vec, v1, v2, r2, collide)
* This function checks if two stars collide
* It assumes that the compact object is a point mass, that the
* velocities are constant and the secondary star is a sphere, radius r2
*
* Method
* ------
* This result comes from constructing two vectors:
*    r1 = v1 * t, the position of the primary star
*    r2 = sep_vec + v2 * t, the position of the secondary star
* The difference between these vectors is
*    d = sep_vec + (v2 - v1) * t
* If the magnitude of d is less than r2, the stars collide.
* So to find the time of collision, we solve for t when minimising d^2
* (squared because then we can use the quadratic formula) and just plug
* it in.
*
* Variables
* ---------
* sep_vec is the separation vector between the two stars (in km)
* v1 is the velocity of the primary star (in km/s)
* v2 is the velocity of the secondary star (in km/s)
* r2 is the radius of the secondary star (in rsun)
* collide is whether the stars collide (logical)

      real*8 sep_vec(3), v1(3), v2(3), v_dif(3), r2, r2km
      real*8 r_dot_v_dif, v_dif_dot, t_min, d_min, d_min_vec(3)
      logical collide
      integer i

      do i = 1, 3
         v_dif(i) = v2(i) - v1(i)
      end do

      call DotProduct(sep_vec, v_dif, r_dot_v_dif)
      call DotProduct(v_dif, v_dif, v_dif_dot)
      t_min = -r_dot_v_dif / v_dif_dot

      if (t_min.lt.0) then
         collide = .false.
      else
         do i = 1, 3
            d_min_vec(i) = sep_vec(i) + v_dif(i) * t_min
         end do
         call VectorMagnitude(d_min_vec, d_min)
         r2km = r2 * rsunkm
         if (d_min.lt.r2km) then
            collide = .true.
         else
            collide = .false.
         end if
      end if 

      RETURN
      END

      SUBROUTINE RandomLogNormal(mean, sigma, result, idum1, twopi)
* This function generates a random number from a log-normal distribution
* following the Box-Muller transform method.
* http://en.wikipedia.org/wiki/Box-Muller_transform

      real*8 mean, sigma, result, twopi
      real*8 u1, u2, z0

      u1 = ran3(idum1)
      u2 = ran3(idum1)

      if (u1.le.0.d0) u1 = 1.0E-10      ! Avoid log(0)

      Z0 = SQRT(-2.0d0 * LOG(u1)) * COS(twopi * u2)
      result = EXP(mean + sigma * Z0)

      RETURN
      END