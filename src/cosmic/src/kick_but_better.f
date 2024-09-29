      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,snstar,
     &                r2,fallback,sigmahold,kick_info,disrupt,bkick)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* TODO: Add documentation of new function
* TODO
* TODO
* TODO                  TODO
* TODO
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
* kick_info[i,17]: random seed at the start of call to kick.f
*
* For cmc kick_info array is zero, not negative.
      integer kw,k,snstar,sn,safety

      real*8 m1,m2,m1n,mbi,mbf,mdif
      real*8 ecc,sep
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mean_anom,ecc_anom,dif,der,del,r
      real*8 u1,u2,vk,vk2,v(4),s,sigmah
      real*8 theta,phi,sin_phi,cos_phi,sin_theta,cos_theta
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
      real*8 bkick(20),r2,jorb
      real*8 ecc_prev, a_prev, mtot_prev, mred_prev
      real*8 natal_kick(3), sep_vec(3), v_rel_vec(3)
      real*8 a_prev_2, a_prev_3, cos_ecc_anom, sin_ecc_anom
      real*8 sqrt_1me2, sep_prev, prefactor, omega
      real*8 h_prev(3), LRL_vec(3)
      integer i
* Output
      logical output,disrupt
*
      real*8 kick_info(2,17)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      safety = 0

* Setup previous values from before the SN
      ecc_prev = ecc
      a_prev = sep

      mtot_prev = m1 + m2
      mred_prev = m1 * m2 / mtot_prev

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
* TODO: Is this not redundant? snstar is already passed as an argument
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
      kick_info(sn,17) = idum1

* set the SNstar of the exploding object in the kick_info array
      kick_info(sn,1) = snstar

* if the system was disrupted from the first supernova, marked as disrupted
      if(kick_info(1,2).eq.1) kick_info(2,2)=1

* sigma is negative for ECSN
      if((sigma.lt.0.d0).and.(kickflag.eq.0))then
         sigma = -1.d0*sigma
* for kick prescriptions other than default, revert to original sigma
      elseif((sigma.lt.0.d0).and.(kickflag.lt.0))then
         sigma = sigmahold
      endif
      sigmah = sigma

* scale down BH kicks if bhsigmafrac is specified
      if(kickflag.eq.0)then
         if(kw.eq.14.or.(kw.eq.13.and.(m1n.ge.mxns)))then
            sigma = sigmah*bhsigmafrac
         endif
      endif

*
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

          if(kickflag.eq.0)then
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
          elseif(kickflag.eq.-1)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 1
             vk = vk * ((m1-m1n)/mean_mej) * (mean_mns/m1n)
             vk2 = vk*vk
          elseif(kickflag.eq.-2)then
* Use kick scaling from Giacobbo & Mapelli 2020, Eq. 2
             vk = vk * ((m1-m1n)/mean_mej)
             vk2 = vk*vk
          elseif(kickflag.eq.-3)then
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

* ----------------------------------------------------------------------
* ----- Done scaling kick magnitudes at this point, now for angles -----
* ----------------------------------------------------------------------

* Before we randomly draw a phi and theta for the natal kick,
* see if a pre-supplied set of phi/theta is passed
      if((natal_kick_array(snstar,2).ge.(-90.d0)).and.
     &       (natal_kick_array(snstar,2).le.(90.d0)))then
          phi = natal_kick_array(snstar,2)*pi/180.d0
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

      if((natal_kick_array(snstar,3).ge.(0.d0)).and.
     &       (natal_kick_array(snstar,3).le.(360.d0)))then
          theta = natal_kick_array(snstar,3)*pi/180.d0
* per supplied kick value we mimic a call to random number generator
           xx = RAN3(idum1)
      else
          theta = twopi*ran3(idum1)
      endif
      sin_theta = SIN(theta)
      cos_theta = COS(theta)

*     save theta and phi (exploding star frame) in the kick_info and
*     natal_kick_array
      kick_info(sn,4) = phi*180/pi
      kick_info(sn,5) = theta*180/pi
      if(using_cmc.eq.0)then
          natal_kick_array(snstar,2) = phi*180/pi
          natal_kick_array(snstar,3) = theta*180/pi
      endif

* create a vector for the natal kick
      natal_kick(1) = vk * cos_phi * cos_theta
      natal_kick(2) = vk * cos_phi * sin_theta
      natal_kick(3) = vk * sin_phi

* Check if the system is already disrupted
      if(a_prev.le.0.d0.and.ecc_prev.lt.0.d0)then
* if the system already disrupted, only apply kick to the current star
        if(snstar.eq.1)then
            kick_info(sn,7) = natal_kick(1)
            kick_info(sn,8) = natal_kick(2)
            kick_info(sn,9) = natal_kick(3)
         elseif(snstar.eq.2)then
            kick_info(sn,11) = natal_kick(1)
            kick_info(sn,12) = natal_kick(2)
            kick_info(sn,13) = natal_kick(3)
         endif
         GOTO 73
      else
* Find the initial separation by randomly choosing a mean anomaly.
* check is user supplied mean anomaly
         if((natal_kick_array(snstar,4).ge.(0.d0)).and.
     &       (natal_kick_array(snstar,4).le.(360.d0)))then

             mean_anom = natal_kick_array(snstar,4) * pi / 180.d0
* per supplied kick value we mimic a call to random number generator
             xx = RAN3(idum1)

* Solve Kepler's equation for the eccentric anomaly from mean anomaly
* https://en.wikipedia.org/wiki/Eccentric_anomaly
             ecc_anom = mean_anom

             if(mean_anom.eq.0.d0) goto 3

  4          dif = ecc_anom - ecc * SIN(ecc_anom) - mean_anom
             if(ABS(dif / mean_anom).le.1.0d-04) goto 3
             der = 1.d0 - ecc * COS(ecc_anom)
             del = dif/der
             ecc_anom = ecc_anom - del
             goto 4

         endif

* TODO: this code seems redundant, why can't we just set the mean anomaly?
         xx = RAN3(idum1)
         mean_anom = xx*twopi
         ecc_anom = mean_anom
 2       dif = ecc_anom - ecc*SIN(ecc_anom) - mean_anom
         if(ABS(dif / mean_anom).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(ecc_anom)
         del = dif/der
         ecc_anom = ecc_anom - del
         goto 2
 3       continue
      endif

* Find the initial separation and relative velocity vectors
      a_prev_2 = a_prev * a_prev
      a_prev_3 = a_prev_2 * a_prev
      cos_ecc_anom = COS(ecc_anom)
      sin_ecc_anom = SIN(ecc_anom)

      omega = SQRT(gmrkm * mtot_prev / a_prev_3)

      sqrt_1me2 = SQRT(1.d0 - ecc_prev * ecc_prev)

      sep_vec(1) = a_prev * (cos_ecc_anom - ecc_prev)
      sep_vec(2) = a_prev * sin_ecc_anom * sqrt_1me2
      sep_vec(3) = 0.d0

      call VectorMagnitude(sep_vec, sep_prev)

      prefactor = omega * a_prev_2 / sep_prev

      v_rel_vec(1) = -prefactor * sin_ecc_anom
      v_rel_vec(2) = prefactor * cos_ecc_anom * sqrt_1me2
      v_rel_vec(3) = 0.d0

* Calculate the specific angular momentum vector
      call CrossProduct(sep_vec, v_rel_vec, h_prev)

* Calculate the Laplace-Runge-Lenz vector
      call CrossProduct(v_rel_vec, h_prev, LRL_vec)
      DO i = 1, 3
         LRL_vec(i) = LRL_vec(i) / (gmrkm * mtot_prev) - sep_vec(i) / sep_prev
      END DO

      call VectorMagnitude(LRL_vec, ecc)

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
      
      RETURN
      END
      

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
