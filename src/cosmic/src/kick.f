      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,sn,
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
* kick_info is a (2,18) array that tracks information about the supernova
* kicks. This allows us to track the total change to the systemic
* velocity and the total change in the orbital plane tilt after both
* supernovae, as well as reproduce systems.
* The first row contains information about the first supernova that
* occurs, the second row the second supernova.
* Note that some values the second row will take into account the
* effect of the first SN (e.g., kick_info[2,10] is the total systemic
* velocity after both supernovae).
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
* kick_info[i,15]: thetaE TODO
* kick_info[i,16]: phiE TODO
* kick_info[i,17]: psiE TODO
* kick_info[i,18]: random seed at the start of call to kick.f
*
* For cmc kick_info array is zero, not negative.
      integer kw,k,sn,safety

      real*8 m1,m2,m1n
      real*8 ecc,ecc_2,sep
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mean_anom,ecc_anom,dif,der,del,r
      real*8 u1,u2,vk,vk2,v(4),s,sigmah
      real*8 theta,phi,sin_phi,cos_phi,sin_theta,cos_theta
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
      real*8 bkick(20),r2,jorb
      real*8 ecc_prev,a_prev,mtot,mtot_prev
      real*8 natal_kick(3), sep_vec(3), v_rel(3), v_rel_prev(3)
      real*8 a_prev_2, a_prev_3, cos_ecc_anom, sin_ecc_anom
      real*8 sqrt_1me2, sep_prev, prefactor, omega
      real*8 h_prev(3),h(3), h_hat(3), h_mag
      real*8 LRL_prev(3), LRL(3), e_hat(3), e_mag
      real*8 v_cm(3), v_sn(3), v_comp(3), v_inf_vec(3), v_inf
      real*8 v_sn_rot(3), v_comp_rot(3), v_cm_rot(3)
      real*8 h_cross_e_hat(3)
      real*8 thetaE, phiE, psiE
      integer i
* Output
      logical output,disrupt
*
      real*8 kick_info(2,18)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      safety = 0

* Setup previous values from before the SN
      ecc_prev = ecc
      a_prev = sep
      mtot_prev = m1 + m2

* ----------------------------------------------------------------------
* -------------- Initialise variables and constants --------------------
* ----------------------------------------------------------------------

* Set up empty arrays and constants
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


* ----------------------------------------------------------------------
* ----- Draw and scale a natal kick magnitude based on input dists -----
* ----------------------------------------------------------------------

* Before we draw the kick from the maxwellian and then scale it
* as desired, let us see if a pre-supplied natal kick maganitude
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

* Check if the system is already disrupted
      if(a_prev.le.0.d0.and.ecc_prev.lt.0.d0)then
* if the system already disrupted, only apply kick to the current star
         disrupt = .true.
* Set that it is disrupted in the kick_info array
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
         GOTO 73
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

      if(mean_anom.eq.0.d0) goto 3

  4   dif = ecc_anom - ecc * SIN(ecc_anom) - mean_anom
      if(ABS(dif / mean_anom).le.1.0d-04) goto 3
      der = 1.d0 - ecc * COS(ecc_anom)
      del = dif/der
      ecc_anom = ecc_anom - del
      goto 4

 3    continue

* ----------------------------------------------------------------------
* ------ Calculate whether system disrupts and CM velocity change ------
* ----------------------------------------------------------------------

* Some helper variables for calculations below
      a_prev_2 = a_prev * a_prev
      a_prev_3 = a_prev_2 * a_prev
      cos_ecc_anom = COS(ecc_anom)
      sin_ecc_anom = SIN(ecc_anom)
      sqrt_1me2 = SQRT(1.d0 - ecc_prev * ecc_prev)
      mtot = m1n + m2
      ecc_2 = ecc * ecc

* Orbital frequency pre-SN
      omega = SQRT(gmrkm * mtot_prev / a_prev_3)

* Separation vector before the supernova
      sep_vec(1) = a_prev * (cos_ecc_anom - ecc_prev)
      sep_vec(2) = a_prev * sin_ecc_anom * sqrt_1me2
      sep_vec(3) = 0.d0
      call VectorMagnitude(sep_vec, sep_prev)

* Relative velocity vector before the supernova
      prefactor = omega * a_prev_2 / sep_prev
      v_rel_prev(1) = -prefactor * sin_ecc_anom
      v_rel_prev(2) = prefactor * cos_ecc_anom * sqrt_1me2
      v_rel_prev(3) = 0.d0

* Specific angular momentum vector pre-SN
      call CrossProduct(sep_vec, v_rel_prev, h_prev)

* Laplace-Runge-Lenz vector pre-SN
      call CrossProduct(v_rel_prev, h_prev, LRL_prev)
      DO i = 1, 3
         LRL_prev(i) = LRL_prev(i) / (gmrkm * mtot_prev)
     &               - sep_vec(i) / sep_prev
      END DO

* Calculate the new systemic velocity of the center of mass
      do i = 1, 3
         v_cm(i) = (-m2 * (m1 - m1n) / mtot_prev / mtot)
     &           * v_rel_prev(i)
     &           + (m1n / mtot * natal_kick(i))
      end do

* New vectors after SN
      v_rel(1) = v_rel_prev(1) + natal_kick(1)
      v_rel(2) = v_rel_prev(2) + natal_kick(2)
      v_rel(3) = v_rel_prev(3) + natal_kick(3)

      call CrossProduct(sep_vec, v_rel, h)
      call CrossProduct(v_rel, h, LRL)
      DO i = 1, 3
         LRL(i) = LRL(i) / (gmrkm * mtot_prev) - sep_vec(i) / sep_prev
      END DO

* Get the Euler angles from previous kick for the rotation matrix
      thetaE = kick_info(1,15) * pi / 180.d0
      phiE = kick_info(1,16) * pi / 180.d0
      psiE = kick_info(1,17) * pi / 180.d0

* Get the new eccentricity
      call VectorMagnitude(LRL, ecc)

* ----------------------------------------------------------------------
* -------- Split based on whether this kick disrupts the system --------
* ----------------------------------------------------------------------

      if(ecc.gt.1.d0)then
* System is now disrupted
         disrupt = .true.
* Set that it is disrupted in the kick_info array
         kick_info(sn,2) = 1
         call VectorHat(LRL, e_mag, e_hat)
         call VectorMagnitude(h, h_mag)
         call VectorHat(h, h_mag, h_hat)
         call CrossProduct(h_hat, e_hat, h_cross_e_hat)

* Velocity at infinity
         v_inf = gmrkm * mtot / h_mag * sqrt(ecc_2 - 1.d0)
         do i = 1, 3
            v_inf_vec(i) = v_inf * ((-1.d0 * e_hat(i) / ecc)
     &          + SQRT(1 - 1.d0 / ecc_2) * h_cross_e_hat(i))
         end do

* Velocity of the star going supernova post-SN
         do i = 1, 3
            v_sn(i) = (m2 / mtot) * v_inf_vec(i) + v_CM(i)
         end do

* Velocity of the companion star post-SN
         do i = 1, 3
            v_comp(i) = -(m1n / mtot) * v_inf_vec(i) + v_CM(i)
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

         call AngleBetweenVectors(h, h_prev, thetaE)
         phiE = ran3(idum1) * twopi
         psiE = ran3(idum1) * twopi

         ecc = -1.d0
         sep = -1.d0
      else
* The system is still bound
* Record the mean anomaly in the arrays
         kick_info(sn,6) = mean_anom * 180 / pi
         if (using_cmc.eq.0) then
            natal_kick_array(sn,4) = mean_anom * 180 / pi
         endif

* Update the total orbital angular momentum
         jorb = m1n * m2 / mtot * h_mag

         if (sn.eq.2) then
            call ChangeBasis(v_cm, thetaE, phiE, psiE, v_cm_rot)
         else
            v_cm_rot = v_cm
         endif

         call AngleBetweenVectors(h, h_prev, thetaE)

      endif

* Set Euler angles in the kick_info array
      kick_info(sn,15) = thetaE * 180 / pi
      kick_info(sn,16) = phiE * 180 / pi
      kick_info(sn,17) = psiE * 180 / pi

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
