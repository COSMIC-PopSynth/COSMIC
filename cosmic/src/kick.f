***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,snstar,
     &                r2,fallback,sigmahold,kick_info,disrupt)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
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
      real*8 vs(3),v1,v2
      real*8 mx1,mx2,r2
      real*8 sigmah,RotInvX
      real*8 signs,sigc,psins,psic,cpsins,spsins,cpsic,spsic
      real*8 csigns
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,sigmahold,bound
      real*8 mean_mns,mean_mej,alphakick,betakick
* Output
      logical output,disrupt
*
      real*8 kick_info(2,17)
      real ran3,xx
      external ran3
*
      output = .false. !useful for debugging...
      safety = 0

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
      if(kick_info(1,1).eq.0) sn=1
      if(kick_info(1,1).gt.0) sn=2

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
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then

* check is user supplied mean anomaly
         if((natal_kick_array(snstar,4).ge.(0.d0)).and.
     &       (natal_kick_array(snstar,4).le.(360.d0)))then

             em = natal_kick_array(snstar,4)*pi/180.d0
* per supplied kick value we mimic a call to random number generator
             xx = RAN3(idum1)
             goto 3

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
*
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
      natal_kick_array(snstar,1) = vk

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
      natal_kick_array(snstar,2) = phi*180/pi
      natal_kick_array(snstar,3) = theta*180/pi

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
        cmu = COS(kick_info(1,15))
        smu = SIN(kick_info(1,15))
        comega = COS(kick_info(1,16))
        somega = SIN(kick_info(1,16))

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
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(0.d0,ecc2)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors, and randomly choose an azimuth angle
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
      smu = SIN(mu)
      omega = twopi*ran3(idum1)
      comega = COS(omega)
      somega = SIN(omega)

* Write angle between initial and current orbital angular momentum vectors
* and eccentric anomaly if system is still bound
      if(sn.eq.1.and.ecc.le.1)then
        kick_info(sn,15) = mu*180/pi
        kick_info(sn,16) = omega*180/pi
        kick_info(sn,6) = em*180/pi
        natal_kick_array(snstar,4) = em*180/pi
      elseif(sn.eq.2)then
* MJZ - Here we calculate the total change in the orbital plane
*       from both SN. Note that these angles mu and omega are in
*       typical spherical coordinates rather than colateral coordinates,
*       so the rotations are slightly different than above.
*       We rotate about z-axis by omega1 then y-axis by mu1
        cmu1 = COS(kick_info(1,15))
        smu1 = SIN(kick_info(1,15))
        comega1 = COS(kick_info(1,16))
        somega1 = SIN(kick_info(1,16))

        x_tilt = cmu1*comega*comega1*smu + cmu*smu1
     &                - cmu1*smu*somega*somega1
        y_tilt = comega1*smu*somega + comega*smu*somega1
        z_tilt = cmu*cmu1 + smu*smu1*somega*somega1
     &               - comega*comega1*smu*smu1

        kick_info(sn,15) = ACOS(z_tilt)*180/pi
        kick_info(sn,16) = ATAN(y_tilt/x_tilt)*180/pi
        kick_info(sn,6) = em*180/pi
        natal_kick_array(snstar,4) = em*180/pi

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
*
            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  kick_info(sn,7) = vs(1)
                  kick_info(sn,8) = vs(2)
                  kick_info(sn,9) = vs(3)
                  kick_info(sn,11) = 0.d0
                  kick_info(sn,12) = 0.d0
                  kick_info(sn,13) = 0.d0
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
*
            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  kick_info(sn,7) = vs(1)
                  kick_info(sn,8) = vs(2)
                  kick_info(sn,9) = vs(3)
                  kick_info(sn,11) = 0.d0
                  kick_info(sn,12) = 0.d0
                  kick_info(sn,13) = 0.d0
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
* In the impossible chance that the system is exactly parabolic...
      if(ecc.eq.1.d0.and.snstar.eq.1)then
         kick_info(sn,7) = vs(1)
         kick_info(sn,8) = vs(2)
         kick_info(sn,9) = vs(3)
         kick_info(sn,11) = -vs(1)
         kick_info(sn,12) = -vs(2)
         kick_info(sn,13) = -vs(3)
      elseif(ecc.eq.1.d0.and.snstar.eq.2)then
         kick_info(sn,7) = -vs(1)
         kick_info(sn,8) = -vs(2)
         kick_info(sn,9) = -vs(3)
         kick_info(sn,11) = vs(1)
         kick_info(sn,12) = vs(2)
         kick_info(sn,13) = vs(3)
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
