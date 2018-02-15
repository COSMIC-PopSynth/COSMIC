***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vk,snstar,
     &                r2,fallback,bkick)
      implicit none
*
* Updated JRH kick routine by PDK (see Kiel & Hurley 2009).
*
* Here theta is the \omega angle within the HTP02 paper (thus phi is phi).
*
* Produces a random SN kick.
* Evolution was added for binaries in which the kick creates
* an eccentricity of greater than unity (i.e. hyperbolic orbit).
*
* In particular, bkick is an array with 12 elements containing information
* of velocities from possibly multiple kicks. For example, if star 2 explodes
* and the binary disruptes (and star 1 has already exploded but the binary 
* survived) then the array will be full and will look like:
* bkick[1] = 1, bkick[2-4] = vx, vy, vz, bkick[5] = 2, bkick[6-8] = vx, vy, vz
* of star 2, bkick[9] = 2, bkick[10-12] = vx, vy, vz of star 1.
* Therefore, to update the Galactic or cluster stellar velocities we will need
* to add bkick[2-4] to both stars and then bkick[6-8] to star 2 and 
* bkick[10-12] to star 1. Here, star 1 is the star that was passed first to
* evolv1.f or evolv2.f.
*
* Add tphys to input params, both above and within evolv2.f
* bkick(5) was made negative in the bse_interface routine
* (originaly it was kick(6)).
* For cmc bkick array is zero, not negative.
*
* Perhaps move wdaflag to COMMON /VALUE4/ with bhflag etc?
*
      integer kw,k,l,snstar
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag,wdaflag
      real*8 tphys,m1,m2,m1n,mbi,mbf,mtilda,mdif
      real*8 ecc,sep,sepn,jorb,ecc2,bb,angle
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      real*8 omega,sino,coso,cosmu,sinmu,x_tilt,y_tilt,z_tilt
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi,alpha,beta,gamma
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha,csig,ssig
      real*8 sthetatest,cthetatest,ct,st,cp,sp,sphitest,cphitest
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 mu,cmu,vs(3),v1,v2,v1a,v1b
      real*8 Ptt,Qtt,Rtt,Stt,mx1,mx2,r2,qdist
      real*8 sigma,sigmah,bhsigmafrac,RotInvX,RotInvZ
      real*8 signs,sigc,psins,psic,cpsins,spsins,cpsic,spsic
      real*8 csigns,ssigns,csigc,ssigc
      real*8 semilatrec,cangleofdeath,angleofdeath,energy
      real*8 fallback,kickscale,bound,phi_old
* Output
      real*8 v1out,v2out,v3out,vkout
      logical output
*
      real*8 bconst,CK,opening_angle
      COMMON /VALUE4/ sigma,bhsigmafrac,bconst,CK,bhflag,opening_angle
      real*8 mxns,neta,bwind,hewind
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      real*8 bkick(16)
*      COMMON /VKICK/ bkick
      real ran3,xx
      external ran3
      integer*8 id1_pass,id2_pass
      real*8 merger
      COMMON /cmcpass/ merger,id1_pass,id2_pass
*
      output = .false. !useful for debugging...
*       write(91,49)kw,m1,m1n,m2,ecc,sep,snstar,fallback,
*    &               bhflag,sigma,mxns,id1_pass,id2_pass
      v1out = 0.d0
      v2out = 0.d0
      v3out = 0.d0
      vkout = 0.d0
* Scaling owing to ECSN.
      kickscale = 0.d0
*
*
      if(sigma.lt.0.d0)then 
         sigma = -1.d0*sigma
         kickscale = 10.d0
      endif
      sigmah = sigma
*Test: Checking if we can make customized sigma for blackholes only
      if(kw.eq.14.or.(kw.eq.13.and.(m1n.ge.mxns)))then
          sigma = sigmah*bhsigmafrac
      endif
      if(output)then

        if(kw.eq.14) write(48,49)kw,m1,m1n,m2,ecc,sep,snstar,fallback,
     &               bhflag,sigma,mxns,id1_pass,id2_pass
      endif
*      if(kw.lt.0)then
*         sigma = 20.d0
*         kw = ABS(kw)
*      endif
*      if(kw.eq.14.and.bhflag.eq.2)then
* Voss & Tauris (2003) method of limiting BH kick momentum.
* Here 3.d0 is the maximum mass of a NS - in future should be
* generalised, i.e. pass through value of Mns,max, rather than
* assume Mns,max = 3.d0 Msun.
*
* Also should implement kick mag. inv. propto. fall back fraction.
* Done, see fallback below.
*         sigma = (mxns/m1n)*sigma
*      endif
      do k = 1,3
         vs(k) = 0.d0
      enddo
      u1 = 0.d0
      u2 = 0.d0
*
      vk = 0.d0
*      bkick(1) = -1.d0
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
*         write(15,*)'kick 1:',xx,idum
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
*         xx = RAN3(idum) !randomise initial orientation... dont require...
*         if(xx.gt.0.5d0) salpha = -1.d0*salpha
*         xx = RAN3(idum)
*         if(xx.gt.0.5d0) calpha = -1.d0*calpha
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
*         write(15,*)'kick 2:',u1,idum
         u2 = RAN3(idum)
         if(u1.gt.0.9999d0) u1 = 0.9999d0
         if(u2.gt.1.d0) u2 = 1.d0
*         write(15,*)'kick 3:',u2,idum
* Generate two velocities from polar coordinates S & THETA.
         s = -2.d0*LOG(1.d0 - u1)
         s = sigma*SQRT(s)
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
      vk = SQRT(vk2)
* Limit BH kick with fallback mass fraction.
*      if(kw.eq.14)then
*Limit BH kick with fallback only if wanted
*      write(20,*)'BH FORM', m1,vk,fallback,kw
      if(kickscale.gt.0.d0)then
         vk = vk/kickscale
         vk2 = vk2/kickscale/kickscale
      endif
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

      sigma = sigmah
* CLR - Allow for a restricted opening angle for SN kicks
*       Only relevant for binaries, obviously
      bound = SIN((90.d0 - opening_angle)*pi/180.d0)
      sphi = (1.d0-bound)*ran3(idum) + bound
      phi = ASIN(sphi)
      cphi = COS(phi)
* CLR - if the orbit has already been kicked, then any polar kick
*       needs to be tilted as well (since L_hat and S_hat are no longer
*       aligned).  Here we take the random kick from above and rotate it
*       about the X axis by the angle mu from the last SN kick
      if(bkick(1).gt.0.d0.and.bkick(5).le.0.d0)then
        cosmu = bkick(13)
        sinmu = sqrt(1 - bkick(13)*bkick(13))
        z_tilt=cosmu*sphi + cphi*sinmu*stheta
        phi = ASIN(z_tilt)
        sphi = z_tilt
        cphi = COS(phi)
*       write(15,*) vk,theta,phi,bkick(13)
      endif
      theta = twopi*ran3(idum)
      stheta = SIN(theta)
      ctheta = COS(theta)

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
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      if(bkick(1).le.0.d0)then
        bkick(13) = cmu
      elseif(bkick(5).le.0.d0)then
        bkick(14) = cmu
      endif
      mu = ACOS(cmu)
* Determine if orbit becomes hyperbolic.
 90   continue
      mx1 = vk*m1n/(m1n+m2)
      mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
      vs(1) = mx1*ctheta*cphi + mx2*salpha
      vs(2) = mx1*stheta*cphi + mx2*calpha
      vs(3) = mx1*sphi
*
* Introduce random orientation of binary system to the Galaxy.
* Not completed here but within binkin/cmc...
*
*      alpha = twopi*ran3(idum)
*      gamma = twopi*ran3(idum)
*      beta = pi*ran3(idum)
*
      if(ecc.gt.1.d0)then
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
*         if(kw.eq.14) write(*,*)'0.5: ',angleofdeath,ctheta
* Find which hemisphere (in x) the NS/companion originated in...
         psins = angleofdeath 
         if((stheta*cphi - calpha).gt.0.d0) psins = -psins
         psic = psins
* Accounting for rotation in x-y plane due to SN event.
         cpsins = COS(signs-psins)
         spsins = SIN(signs-psins)
         cpsic = COS(sigc-psic)
         spsic = SIN(sigc-psic)
*         if(kw.eq.14) write(*,*)'1: ',semilatrec,r,cangleofdeath,
*     &        angleofdeath,sigc,psic,cpsic,spsic
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
         if(bkick(1).le.0.d0)then
            bkick(1) = float(snstar)
*     
            bkick(2) = ((m1n/mbf)*vk*ctheta*cphi + 
     &           (mdif*m2)/(mbi*mbf)*vr*salpha) -
     &           v1*(spsins)*RotInvX*(m2/mbf)
*     
            bkick(3) = (m1n/mbf)*vk*stheta*cphi + 
     &           (mdif*m2)/(mbi*mbf)*vr*calpha - 
     &           v1*(cpsins)*(m2/mbf)
*     
            bkick(4) = (m1n/mbf)*vk*sphi
*     
            bkick(5) = float(snstar)
*     
            bkick(6) = ((mdif*m2)/(mbi*mbf)*vr*salpha + 
     &           (m1n/mbf)*vk*ctheta*cphi) + 
     &           v1*(spsic)*(m1n/mbf)*RotInvX
*     
            bkick(7) = (mdif*m2)/(mbi*mbf)*vr*calpha + 
     &           (m1n/mbf)*vk*stheta*cphi + 
     &           v1*(cpsic)*(m1n/mbf)
*
            bkick(8) = (m1n/mbf)*vk*sphi
*
            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  bkick(2) = vs(1)
                  bkick(3) = vs(2)
                  bkick(4) = vs(3)
                  bkick(6) = 0.d0
                  bkick(7) = 0.d0
                  bkick(8) = 0.d0
                  m2 = -1.d0*m2
               endif
            endif
            v1out = bkick(2)
            v2out = bkick(3)
            v3out = bkick(4)
         elseif(bkick(1).gt.0.d0.and.bkick(5).le.0.d0)then
            bkick(5) = float(snstar)
*     
            bkick(6) = ((m1n/mbf)*vk*ctheta*cphi + 
     &           (mdif*m2)/(mbi*mbf)*vr*salpha) -
     &           v1*(spsins)*RotInvX*(m2/mbf)
*     
            bkick(7) = (m1n/mbf)*vk*stheta*cphi + 
     &           (mdif*m2)/(mbi*mbf)*vr*calpha - 
     &           v1*(cpsins)*(m2/mbf)
*
            bkick(8) = (m1n/mbf)*vk*sphi
*     
            bkick(9) = float(snstar)
*     
            bkick(10) = ((mdif*m2)/(mbi*mbf)*vr*salpha + 
     &           (m1n/mbf)*vk*ctheta*cphi) + 
     &           v1*(spsic)*(m1n/mbf)*RotInvX
*     
            bkick(11) = (mdif*m2)/(mbi*mbf)*vr*calpha + 
     &           (m1n/mbf)*vk*stheta*cphi + 
     &           v1*(cpsic)*(m1n/mbf)
*     
            bkick(12) = m1n/mbf*vk*sphi
*
            if(psins.lt.0.d0)then
               if(r2.gt.sepn*(ecc - 1.d0))then
                  bkick(6) = vs(1)
                  bkick(7) = vs(2)
                  bkick(8) = vs(3)
                  bkick(10) = 0.d0
                  bkick(11) = 0.d0
                  bkick(12) = 0.d0
                  m2 = -1.d0*m2
               endif
            endif
            v1out = bkick(6)
            v2out = bkick(7)
            v3out = bkick(8)
*     
         endif
         ecc = MIN(ecc,99.99d0)
      endif
* Calculate the components of the velocity of the new centre-of-mass.
* 90   mx1 = vk*m1n/(m1n+m2)
* 4    continue
      if(ecc.lt.1.d0)then
*         if(ecc.eq.1.d0.or.ecc.lt.0.d0) m2 = -1.d0 * m2
* 1st time with kick.
         if(bkick(1).le.0.d0)then
            bkick(1) = float(snstar)
            bkick(2) = vs(1)
            bkick(3) = vs(2)
            bkick(4) = vs(3)
            v1out = bkick(2)
            v2out = bkick(3)
            v3out = bkick(4)
* 2nd time with kick.
         elseif(bkick(5).le.0.d0)then
            bkick(5) = float(snstar)
            bkick(6) = vs(1)
            bkick(7) = vs(2)
            bkick(8) = vs(3)
            v1out = bkick(6)
            v2out = bkick(7)
            v3out = bkick(8)
* 2nd time with kick if already disrupted.
         elseif(bkick(5).gt.0.d0)then
            bkick(9) = float(snstar)
            bkick(10) = vs(1)
            bkick(11) = vs(2)
            bkick(12) = vs(3)
            v1out = bkick(10)
            v2out = bkick(11)
            v3out = bkick(12)
         endif
      endif
*
      if(ecc.eq.1.d0)then
* 1st time with kick.
         if(bkick(1).le.0.d0)then
            bkick(1) = float(snstar)
            bkick(2) = vs(1)
            bkick(3) = vs(2)
            bkick(4) = vs(3)
            bkick(5) = float(snstar)
            bkick(6) = -vs(1)
            bkick(7) = -vs(2)
            bkick(8) = -vs(3)
            v1out = bkick(2)
            v2out = bkick(3)
            v3out = bkick(4)
* 2nd time with kick.
         elseif(bkick(1).gt.0.d0.and.bkick(5).le.0.d0)then
            bkick(5) = float(snstar)
            bkick(6) = vs(1)
            bkick(7) = vs(2)
            bkick(8) = vs(3)
            bkick(9) = float(snstar)
            bkick(10) = -vs(1)
            bkick(11) = -vs(2)
            bkick(12) = -vs(3)
            v1out = bkick(6)
            v2out = bkick(7)
            v3out = bkick(8)
         endif
      endif
* Randomly rotate system
      CALL randomness3(idum,bkick(2),bkick(3),bkick(4),bkick(6),
     &            bkick(7),
     &            bkick(8),bkick(10),bkick(11),bkick(12))
*
      if(ecc.gt.99.9d0) ecc = 99.9d0
      if(output)then
         if(sep.le.0.d0.or.ecc.ge.1.d0)then
            vkout = sqrt(v1out*v1out+v2out*v2out+v3out*v3out)
            write(44,43)kw,m1,m1n,sigma,vk,v1out,v2out,v3out,vkout,
     &                  (bkick(l),l=1,12),id1_pass,id2_pass
         else
            vkout = sqrt(v1out*v1out+v2out*v2out+v3out*v3out)
            write(45,47)kw,m1,m1n,sigma,vk,v1out,v2out,v3out,
     &                  vkout,id1_pass,id2_pass
         endif
      endif
 43   FORMAT(i3,1p,8e12.4,1x,12e12.4,1x,i10,i10)
 47   FORMAT(i3,1p,8e12.4,1x,i10,i10)
 49   FORMAT(i3,1p,5e12.4,1x,0p,i3,f12.4,1x,i3,f12.4,f12.4,1x,i10,i10)
*
      RETURN
      END
***
      SUBROUTINE randomness3(idum,vx1,vy1,vz1,vx2,vy2,vz2,
     &                            vx3,vy3,vz3)
*
      implicit none
*
      INTEGER idum
      real ran3
      external ran3
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      REAL*8 vx1,vy1,vz1,alpha,gamma,beta,pi,twopi,vx2,vy2,vz2
      REAL*8 cg,sg,ca,sa,cb,sb,vx1s,vy1s,vz1s,vx2s,vy2s,vz2s
      REAL*8 vx3,vy3,vz3,vx3s,vy3s,vz3s
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
      vx3s = vx3
      vy3s = vy3
      vz3s = vz3
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
**
      vx3 = vx3s*(cb*ca)
      vx3 = vx3 + vy3s*(cb*sa) 
      vx3 = vx3 - vz3s*sb
*     
      vy3 = vx3s*(sg*sb*ca - cg*sa)
      vy3 = vy3 + vy3s*(sg*sb*sa + cg*ca)
      vy3 = vy3 + vz3s*(cb*sg)
*     
      vz3 = vx3s*(cg*sb*ca + sg*sa)
      vz3 = vz3 + vy3s*(cg*sb*sa - sg*ca)
      vz3 = vz3 + vz3s*(cb*cg)
*
*
      RETURN
      END
*
