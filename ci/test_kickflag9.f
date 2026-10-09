      program test_kickflag9
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      include 'const_bse.h'
      integer n,i,j,flag,slot,count
      parameter(n=30000)
      real*8 draws(n),vk,info(2,19),a,e,p,expected,threshold
      real*8 fractions(3),factor,raw,relative(3),recoil(3)
      real*8 g,rsun,mass_pre,mass_post,delta_mass
      logical disrupted
      parameter(g=1.3271244d11,rsun=6.96d5)
      real*8 mixture_cdf
      external mixture_cdf
      data fractions/0.2d0,0.5d0,1.d0/

* Check the independently normalized mixture CDF on a log-speed grid.
      idum1 = -20261009
      do i=1,n
         call single_kick(13,1.4d0,3,1.d0,-100.d0,1,
     &                    vk,info)
         if(.not.ieee_is_finite(vk)) error stop 'nonfinite NS kick'
         if(vk.le.0.05d0.or.vk.ge.1000.d0) error stop 'NS kick bounds'
         draws(i) = vk
      enddo
      threshold = sqrt(log(2000.d0)/(2.d0*n))
      do j=0,60
         a = exp(log(0.05d0)+j*log(1000.d0/0.05d0)/60.d0)
         count = 0
         do i=1,n
            if(draws(i).le.a) count = count+1
         enddo
         p = dble(count)/n
         expected = mixture_cdf(a)
         if(abs(p-expected).ge.threshold) error stop 'mixture CDF'
      enddo

* Same-seed BH draws must use the same base law with scaling disabled.
      idum1 = -20261009
      do i=1,n
         call single_kick(14,7.d0,3,1.d0,-100.d0,1,vk,info)
         if(vk.ne.draws(i)) error stop 'NS/BH base law mismatch'
      enddo

* Independent scale-factor checks: fallback is fixed at 0.4 here.
      idum1 = -314159
      call single_kick(14,7.d0,3,1.d0,-100.d0,1,raw,info)
      do flag=0,4
         factor = 1.d0
         if(flag.eq.0) factor = 0.d0
         if(flag.eq.1.or.flag.eq.4) factor = 0.6d0
         if(flag.eq.2) factor = 3.d0/7.d0
         do j=1,3
            idum1 = -314159
            call single_kick(14,7.d0,flag,fractions(j),
     &                       -100.d0,1,vk,info)
            expected = raw*fractions(j)*factor
            if(abs(vk-expected).gt.1.d-12) error stop 'BH scaling'
         enddo
      enddo

* Supplied magnitudes override sampling, including the second SN slot.
      do slot=1,2
         idum1 = -271828
         call single_kick(13,1.4d0,3,1.d0,123.4d0,slot,vk,info)
         if(vk.ne.123.4d0) error stop 'supplied NS kick'
         call single_kick(14,7.d0,3,0.2d0,123.4d0,slot,vk,info)
         if(vk.ne.123.4d0) error stop 'supplied BH kick'
         call single_kick(14,7.d0,4,0.2d0,123.4d0,slot,vk,info)
         if(abs(vk-74.04d0).gt.1.d-12) error stop 'supplied fallback'
      enddo

* Circular zero-kick limit: a'=3a and e=DeltaM/M_post=2/3.
      call configure_kick(3,1.d0)
      natal_kick_array(1,1:4) = 0.d0
      info = 0.d0
      idum1 = -161803
      a = 100.d0
      e = 0.d0
      disrupted = .false.
      call orbit_kick(a,e,vk,info,disrupted)
      if(disrupted) error stop 'Blaauw limit disrupted'
      if(abs(a-300.d0).gt.1.d-9) error stop 'Blaauw separation'
      if(abs(e-2.d0/3.d0).gt.1.d-12) error stop 'Blaauw eccentricity'
      mass_pre = 4.d0
      mass_post = 2.4d0
      delta_mass = 1.6d0
      relative = 0.d0
      relative(2) = sqrt(g*mass_pre/(100.d0*rsun))
      recoil = -delta_mass*relative/(mass_pre*mass_post)
      if(maxval(abs(info(1,7:9)-recoil)).gt.1.d-10)
     &   error stop 'Blaauw systemic recoil'

* A high supplied kick must unbind the same initially circular orbit.
      call configure_kick(3,1.d0)
      natal_kick_array(1,1:4) =
     &   (/250.d0,30.d0,60.d0,0.d0/)
      info = 0.d0
      a = 100.d0
      e = 0.d0
      disrupted = .false.
      call orbit_kick(a,e,vk,info,disrupted)
      if(.not.disrupted.or.e.le.1.d0.or.a.ge.0.d0)
     &   error stop 'high-kick disruption'
      if(.not.all(ieee_is_finite(info))) error stop 'nonfinite orbit'

* Exercise both SN slots across bound and disrupted binary histories.
      call binary_histories(1000)

      write(*,*) 'kickflag=9 native checks passed: ',n,
     &           ' NS and ',n,' paired BH draws'
      end

      subroutine binary_histories(n)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      include 'const_bse.h'
      integer n,i,bound_first
      real*8 info(2,19),a,e,vk,jorb,fb
      logical disrupted
      bound_first = 0
      do i=1,n
         call configure_kick(3,1.d0)
         idum1 = -i
         info = 0.d0
         a = 100.d0
         e = 0.4d0
         jorb = 0.d0
         fb = 0.d0
         disrupted = .false.
         call kick(13,3.d0,2.d0,1.4d0,3.d0,e,a,jorb,vk,1,
     &             0.d0,fb,265.d0,info,disrupted,0.d0)
         if(.not.disrupted) bound_first = bound_first+1
         if(.not.all(ieee_is_finite(info))) error stop 'first SN finite'
         call kick(13,3.d0,2.d0,1.4d0,1.4d0,e,a,jorb,vk,2,
     &             0.d0,fb,265.d0,info,disrupted,1.d0)
         if(.not.all(ieee_is_finite(info)))
     &   error stop 'second SN finite'
         if(.not.disrupted.and.(a.le.0.d0.or.e.ge.1.d0))
     &      error stop 'bound binary state'
      enddo
      if(bound_first.eq.0.or.bound_first.eq.n)
     &   error stop 'missing binary outcome'
      end

      subroutine configure_kick(flag,fraction)
      implicit none
      include 'const_bse.h'
      integer flag
      real*8 fraction
      kickflag = 9
      bhflag = flag
      bhsigmafrac = fraction
      sigma = 265.d0
      mxns = 3.d0
      polar_kick_angle = 90.d0
      using_cmc = 0
      natal_kick_array = -100.d0
      natal_kick_array(:,5) = 0.d0
      end

      subroutine single_kick(kw,mnew,flag,fraction,supplied,slot,
     &                       vk,info)
      implicit none
      include 'const_bse.h'
      integer kw,flag,slot
      real*8 mnew,fraction,supplied,vk,info(2,19),a,e,jorb,fb
      logical disrupted
      call configure_kick(flag,fraction)
      natal_kick_array(slot,1) = supplied
      info = 0.d0
      a = 0.d0
      e = -1.d0
      jorb = 0.d0
      fb = 0.4d0
      disrupted = .false.
      call kick(kw,8.d0,2.d0,mnew,0.d0,e,a,jorb,vk,slot,
     &          0.d0,fb,265.d0,info,disrupted,0.d0)
      end

      subroutine orbit_kick(a,e,vk,info,disrupted)
      implicit none
      real*8 a,e,vk,info(2,19),jorb,fb
      logical disrupted
      jorb = 0.d0
      fb = 0.d0
      call kick(13,3.d0,2.d0,1.4d0,1.d0,e,a,jorb,vk,1,
     &          0.d0,fb,265.d0,info,disrupted,0.d0)
      end

      real*8 function mixture_cdf(v)
      implicit none
      real*8 v,normal_cdf,lo,hi,x
      external normal_cdf
      lo = normal_cdf((log(0.05d0)-1.87d0)/0.55d0)
      hi = normal_cdf((log(1000.d0)-1.87d0)/0.55d0)
      x = normal_cdf((log(v)-1.87d0)/0.55d0)
      mixture_cdf = 0.126d0*(x-lo)/(hi-lo)
      lo = normal_cdf((log(0.05d0)-5.62d0)/0.71d0)
      hi = normal_cdf((log(1000.d0)-5.62d0)/0.71d0)
      x = normal_cdf((log(v)-5.62d0)/0.71d0)
      mixture_cdf = mixture_cdf+0.874d0*(x-lo)/(hi-lo)
      end

      real*8 function normal_cdf(x)
      implicit none
      real*8 x
      normal_cdf = 0.5d0*(1.d0+erf(x/sqrt(2.d0)))
      end
