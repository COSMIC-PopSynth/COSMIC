***
      SUBROUTINE assign_remnant(zpars,mc,mcbagb,mass,kidx,mt,kw,bhspin)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      common /fall/fallback
      REAL*8 fallback
      real*8 zpars(20)

      real*8 avar,bvar
      real*8 mc,mcbagb,mass,mt,mc_tot,met
      real*8 frac,kappa,sappa,alphap,polyfit
      real*8 m_proto,m_FeNi,m_fb,bhspin,mrem,mch,dMppi
      real*8 fmix, mcritnsbh, mtemp1, mtemp2
      integer kw,kidx

* Inputs
*       zpars      : Array of metallicity dependent parameters
*       mc         : CO core mass before SN
*       mcbagb     : Core mass at the base of the AGB
*       mass       : Previous epoch mass of the star
*       kidx       : Index of the star in the pisn track arrays

* Outputs
*       mt         : Remnant mass after SN
*       kw         : Stellar (remnant) type
*       bhspin     : Dimensionless spin parameter of BH remnant

* total core mass before SN (CO + He layers)
      mc_tot = mc_co(kidx) + mc_he(kidx)

* Set the Chandrasekhar mass
      mch = 1.44d0 !set here owing to AIC ECSN model.

* Check if remnant is below Chandrasekhar mass
      if(mc_co(kidx).lt.mch)then
         if(ifflag.ge.1)then
*
* Invoke WD IFMR from HPE, 1995, MNRAS, 272, 800.
*
            if(zpars(14).ge.1.0d-08)then
               mc = MIN(0.36d0+0.104d0*mass,0.58d0+0.061d0*mass)
               mc = MAX(0.54d0+0.042d0*mass,mc)
               if(mass.lt.1.d0) mc = 0.46d0
            else
               mc = MIN(0.29d0+0.178d0*mass,0.65d0+0.062d0*mass)
               mc = MAX(0.54d0+0.073d0*mass,mc)
            endif
            mc = MIN(mch,mc)
         endif
         
         mt = mc
         if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
            kw = 11
         elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then !double check what this should be. should be ecsn_mlow. Remember need to add option if ecsn = 0 (i.e. no ECSN!!!)
*
* Zero-age Carbon/Oxygen White Dwarf
*
            kw = 11
         elseif(ecsn.gt.0.d0.and.mcbagb.ge.ecsn_mlow.and.
     &          mcbagb.le.ecsn.and.mc.lt.1.08d0)then
            kw = 11
*               elseif(mcbagb.ge.1.6d0.and.mcbagb.le.2.5d0.and.
*                      mc.lt.1.08d0)then !can introduce this into code at some point.
*                  kw = 11

         else
*
* Zero-age Oxygen/Neon White Dwarf
*
            kw = 12
         endif
         mass = mt
*
      else
         if(ecsn.gt.0.d0.and.mcbagb.lt.ecsn_mlow)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
            kw = 15
         elseif(ecsn.eq.0.d0.and.mcbagb.lt.1.6d0)then
*
* Star is not massive enough to ignite C burning.
* so no remnant is left after the SN
*
            kw = 15
         else
* Beginning of supernova block
*
* Chris Belczynski Evolutionary Roads Weak PPISN
* This has to happen before the SNa, because it modifies
* the properties of the star during explosion
            if(pisn.eq.-5.and.mt.ge.45d0)then
              if(mcbagb.ge.65d0) then
                mt = 0.d0
                kw = 15
              else
* PPISN
                if(mcbagb.ge.60d0) then
                  mtemp1 = 938d0 - (14.3d0*mcbagb)
                elseif(mcbagb.ge.40d0) then
                  mtemp1 = 55.6d0
                else
                  mtemp1 = 6.0d0 + (0.83d0*mcbagb)
                endif
* Update mass
                if(mt.gt.mtemp1) then
                  mt = mtemp1
                endif
                if(mcbagb.gt.mtemp1) then
                  mcbagb = mtemp1
                endif
                if(mc.gt.mtemp1) then
                  mc = mtemp1
                endif
              endif
            endif
* Carry on with the Supernovae
*
* Use remnant mass given by Hurley+2000
            if(remnantflag.eq.0)then
               mt = 1.17d0 + 0.09d0*mc_co(kidx)
            elseif(remnantflag.eq.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407. Equation 1
*
*   TW: Belczynski+02 states that FeNi core mass comes from Woosley 1986 (Nucleosynthesis and Stellar Evolution)
*   As far I can tell, this is a linear fit to Table 6. But that table is He core mass, not CO core mass.
*   I also get a different slope, so I imagine there was a conversion to CO core mass somewhere in their fit?
               if(mc_co(kidx).lt.2.5d0)then
                  m_FeNi = 0.161767d0*mc_co(kidx) + 1.067055d0
               else
                  m_FeNi = 0.314154d0*mc_co(kidx) + 0.686088d0
               endif
               if(mc_co(kidx).le.5.d0)then
                  mt = m_FeNi
                  fallback = 0.d0
               elseif(mc_co(kidx).lt.7.6d0)then
*   TW: the fallback fraction is assumed to linearly increase from 0 to 1 between MCO=5 and MCO=7.6, hence the 2.6
                  fallback = (mc_co(kidx) - 5.d0)/2.6d0
                  mt = m_FeNi + fallback * (mt - m_FeNi)
               elseif(mc_co(kidx).gt.7.60)then
                  fallback = 1.d0
               endif
            elseif(remnantflag.eq.2)then
*
* Use NS/BH masses given by Belczynski+08. PK.
*
               ! calculate m_FeNi following Eq. 1 (fit to Timmes+1996)
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  m_FeNi = 1.38d0
               elseif(mc_co(kidx).lt.4.82d0)then
                  m_FeNi = 1.5d0
               elseif(mc_co(kidx).ge.4.82d0
     &                .and.mc_co(kidx).lt.6.31d0)then
                  m_FeNi = 2.11d0
               elseif(mc_co(kidx).ge.6.31d0
     &                .and.mc_co(kidx).lt.6.75d0)then
                  m_FeNi = 0.69*mc_co(kidx) - 2.26d0
               elseif(mc_co(kidx).ge.6.75d0)then
                  m_FeNi = 0.37*mc_co(kidx) - 0.07d0
               endif
               ! now calculate the remnant mass after fallback (Eq. 2)
               if(mc_co(kidx).le.5.d0)then
                  mt = m_FeNi
                  fallback = 0.d0
               elseif(mc_co(kidx).lt.7.6d0)then
                  fallback = (mc_co(kidx) - 5.d0) / 2.6d0
                  mt = m_FeNi + fallback * (mt - m_FeNi)
               elseif(mc_co(kidx).gt.7.60)then
                  fallback = 1.d0
               endif
            elseif(remnantflag.eq.3)then
*
* Use the "Rapid" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*              We use the updated proto-core mass from Giacobbo & Mapelli 2020
               m_proto = 1.1d0

               ! Calculate remnant mass from Eq. 16 + 17
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  mt = 1.38d0   ! ECSN fixed mass, no fallback
               elseif(mc_co(kidx).lt.2.5d0)then
                  fallback = 0.2d0 / (mt - m_proto)
                  mt = m_proto + 0.2d0
               elseif(mc_co(kidx).lt.6.d0)then
                  m_fb = 0.286d0 * mc_co(kidx) - 0.514d0
                  fallback = m_fb / (mt - m_proto)
                  mt = m_proto + m_fb
               elseif(mc_co(kidx).lt.7.d0)then
                  fallback = 1.d0
               elseif(mc_co(kidx).lt.11.d0)then
                  avar = 0.25d0 - (1.275 / (mt - m_proto))
                  bvar = 1.d0 - 11.d0*avar
                  fallback = avar*mc_co(kidx) + bvar
                  mt = m_proto + fallback*(mt - m_proto)
               elseif(mc_co(kidx).ge.11.d0)then
                  fallback = 1.d0
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            elseif(remnantflag.eq.4)then
*
* Use the "Delayed" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*              Calculate the proto-core mass following Eq. 18
               if(mc_co(kidx).lt.3.5d0)then
                  m_proto = 1.2d0
               elseif(mc_co(kidx).lt.6.d0)then
                  m_proto = 1.3d0
               elseif(mc_co(kidx).lt.11.d0)then
                  m_proto = 1.4d0
               elseif(mc_co(kidx).ge.11.d0)then
                  m_proto = 1.6d0
               endif

               ! Calculate remnant mass from Eq. 19 + 20
               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                             mcbagb.ge.ecsn_mlow)then
                  mt = 1.38d0   ! ECSN fixed mass, no fallback
               elseif(mc_co(kidx).lt.2.5d0)then
                  fallback = 0.2d0 / (mt - m_proto)
                  mt = m_proto + 0.2d0
               elseif(mc_co(kidx).lt.3.5d0)then
                  m_fb = 0.5d0 * mc_co(kidx) - 1.05d0
                  fallback = m_fb / (mt - m_proto)
                  mt = m_proto + m_fb
               elseif(mc_co(kidx).lt.11.d0)then
                  avar = 0.133d0 - (0.093d0 / (mt - m_proto))
                  bvar = 1.d0 - 11.d0*avar
                  fallback = avar*mc_co(kidx) + bvar
                  mt = m_proto + fallback*(mt - m_proto)
               elseif(mc_co(kidx).ge.11.d0)then
                  fallback = 1.d0
               endif
*              if the user requests it, limit the final remnant mass to
*              is the total **core** mass, not the total stellar mass
               if(fryer_mass_limit.eq.1)then
                  mt = min(mt, mc_tot)
               endif
            elseif(remnantflag.eq.5)then
               call assign_remnant_mandel_muller(mc, mc_tot, mt)
            elseif(remnantflag.eq.6)then
                met = 10**(LOG10(zpars(14))/0.4)
               call assign_remnant_maltsev(mc,mc_tot,met,kidx,kw,mt)
            elseif(remnantflag.eq.7)then
*
* Use the Fryer et al. 2022 SN Prescription
*
*                    For this, we just set the proto-core mass to one
               if(mc.le.3.5d0)then
                  m_proto = 1.2d0
               elseif(mc.le.6.d0)then
                  m_proto = 1.3d0
               elseif(mc.le.11.d0)then
                  m_proto = 1.4d0
               elseif(mc.gt.11.d0)then
                  m_proto = 1.6d0
               endif

               if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &              mcbagb.ge.ecsn_mlow)then
                  mt = 1.38d0   ! ECSN fixed mass, no fallback
               else
* Parameters of Fryer2022 model
                  fmix=1.0
                  mcritnsbh=5.75
* We need mt in multiple places, so temp1 will be the working mt
                  mtemp1=mt
* mtemp2 is the calculated value of the remnant mass
                  mtemp2=1.2 + (0.05*fmix) + 
     &                (0.01*((mc/fmix)**2)) +
     &                EXP(fmix*(mc-mcritnsbh))
* We don't care about mtemp2 if it's less than zero
                  if(mtemp2.lt.0.)then
                      mtemp1 = 0.
                      kw=15
* We only care about mtemp2 if it is less than the total
*   mass of the star
                  elseif(mtemp2.lt.mt)then
                      mtemp1 = mtemp2
* If mtemp2 is less, we also want to estimate the fallback fraction
                      fallback=(mtemp1-m_proto) /(mt-m_proto)
                      mt = m_proto + fallback*(mtemp1 - m_proto)
                  endif
               endif
               mc = mt
            endif
            
* Assign the BH spin based on the chosen prescription
            call assign_remnant_spin(mc, bhspin)

            ! convert from baryonic to gravitational mass
            call baryonic_to_gravitational_mass(mt, mrem)

* Determine whether a zero-age NS or BH is formed
            if(mrem.le.mxns)then
               mt = mrem
               mc = mt
               kw = 13
            else
               mt = mrem
               mc = mt
               kw = 14

* CLR - (Pulsational) Pair-Instability Supernova

* Belczynski+2016 prescription: just shrink any BH with a He core mass
* between 45 and 65 solar masses (provided the pisn flag is set at 45),
* and blow up anything between 65 and 135 solar masses.
* Cheap, but effective
               if(pisn.gt.0)then
                  if(mc_tot.ge.pisn.and.mc_tot.lt.65.d0)then
                     mt = pisn

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                     pisn_track(kidx)=6
                  elseif(mc_tot.ge.65.d0.and.mc_tot.lt.135.d0)then
                     mt = 0.d0
                     mc = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
* The Spera+Mapelli2017 prescription is a tad more sophisticated:
* complex fitting formula to Stan Woosley's PSN models.  HOWEVER, these
* were done using the ZAMS mass/core mass/remnant mass relationships for
* SEVN, not BSE.  In other words, I woud be careful using this (and in
* practice, it doesn't vary that much from Belczynski's prescription,
* since the He core masses are the same in both)
               elseif(pisn.eq.-1)then
                  frac = mc_tot/mass
                  kappa = 0.67d0*frac + 0.1d0
                  sappa = 0.5226d0*frac - 0.52974d0
                  if(mc_tot.le.32.d0)then
                     alphap = 1.0d0
                  elseif(frac.lt.0.9d0.and.mc_tot.le.37.d0)then
                     alphap = 0.2d0*(kappa-1.d0)*mc_tot +
     &                        0.2d0*(37.d0 - 32.d0*kappa)
                     pisn_track(kidx)=6
                  elseif(frac.lt.0.9d0.and.mc_tot.le.60.d0)then
                     alphap = kappa
                     pisn_track(kidx)=6
                  elseif(frac.lt.0.9d0.and.mc_tot.lt.64.d0)then
                     alphap = kappa*(-0.25d0)*mc_tot + kappa*16.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.37.d0)then
                     alphap = sappa*(mc_tot - 32.d0) + 1.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.56.d0.and.
     &                   sappa.lt.-0.034168d0)then
                     alphap = 5.d0*sappa + 1.d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.le.56.d0.and.
     &                   sappa.ge.-0.034168d0)then
                     alphap = (-0.1381d0*frac + 0.1309d0)*
     &                        (mc_tot - 56.d0) + 0.82916d0
                     pisn_track(kidx)=6
                  elseif(frac.ge.0.9d0.and.mc_tot.lt.64.d0)then
                     alphap = -0.103645d0*mc_tot + 6.63328d0
                     pisn_track(kidx)=6
                  elseif(mc_tot.ge.64.d0.and.mc_tot.lt.135.d0)then
                     alphap = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  elseif(mc_tot.ge.135.d0)then
                     alphap = 1.0d0
                  endif
                  mt = alphap*mt

                  ! convert from baryonic to gravitational mass
                  call baryonic_to_gravitational_mass(mt, mrem)
                  mt = mrem

* Fit (8th order polynomial) to Table 1 in Marchant+2018.
               elseif(pisn.eq.-2)then
                  if(mc_tot.ge.31.99d0.and.mc_tot.le.61.10d0)then
                     polyfit = -6.29429263d5
     &                      + 1.15957797d5*mc_tot
     &                      - 9.28332577d3*mc_tot**2d0
     &                      + 4.21856189d2*mc_tot**3d0
     &                      - 1.19019565d1*mc_tot**4d0
     &                      + 2.13499267d-1*mc_tot**5d0
     &                      - 2.37814255d-3*mc_tot**6d0
     &                      + 1.50408118d-5*mc_tot**7d0
     &                      - 4.13587235d-8*mc_tot**8d0
                     mt = polyfit
                     pisn_track(kidx)=6

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                  elseif(mc_tot.gt.61.10d0.and.
     &                   mc_tot.lt.124.12d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif

* Fit (8th order polynomial) to Table 5 in Woosley2019.
               elseif(pisn.eq.-3)then
                  if(mc_tot.ge.29.53d0.and.mc_tot.le.60.12d0)then
                     polyfit = -3.14610870d5
     &                      + 6.13699616d4*mc_tot
     &                      - 5.19249710d3*mc_tot**2d0
     &                      + 2.48914888d2*mc_tot**3d0
     &                      - 7.39487537d0*mc_tot**4d0
     &                      + 1.39439936d-1*mc_tot**5d0
     &                      - 1.63012111d-3*mc_tot**6d0
     &                      + 1.08052344d-5*mc_tot**7d0
     &                      - 3.11019088d-8*mc_tot**8d0
                     mt = polyfit
                     pisn_track(kidx)=6

                     ! convert from baryonic to gravitational mass
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

                  elseif(mc_tot.gt.60.12d0.and.
     &                   mc_tot.lt.135.d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
* Apply the PPISN prescription from Renzo+2022 (https://ui.adsabs.harvard.edu/abs/2022RNAAS...6...25R/abstract)
* with the adaptations from Hendriks+2023 (https://scixplorer.org/abs/2023MNRAS.526.4130H/abstract)
* This is a top-down prescription, where we subtract mass from the total core mass
               elseif(pisn.eq.-4)then
                  if(mc_co(kidx).ge.38.d0+ppi_co_shift
     &               .and.mc_co(kidx).le.114.d0)then
*       Calculate DeltaM_PPI using Eq.6 from Hendriks+2023 (equivalently Eq.2 from Renzo+2022)
                     met = 10**(LOG10(zpars(14))/0.4)
                     dMppi = (0.0006d0 * LOG10(met) + 0.0054)
     &                      * (mc_co(kidx) - ppi_co_shift - 34.8d0)**3
     &                      - 0.0013 * (mc_co(kidx)
     &                                  - ppi_co_shift - 34.8d0)**2
     &                      + ppi_extra_ml
*       Set the remnant mass equal to the total core mass minus the PPI mass loss.
*       We use core mass not total mass because envelopes are expected to be removed by the first PPI pulse (e.g. Renzo+2020b)
                     mt = mc_tot - dMppi
                     
                     call baryonic_to_gravitational_mass(mt, mrem)
                     mt = mrem

*       If the remnant mass is reduced below 10 Msun, assume a full PISN with no remnant
                     if(mt.lt.10.d0)then
                        mt = 0.0d0
                        kw = 15
                        pisn_track(kidx)=7
*       Otherwise we have a PPISN
                     else
                        pisn_track(kidx)=6
                     endif
*       For very large cores, we assume a full PISN with no remnant
                  elseif(mc_co(kidx).gt.114.d0)then
                     mt = 0.d0
                     kw = 15
                     pisn_track(kidx)=7
                  endif
               endif

               mc = mt
* Store the initial BH mass for calculating the ISCO later
               if(Mbh_initial.eq.0)then
                  Mbh_initial = mt
               endif
            endif
         endif
      endif
*
      end


      SUBROUTINE baryonic_to_gravitational_mass(mt, mrem)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      real*8 mt, mrem

      ! remnantflag 0 and 1 already calculate gravitational mass
      if(remnantflag.le.1)then
         mrem = mt
      else
         ! negative values set the absolute maximum mass loss
         if(rembar_massloss.ge.0d0)then
            ! calculate Mrem from mt using quadratic formula
            ! mt - mrem = 0.075 mrem^2 (Lattimer & Yahil 1989, Timmes+1996)
            mrem = 6.6666667d0*(SQRT(1.d0 + 0.3d0*mt) - 1.d0)

            ! limit to maximum mass loss
            mrem = MAX(mrem, mt - rembar_massloss)

         ! positive values set the fractional mass loss
         else
            mrem = (1.d0 + rembar_massloss) * mt
         endif
      endif

      end


      SUBROUTINE assign_remnant_mandel_muller(mc, mc_tot, mt)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real ran3
      EXTERNAL ran3
      EXTERNAL RandomTruncatedNormal

      real*8 mc, mc_tot, mt

      real*8 mm_m1, mm_m2, mm_m3, mm_m4, min_ns_mass
      real*8 pBH, pCF
      real*8 u_pBH, u_pCF, u_rem
      real*8 ns_mu, ns_sigma

* Use the Mandel & Mueller 2020 prescription
*
      pBH = 0.0d0
      mm_m1 = 2.0d0
      mm_m2 = 3.0d0
      mm_m3 = 7.0d0
      mm_m4 = 8.0d0
      min_ns_mass = 1.13d0

* Determine probability of forming a BH based on core mass
      if(mc.lt.mm_m1)then
         pBH = 0.0d0
      elseif(mc.ge.mm_m1.and.mc.lt.mm_m3)then
         pBH = (mc - mm_m1)/(mm_m3 - mm_m1)
      else
         pBH = 1.0d0
      endif

* Draw random numbers for BH/NS decision, complete fallback decision
* and remnant mass assignment
      u_pBH = ran3(idum1)
      u_pCF = ran3(idum1)
      u_rem = ran3(idum1)

* Determine remnant type based on pBH
      if(u_pBH.le.pBH)then
* BH formed
         pCF = 0.0d0
         if (mc.ge.mm_m1.and.mc.lt.mm_m4) then
            pCF = (mc - mm_m1)/(mm_m4 - mm_m1)
         else
            pCF = 1.0d0
         endif

         if (u_pCF.le.pCF) then
* Complete fallback occurred, remnant mass equals pre-SN core mass
            mt = mc_tot
         else
* Partial fallback occurred, remnant mass drawn from Normal,
* but truncated to always be between max NS mass and CO core mass
            call RandomTruncatedNormal(0.8d0 * mc, 0.5d0 * 0.5d0,
     &                                 idum1, mxns, mc, mt)
         endif
      else
* NS formed, determine mu and sigma for random normal draw
* normal is truncated to be between min/max NS mass
         if (mc.lt.mm_m1) then
             ns_mu = 1.2d0
             ns_sigma = 0.02d0
         elseif (mc.ge.mm_m1.and.mc.lt.mm_m2) then
             ns_mu = 1.4d0 + 0.5d0 * (mc - mm_m1) / (mm_m2 - mm_m1)
             ns_sigma = 0.05d0
         else
             ns_mu = 1.4d0 + 0.4d0 * (mc - mm_m2) / (mm_m3 - mm_m2)
             ns_sigma = 0.05d0
         endif
         call RandomTruncatedNormal(ns_mu, ns_sigma, idum1, min_ns_mass,
     &                              mxns, mt)
      endif

      end


      SUBROUTINE assign_remnant_maltsev(mc, mc_tot, met, kidx, kw, mt)
* Assign remnant mass using the Maltsev et al. 2025 prescription,
* with additional details from Willcox et al. 2025
*
* Inputs:
*       mc         : CO core mass before SN
*       mc_tot     : Total core mass before SN (CO + He layers)
*       met        : Metallicity of the star
*       kidx       : Index of the star
*       kw         : Stellar (remnant) type
*
* Outputs:
*       mt         : Remnant mass after SN


      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      common /fall/fallback
      REAL*8 fallback

      real ran3
      EXTERNAL ran3

      real*8 mc, mc_tot, mt, met, u_NS
      real*8 log10Z_bounded
      integer mt_type, kidx, kw
      integer first_mt_type_as_donor
      EXTERNAL first_mt_type_as_donor

* Use the Maltsev+25 prescription with additional details from Willcox+25

      real*8 M1, M2, M3
      real*8 M1S, M2S, M3S, M1S_Z01, M2S_Z01, M3S_Z01
      real*8 M1A, M2A, M3A, M1A_Z01, M2A_Z01, M3A_Z01
      real*8 M1B, M2B, M3B, M1B_Z01, M2B_Z01, M3B_Z01
      real*8 M1C, M2C, M3C, M1C_Z01, M2C_Z01, M3C_Z01
      real*8 M_min, M_max, M_NS
      real*8 log10_1, log10_1_div_10, log10_1_div_50

* Taken from summary in Table A1 of Willcox+25
      PARAMETER(M1S=6.6d0, M2S=7.2d0, M3S=13.0d0,
     &          M1A=7.4d0, M2A=8.4d0, M3A=15.4d0,
     &          M1B=7.7d0, M2B=8.3d0, M3B=15.2d0,
     &          M1C=6.6d0, M2C=7.1d0, M3C=13.2d0,
     &          M1S_Z01=6.1d0, M2S_Z01=6.6d0, M3S_Z01=12.9d0,
     &          M1A_Z01=7.0d0, M2A_Z01=7.4d0, M3A_Z01=13.7d0,
     &          M1B_Z01=6.9d0, M2B_Z01=7.9d0, M3B_Z01=13.7d0,
     &          M1C_Z01=6.3d0, M2C_Z01=7.1d0, M3C_Z01=12.3d0,
     &          M_min=5.62d0, M_max=16.18d0,
     &          log10_1=0, log10_1_div_10=-1, log10_1_div_50=-1.69897d0,
     &          M_NS=1.4d0)

* Save some computation by immediately producing either a NS or a direct
collapse BH if the CO core mass is outside the Maltsev+25 range
      if(mc.lt.M_min)then
         fallback = 0.d0
         mt = M_NS
         return
      elseif(mc.gt.M_max)then
         fallback = 1.d0
         mt = mc_tot
         return
      endif

* Determine the mass transfer type of the donor star at first mass transfer
      mt_type = first_mt_type_as_donor(kidx)

* If star has not undergone mass transfer as donor, but has
* self-stripped (kw in [7,8,9]), assume case B mass transfer
      if(mt_type.eq.-1.and.(kw.eq.7.or.kw.eq.8.or.kw.eq.9))then
         mt_type = 1
      endif

* Normalize metallicity to solar
      met = met / zsun

* Bound log10Z based on maltsev_mode choice from user
      if (maltsev_mode.eq.0) then
        log10Z_bounded = log10(met)
      elseif (maltsev_mode.eq.1) then
        log10Z_bounded = min(max(log10(met), log10_1_div_50), log10_1)
      elseif (maltsev_mode.eq.2) then
        log10Z_bounded = min(max(log10(met), log10_1_div_10), log10_1)
      endif

* Determine the mass boundaries based on MT type and metallicity
      if(mt_type.eq.-1)then
         M1 = M1S + (M1S - M1S_Z01) * log10Z_bounded
         M2 = M2S + (M2S - M2S_Z01) * log10Z_bounded
         M3 = M3S + (M3S - M3S_Z01) * log10Z_bounded
      elseif(mt_type.eq.0)then
         M1 = M1A + (M1A - M1A_Z01) * log10Z_bounded
         M2 = M2A + (M2A - M2A_Z01) * log10Z_bounded
         M3 = M3A + (M3A - M3A_Z01) * log10Z_bounded
      elseif(mt_type.eq.1)then
         M1 = M1B + (M1B - M1B_Z01) * log10Z_bounded
         M2 = M2B + (M2B - M2B_Z01) * log10Z_bounded
         M3 = M3B + (M3B - M3B_Z01) * log10Z_bounded
      elseif(mt_type.eq.2)then
         M1 = M1C + (M1C - M1C_Z01) * log10Z_bounded
         M2 = M2C + (M2C - M2C_Z01) * log10Z_bounded
         M3 = M3C + (M3C - M3C_Z01) * log10Z_bounded
      endif

* Determine fallback and remnant mass based on CO core mass
*   Direct collapse:    if m1 <= mc <= m2 or mc >= m3
*   Partial fallback:   maltsev_pf_prob fraction of the time if m2 < mc < m3
*   NS:                 otherwise
      u_NS = ran3(idum1)
      if ((mc.ge.M1.and.mc.le.M2).or.(mc.ge.M3)) then
         fallback = 1.0d0
         mt = mc_tot
      elseif (mc.gt.M2.and.mc.lt.M3.and.u_NS.lt.maltsev_pf_prob) then
         fallback = maltsev_fallback
         mt = (mc_tot - M_NS) * fallback + M_NS
      else
         fallback = 0.0d0
         mt = M_NS
      endif

      return
      end


      SUBROUTINE assign_remnant_spin(mc, bhspin)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'

      real*8 ran3, mc, bhspin
      EXTERNAL ran3

* Set all BH spins equal to bhspinmag
      if(bhspinflag.eq.0)then
         bhspin = bhspinmag
* Randomly assign BH spins between 0 and bhspinmag
      elseif(bhspinflag.eq.1)then
         bhspin = ran3(idum1) * bhspinmag
* Assign BH spins based on Belczynski+17 prescription
      elseif(bhspinflag.eq.2)then
         if(mc.le.13.d0)then
            bhspin = 0.9d0
         elseif(mc.lt.27.d0)then
            bhspin = -0.064d0*mc + 1.736d0
         else
            bhspin = 0.0d0
         endif
      endif

      end

      integer function first_mt_type_as_donor(star)
* Find the first occurrence of mass transfer where a star is the donor
* and return the type
*
*     inputs:
*       star : integer (1 or 2) indicating which star to check
*
*     return value:
*       mt_type : integer indicating the stellar type at first mass
*                 transfer as donor for the specified star
*                 -1 if no mass transfer as donor found
*                 0 if case A mass transfer found
*                 1 if case B mass transfer found
*                 2 if case C mass transfer found
*
      INCLUDE 'const_bse.h'

      integer star
      integer i, col, kstar
      integer kstar1_col, kstar2_col, evol_type_col

* This part is tricky, because we allow users to change bpp_columns,
* nothing is in a fixed order. So we need to find them!
* Determine which column to find kstar1, kstar2, evolve_type based on
* col_inds_bpp. kstar1 col is the index of col_inds_bpp where value = 4,
* kstar2 same with 5, evol_type same with 11
      kstar1_col = -1
      kstar2_col = -1
      evol_type_col = -1
      do 5 i = 1, n_col_bpp
         if (col_inds_bpp(i).eq.4) then
            kstar1_col = i
         else if (col_inds_bpp(i) .eq. 5) then
            kstar2_col = i
         else if (col_inds_bpp(i) .eq. 11) then
            evol_type_col = i
         endif
   5  continue

      if (kstar1_col.eq.-1.or.kstar2_col.eq.-1..or.
     &    evol_type_col.eq.-1) then
         WRITE(*,*) 'Error in first_mt_type_as_donor: could not find '//
     &      'necessary columns in bpp (kstar1, kstar2, evol_type)'
         first_mt_type_as_donor = -1
         return
      endif

* If either star is a massless remnant (i.e. this is a merger product
* that's becoming a remnant), treat as no mass transfer and return -1
      if (int(bpp(bpp_ind,kstar1_col)).eq.15.or.
     &    int(bpp(bpp_ind,kstar2_col)).eq.15) then
            first_mt_type_as_donor = -1
            return
      endif

* work out which kstar column to use based on which star this is
      if (star.eq.1) then
         col = kstar1_col
      else if (star.eq.2) then
         col = kstar2_col
      endif

      kstar = -1
      do 10 i = 1, bpp_ind
         if (int(bpp(i,evol_type_col)).eq.3
     &   .or.int(bpp(i,evol_type_col)).eq.7) then
            kstar = int(bpp(i,col))
            goto 20
         endif
   10 continue

   20 if (kstar.eq.-1) then
         first_mt_type_as_donor = -1
      else if (kstar.ge.0.and.kstar.le.1.or.kstar.eq.7) then
         first_mt_type_as_donor = 0
      else if (kstar.eq.2.or.kstar.eq.8) then
         first_mt_type_as_donor = 1
      else if (kstar.ge.3.and.kstar.le.6.or.kstar.eq.9) then
         first_mt_type_as_donor = 2
      else
         first_mt_type_as_donor = -1
      endif

      return
      end
