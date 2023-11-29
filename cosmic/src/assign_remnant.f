***
      SUBROUTINE assign_remnant(zpars,mc,mcbagb,mass,mt,kw,bhspin,kidx)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      common /fall/fallback
      REAL*8 fallback
      REAL ran3
      EXTERNAL ran3
      real*8 zpars(20)

      real*8 avar,bvar
      real*8 mc,mcbagb,mass,mt
      real*8 frac,kappa,sappa,alphap,polyfit
      real*8 mcx, bhspin,mrem,mch
      integer kw,kidx

      
* input mc(or mcmax),mass, mcbagb
* output mt, kw, mc, mass,bhspin
* common: mch, mxns,
* ifflag,ecsn (or mc1), ecsn_low(mc2), bhspinflag,
*  bhspinmag,  Mbh_initial

            mch = 1.44d0 !set here owing to AIC ECSN model.
*
* Mc,x represents the C core mass and we now test whether it
* exceeds either the total mass or the maximum allowed core mass.
*
*         if(mcmax-mcx.lt.tiny)then
*            aj = 0.d0
*            mc = mcmax
            if(mc.lt.mch)then
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
     &                mcbagb.le.ecsn.and.mc.lt.1.08d0)then
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
                  if(remnantflag.eq.0)then
                     mt = 1.17d0 + 0.09d0*mc
                  elseif(remnantflag.eq.1)then
*
* Use NS/BH mass given by Belczynski et al. 2002, ApJ, 572, 407.
*
                     if(mc.lt.2.5d0)then
                        mcx = 0.161767d0*mc + 1.067055d0
                     else
                        mcx = 0.314154d0*mc + 0.686088d0
                     endif
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (mc - 5.d0)*(mt - mcx)/2.6d0
                        fallback = (mc - 5.d0)/2.6d0
                     elseif(mc.gt.7.60)then
                        fallback = 1.d0
                     endif
                  elseif(remnantflag.eq.2)then
*
* Use NS/BH masses given by Belczynski+08. PK.
*
                     !First calculate the proto-core mass
                     if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                    mcbagb.ge.ecsn_mlow)then
                        mcx = 1.38d0
*                     elseif(mc.lt.4.29d0)then
                     elseif(mc.lt.4.82d0)then
                        mcx = 1.5d0
*                     elseif(mc.ge.4.29d0.and.mc.lt.6.31d0)then
                     elseif(mc.ge.4.82d0.and.mc.lt.6.31d0)then
                        mcx = 2.11d0
                     elseif(mc.ge.6.31d0.and.mc.lt.6.75d0)then
                        mcx = 0.69*mc - 2.26d0
                     elseif(mc.ge.6.75d0)then
                        mcx = 0.37*mc - 0.07d0
                     endif
                     !now calculate the remnant mass after fallback
                     if(mc.le.5.d0)then
                        mt = mcx
                        fallback = 0.d0
                     elseif(mc.lt.7.6d0)then
                        mt = mcx + (0.378d0*mc - 1.889d0)*(mt - mcx)
                        fallback = (0.378d0*mc - 1.889d0)
                     elseif(mc.gt.7.60)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(remnantflag.eq.3)then
*
* Use the "Rapid" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    We use the updated proto-core mass from Giacobbo & Mapelli 2020
                     mcx = 1.1d0
                     if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                    mcbagb.ge.ecsn_mlow)then
                        mt = 1.38d0   ! ECSN fixed mass, no fallback
                     elseif(mc.le.2.5d0)then
                        fallback = 0.2d0 / (mt - mcx)
                        mt = mcx + 0.2d0
                     elseif(mc.le.6.d0)then
                        fallback = (0.286d0*mc - 0.514d0) / (mt - mcx)
                        mt = mcx + 0.286d0*mc - 0.514d0
                     elseif(mc.le.7.d0)then
                        fallback = 1.d0
                     elseif(mc.le.11.d0)then
                        avar = 0.25d0 - (1.275 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.gt.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  elseif(remnantflag.eq.4)then
*
* Use the "Delayed" SN Prescription (Fryer et al. 2012, APJ, 749,91)
*
*                    For this, we just set the proto-core mass to one
                     if(mc.le.3.5d0)then
                        mcx = 1.2d0
                     elseif(mc.le.6.d0)then
                        mcx = 1.3d0
                     elseif(mc.le.11.d0)then
                        mcx = 1.4d0
                     elseif(mc.gt.11.d0)then
                        mcx = 1.6d0
                     endif

                     if(ecsn.gt.0.d0.and.mcbagb.le.ecsn.and.
     &                    mcbagb.ge.ecsn_mlow)then
                        mt = 1.38d0   ! ECSN fixed mass, no fallback
                     elseif(mc.lt.2.5d0)then
                        fallback = 0.2d0 / (mt - mcx)
                        mt = mcx + 0.2
                     elseif(mc.lt.3.5d0)then
                        fallback = (0.5d0 * mc - 1.05d0) / (mt - mcx)
                        mt = mcx + 0.5d0 * mc - 1.05d0
                     elseif(mc.lt.11.d0)then
                        avar = 0.133d0 - (0.093d0 / (mt - mcx))
                        bvar = 1.d0 - 11.d0*avar
                        fallback = avar*mc + bvar
                        mt = mcx + fallback*(mt - mcx)
                     elseif(mc.ge.11.d0)then
                        fallback = 1.d0
                     endif
                     if(bhspinflag.eq.0)then
                            bhspin = bhspinmag
                     elseif(bhspinflag.eq.1)then
                            bhspin = ran3(idum1) * bhspinmag
                     elseif(bhspinflag.eq.2)then
                         if(mc.le.13.d0)then
                             bhspin = 0.9d0
                         elseif(mc.lt.27.d0)then
                             bhspin = -0.064d0*mc + 1.736d0
                         else
                             bhspin = 0.0d0
                         endif
                     endif
                     mc = mt
                  endif

* Specify the baryonic to gravitational remnant mass prescription
* MJZ 04/2020

* Determine gravitational mass using Lattimer & Yahil 1989 for remnantflag>1
                  if(remnantflag.le.1)then
                     mrem = mt
                  else
                     mrem = 6.6666667d0*(SQRT(1.d0+0.3d0*mt)-1.d0)
* If rembar_massloss >= 0, limit the massloss by rembar_massloss
                     if(rembar_massloss.ge.0d0)then
                        if((mt-mrem).ge.rembar_massloss)
     &                               mrem = mt-rembar_massloss
                     endif
                  endif

* Determine whether a zero-age NS or BH is formed
                  if(mrem.le.mxns)then
                     mt = mrem
                     kw = 13
                  else
                     kw = 14

* CLR - (Pulsational) Pair-Instability Supernova

* Belczynski+2016 prescription: just shrink any BH with a He core mass
* between 45 and 65 solar masses (provided the pisn flag is set at 45),
* and blow up anything between 65 and 135 solar masses.
* Cheap, but effective
                     if(pisn.gt.0)then
                        if(mcbagb.ge.pisn.and.mcbagb.lt.65.d0)then
                           mt = pisn
                           mc = pisn
                           pisn_track(kidx)=6
                        elseif(mcbagb.ge.65.d0.and.mcbagb.lt.135.d0)then
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
                        frac = mcbagb/mt
                        kappa = 0.67d0*frac + 0.1d0
                        sappa = 0.5226d0*frac - 0.52974d0
                        if(mcbagb.le.32.d0)then
                           alphap = 1.0d0
                        elseif(frac.lt.0.9d0.and.mcbagb.le.37.d0)then
                           alphap = 0.2d0*(kappa-1.d0)*mcbagb +
     &                              0.2d0*(37.d0 - 32.d0*kappa)
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mcbagb.le.60.d0)then
                           alphap = kappa
                           pisn_track(kidx)=6
                        elseif(frac.lt.0.9d0.and.mcbagb.lt.64.d0)then
                           alphap = kappa*(-0.25d0)*mcbagb + kappa*16.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.37.d0)then
                           alphap = sappa*(mcbagb - 32.d0) + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.56.d0.and.
     &                         sappa.lt.-0.034168d0)then
                           alphap = 5.d0*sappa + 1.d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.le.56.d0.and.
     &                         sappa.ge.-0.034168d0)then
                           alphap = (-0.1381d0*frac + 0.1309d0)*
     &                              (mcbagb - 56.d0) + 0.82916d0
                           pisn_track(kidx)=6
                        elseif(frac.ge.0.9d0.and.mcbagb.lt.64.d0)then
                           alphap = -0.103645d0*mcbagb + 6.63328d0
                           pisn_track(kidx)=6
                        elseif(mcbagb.ge.64.d0.and.mcbagb.lt.135.d0)then
                           alphap = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        elseif(mcbagb.ge.135.d0)then
                           alphap = 1.0d0
                        endif
                        mt = alphap*mt

* Fit (8th order polynomial) to Table 1 in Marchant+2018.
                     elseif(pisn.eq.-2)then
                        if(mcbagb.ge.31.99d0.and.mcbagb.le.61.10d0)then
                           polyfit = -6.29429263d5
     &                            + 1.15957797d5*mcbagb
     &                            - 9.28332577d3*mcbagb**2d0
     &                            + 4.21856189d2*mcbagb**3d0
     &                            - 1.19019565d1*mcbagb**4d0
     &                            + 2.13499267d-1*mcbagb**5d0
     &                            - 2.37814255d-3*mcbagb**6d0
     &                            + 1.50408118d-5*mcbagb**7d0
     &                            - 4.13587235d-8*mcbagb**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mcbagb.gt.61.10d0.and.
     &                         mcbagb.lt.124.12d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif

* Fit (8th order polynomial) to Table 5 in Woosley2019.
                     elseif(pisn.eq.-3)then
                        if(mcbagb.ge.29.53d0.and.mcbagb.le.60.12d0)then
                           polyfit = -3.14610870d5
     &                            + 6.13699616d4*mcbagb
     &                            - 5.19249710d3*mcbagb**2d0
     &                            + 2.48914888d2*mcbagb**3d0
     &                            - 7.39487537d0*mcbagb**4d0
     &                            + 1.39439936d-1*mcbagb**5d0
     &                            - 1.63012111d-3*mcbagb**6d0
     &                            + 1.08052344d-5*mcbagb**7d0
     &                            - 3.11019088d-8*mcbagb**8d0
                           mt = polyfit
                           pisn_track(kidx)=6
                        elseif(mcbagb.gt.60.12d0.and.
     &                         mcbagb.lt.135.d0)then
                           mt = 0.d0
                           kw = 15
                           pisn_track(kidx)=7
                        endif
                     endif

* Convert baryonic mass to gravitational mass
* MJZ 04/2020
                     if(remnantflag.le.1)then
                        mrem = mt
                     else
* If rembar_massloss >= 0, limit the massloss by rembar_massloss
                        if(rembar_massloss.ge.0d0)then
                           mrem = mt-rembar_massloss
* If -1 < rembar_massloss < 0, assume this fractional mass loss
                        else
                           mrem = (1.d0+rembar_massloss)*mt
                        endif
                     endif
                     mt = mrem
* Store the initial BH mass for calculating the ISCO later
                     if(Mbh_initial.eq.0)then
                        Mbh_initial = mt
                     endif
                  endif
               endif
            endif
*
      end
