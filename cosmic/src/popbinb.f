
      SUBROUTINE runPopbin(binaryFile, binID, binTYPE, nCore, 
     &                      nCoreWrite)
      
      INCLUDE 'const_bse.h'
      
      integer i,j,k,jj,kk,ll,mm, timestep,timestepCalc
      integer kw,kw2,kwx,kwx2,kstar(2)
      integer i1,i2,kdum,nSteps,fileNum
*
      real*8 m1,m2,tmax,timeEvolved,maxEpochTime
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 sep0,tb0,tb,ecc0,ecc,aursun,yeardy,yearsc,tol
      PARAMETER(aursun=214.95d0,yeardy=365.25d0,yearsc=3.1557d+07)
      PARAMETER(tol=1.d-07)
      real*8 mx,mx2,tbx,eccx
      CHARACTER*8 label(14)
      CHARACTER str1,str2,str3,str4
      REAL*8 B_0(2),bacc(2),tacc(2),bkick(12)

      
* DEFINE PARAMETERS FOR POPULATION SIM:
      integer binaryCountMax,time,key,Ncores,NbinEpoch,Nepochs,GxFlag
      integer HeHeFlag,HeCOFlag,HeONeFlag,COCOFlag,saveFileFlag
      integer COONeflag,ONeONeflag,NSHeflag,NSCOflag,NSONflag,LISAflag  
      integer NSNSflag,NSBHflag,BHBHflag,BHHeflag, BHCOflag, BHONeflag
      parameter(binaryCountMax=1000000000)
      real*8 EpochLength,SFR,startEpoch,met,binwidth,birthTime
      real*8 localTime(80),tBorn1,tBorn2
      character*46 binaryFile,SFH
      integer nCore,nBinEvolve,binID,bID,mID,aID,diID
      character*4 binTYPE
      character*2 nCoreWrite 
      character*46 binBorn, binMerged, binAlive, mtOnset
      real*8 porb,sep,m1i,m2i,lum1,lum2,rad1,rad2,Temp1,Temp2
      integer popCount
      character*6 popCountWrite


          
*
************************************************************************
* BSE parameters:
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (1.0 normally).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.5).
*
* ceflag = 1 activates spin-energy correction in common-envelope (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (0). 
* bhflag > 0 allows velocity kick at BH formation (0). 
* nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
* mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*

*
* Open the input file - list of binary initial parameters. 
*
      OPEN(10,file='BSEparams.in',status='old')
*      
* Read out the header lines
*      
      DO i=1,47
         READ(10,*)
      ENDDO
*      
* Read in the BSE parameters      
*
      READ(10,*)neta,bwind,hewind,alpha1,lambda  
      READ(10,*)ceflag,tflag,ifflag,wdflag
      READ(10,*)bhflag, nsflag, mxns, pts1, pts2, pts3
      READ(10,*)sigma, beta, xi, acc2, epsnov, eddfac, gamma
      READ(10,*)bconst, CK, merger, windflag, fbkickswitch
      

*
* Set the seed for the random number generator. 
*
      idum = 3234
      if(idum.gt.0) idum = -idum
      
*
* Set the collision matrix.
*
      CALL instar

*
* Read the number of binaries to evolve from the Gx Params file.
*
      OPEN(11,file='Params.in',status='old')
*      
* Read out the header lines
*      
      DO i=1,33
         READ(11,*)
      ENDDO    
*      
* Read in the simulation parameters      
*       
      READ(11,*)NbinEpoch,Nepochs,EpochLength,Ncores,Ngalaxies,GxFlag
      READ(11,*)met,SFH,SFR,binwidth,LISAflag,saveFileFlag
      READ(11,*)HeHeFlag,HeCOFlag,HeONeFlag,COCOFlag,NSHeflag,
     &          NSCOflag,NSONeflag
      READ(11,*)NSNSflag,NSBHflag,BHBHflag,BHHeflag,
     &          BHCOflag,BHONeflag
      
*
* Compute the number of binaries per binaryFile
*
      nBinEvolve = INT(NbinEpoch*(Nepochs))
      
*
* Open the initial binary parameter file
* 
      OPEN(12,file=binaryFile,status='old')
 
      
*     
* create the files to output the binaries that are born, merged, and alive at epoch end
*      
      bID = binID+60
      binBorn='bulgeTmp/lambda10/born'//binTYPE//
     &         '_'//nCoreWrite//'.dat'
      OPEN(bID,file=binBorn, status='unknown')
      WRITE(bID,*)'#TIME    K1   K2      '
      
      mID = binID+70
      binMerged='bulgeTmp/lambda10/merged'//binTYPE//
     &           '_'//nCoreWrite//'.dat'
      OPEN(mID,file=binMerged, status='unknown')
      WRITE(mID,*)'#TIME    K1   K2      '
      

      diID = binID+90
      mtOnset = 'bulgeTmp/lambda10/MTonset'//
     &           binTYPE//'_'//nCoreWrite//'.dat'
      OPEN(diID, file=mtOnset, status='unknown')
      WRITE(*,*)mtOnset
      WRITE(diID, *)'#BIN#,TIME,M1,M2,PERIOD,SEP,ECC'

*
* Initialize the binary count
*
      popCount = 0

      WRITE(*,*)'The number of binaries is:',nBinEvolve
*
* Loop over all binaries in binaries.in file and evolve  
*           
      DO i = 1,(nBinEvolve)
*
* Read in parameters and set coefficients which depend on metallicity. 
*
         READ(12,*)startEpoch,m1,m2,tb,ecc,z,tmax,maxEpochTime
         CALL zcnsts(z,zpars)
           
         label(1) = 'INITIAL '
         label(2) = 'KW CHNGE'
         label(3) = 'BEG RCHE'
         label(4) = 'END RCHE'
         label(5) = 'CONTACT '
         label(6) = 'COELESCE'
         label(7) = 'COMENV  '
         label(8) = 'GNTAGE  '
         label(9) = 'NO REMNT'
         label(10) = 'MAX TIME'
         label(11) = 'DISRUPT '
         label(12) = 'BEG SYMB'
         label(13) = 'END SYMB'
         label(14) = 'BEG BSS'
*
         ecc0 = ecc
         tb0 = tb/yeardy
         sep0 = aursun*(tb0*tb0*(mass(1) + mass(2)))**(1.d0/3.d0)
         tb0 = tb
*
* Initialize the binary. 
*
         kstar(1) = 1
         mass0(1) = m1
         mass(1) = m1
         massc(1) = 0.0
         ospin(1) = 0.0
         epoch(1) = 0.0
*
         kstar(2) = 1
         mass0(2) = m2
         mass(2) = m2
         massc(2) = 0.0
         ospin(2) = 0.0
         epoch(2) = 0.0
*
         tphys = 0.0
         tphysf = tmax
         dtp = tphysf
         
*
* Initialize the kick matrix
*

          DO jj = 1,12
             bkick(jj) = 0
          ENDDO
             
*
* Initialize the magnetic braking parameters. These are always 0 initially
*
          
          B_0(1) = 0.0
          B_0(2) = 0.0
          bacc(1) = 0.0
          bacc(2) = 0.0
          tacc(1) = 0.0
          tacc(2) = 0.0
          
*
* Evolve the binary. 
*
         CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &               menv,renv,ospin,epoch,tms,B_0,bacc,tacc,
     &               tphys,tphysf,dtp,z,zpars,tb,ecc,bkick)

*


************************************************************************
* Output:
* First check that bcm is not empty.
*
         if(bcm(1,1).lt.0.0) goto 50
*
* The bcm array stores the stellar and orbital parameters at the 
* specified output times. The parameters are (in order of storage):
*
*    Time, 
*    [stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(Teff), core mass, core radius, mass of any convective 
*    envelope, radius of the envelope, epoch, spin, mass loss rate and 
*    ratio of radius to roche lobe radius (repeated for secondary)],
*    period, separation, eccentricity.
*
         mm=0
         j = 0
 30      j = j + 1
         if(bcm(j,1).lt.0.0)then
            bcm(j-1,1) = bcm(j,1)
            j = j - 1
         endif
         kw = INT(bcm(j,2))
         kw2 = INT(bcm(j,16))
         time = bcm(j,1)
         key = INT(bpp(j,10))
         m1i = bcm(j,3)
         m2i = bcm(j,17)
         m1 = bcm(j,4)
         m2 = bcm(j,18)
         porb = bcm(j,30)
         sep = bcm(j,31)
         ecc = bcm(j,32)
         lum1 = bcm(j,5)
         rad1 = bcm(j,6)
         Temp1 = bcm(j,7)
         lum2 = bcm(j,19)
         rad2 = bcm(j,20)
         Temp2 = bcm(j,21)
         
          
*
* The bpp array acts as a log, storing parameters at each change
* of evolution stage.
*
 50      j = 0
         ll=0
         kk=0
         mm=0
     
 52      j = j + 1
         
         if(bpp(j,1).lt.0.0) goto 60
         kstar(1) = INT(bpp(j,4))
         kstar(2) = INT(bpp(j,5))
         kw = INT(bpp(j,10))
         localTime(j) = bpp(j,1)+startEpoch
         tBorn1 = bpp(j,11)+startEpoch
         tBorn2 = bpp(j,12)+startEpoch
         
         if(binID.eq.0)then
            if(kstar(1).eq.10.and.kstar(2).eq.10)then
               if(bpp(j,7).ne.-1.00.and.bpp(j,7).lt.1.00)then
                  if(kk.eq.0)then
                     WRITE(bID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(bID)
                     birthTime = localTime(j)
                     kk = kk+1
                  end if
               end if
               if(INT(bpp(j+1,10)).eq.6.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(mID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(mID)
                  end if
               end if
               if(INT(bpp(j,10)).eq.3.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(diID,129)i,localTime(j),bpp(j,2),bpp(j,3),
     &                             bpp(j,6),bpp(j,7),bpp(j,8),
     &                             bpp(j,9),tBorn1,tBorn2,birthTime
                     FLUSH(diID)
                  end if
               end if
            end if
         else if(binID.eq.1)then
            if((kstar(1).eq.10.and.kstar(2).eq.11)
     &         .or.(kstar(2).eq.10.and.kstar(1).eq.11))then
               if(bpp(j,7).ne.-1.00.and.bpp(j,7).lt.1.00)then
                  if(ll.eq.0)then
                     WRITE(bID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(bID)
                     birthTime = localTime(j)
                     ll = ll+1
                  end if
               end if
               if(INT(bpp(j+1,10)).eq.6.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(mID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(mID)
                  end if 
               end if
               if(INT(bpp(j,10)).eq.3.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(diID,129)i,localTime(j),bpp(j,2),bpp(j,3),
     &                          bpp(j,6),bpp(j,7),bpp(j,8),bpp(j,9),
     &                          tBorn1,tBorn2,birthTime
                     FLUSH(diID)
                  end if
               end if
            end if
         else if(binID.eq.2)then
            if((kstar(1).eq.10.and.kstar(2).eq.12)
     &         .or.(kstar(2).eq.10.and.kstar(1).eq.12))then
               if(bpp(j,7).ne.-1.00.and.bpp(j,7).lt.1.00)then
                  if(ll.eq.0)then
                     WRITE(bID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(bID)
                     birthTime = localTime(j)
                     ll = ll+1
                  end if
               end if
               if(INT(bpp(j+1,10)).eq.6.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(mID,128)i,localTime(j),kstar(1),kstar(2)
                     FLUSH(mID)
                  end if 
               end if
               if(INT(bpp(j,10)).eq.3.and.bpp(j,7).ne.-1.00)then
                  if(localTime(j).lt.maxEpochTime)then
                     WRITE(diID,129)i,localTime(j),bpp(j,2),bpp(j,3),
     &                          bpp(j,6),bpp(j,7),bpp(j,8),bpp(j,9),
     &                          tBorn1,tBorn2,birthTime
                     FLUSH(diID)
                  end if
               end if
            end if
         end if

        
 128     FORMAT(i9,f16.4,2i5)
 129     FORMAT(i9,10f16.4)
         goto 52
 60      continue
*100     FORMAT(i7,f11.4,2f9.3,2i3,f13.3,f6.2,2f8.3,2x,a8)

      
************************************************************************

      ENDDO
*
      
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(bID)
      CLOSE(mID)
      CLOSE(diID) 
      RETURN      
      
      END SUBROUTINE runPopbin
      
