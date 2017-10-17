#! /usr/bin/env python

# Code: sampleBinaries.py
# Version: 1
# Version changes: SAMPLE INITIAL BINARIES AND EVOLVE USING BSE
#
# Input: coreNumber, nBin to draw per epoch, range of epochs, epochLength,
# Output: No output; writes to file
#
# Edited on: 9 Sep 2015


import numpy as np
import math
import random
import time
from aCOSMIC import _popbinb as _popbin


##################################################################################
# FUNCTION TO SAMPLE THE INITIAL BINARIES
##################################################################################
def sample(nCore, binID, binTYPE, nBin, rangeEpoch, epochLength, SFH, maxTime, startTime, met):

    # SET THE RANDOM SEED
    ##############################################################################    
    np.random.seed()

    # SET TIME TO TRACK COMPUTATION TIME
    ##############################################################################
    start_time = time.time()
    
    # CONSTANTS
    ##############################################################################
    G = 6.67384*math.pow(10, -11.0)
    c = 2.99792458*math.pow(10, 8.0)
    Rsun = 6.955*math.pow(10, 8)
    Msun = 1.9891*math.pow(10,30)
    sec_in_day = 86400.0
   
    # DECLARE NUMPY ARRAYS TO HOLD THE BINARY PARAMETERS
    ##############################################################################
    mass1_pass = []
    mass2_pass = []
    tb_pass = []
    ecc_pass = []
    met_pass = []
    tphysf_pass = [] 
    fixedPopLog = []
    initialParams = []
    epochStart_pass = []
    epochMax_pass = []

    # ASSIGN MINIMUM TOTAL BINARY MASS ACCORDINGLY
    ##############################################################################
    if binID == 0:
        minMass = 0.0
        maxMass = 4.0
    elif binID == 1:
        minMass = 0.0
        maxMass = 8.0
    elif binID == 2:
        minMass = 1.0
        maxMass = 10.0
    elif binID == 3:
        minMass = 0.0
        maxMass = 8.0
    elif binID == 4:
        minMass = 4.0
        maxMass = 20.0
    elif binID == 5:
        minMass = 0.0
        maxMass = 20.0
    elif binID == 6:
        minMass = 10.0
        maxMass = 50.0
    elif binID == 7:
        minMass = 10.0
        maxMass = 50.0
    elif binID == 8:
        minMass = 10.0
        maxMass = 50.0
    elif binID == 9:
        minMass = 10.0
        maxMass = 50.0
    elif binID == 10:
        minMass = 10.0
        maxMass = 150.0
    elif binID == 11:
        minMass = 10.0
        maxMass = 150.0 
    elif binID == 12:
        minMass = 10.0
        maxMass = 150.0
    elif binID == 13:
        minMass = 10.0
        maxMass = 150.0
    elif binID == 14:
        minMass = 6.0
        maxMass = 150.0
        
    kk = 0
    ll = 0
    for i in rangeEpoch:
        
        Nbinsave = 0
        nBinTot = 0
        Nstars = 0
        mass_bin_fixed = 0
        mass_total_fixed = 0

        # SET EVOLUTION TIME OF FIXED POPUPLATION
        ##############################################################################      
        epochNumber = i
        
        # IF rangeEpoch[1] == 1, THIS REPRESENTS THE DELTA FUNCTION BURST
        ##############################################################################      
        
        if SFH == " burst":
            EpochTime = epochLength
            EpochStart = startTime
        # IF NOT, WE USE THE USER SPECIFIED SFR
        ############################################################################## 
        else:    
            EpochTime = epochNumber*epochLength
            EpochStart = maxTime - EpochTime
        while Nbinsave < nBin:
            
            # THE EVOLUTION TIME IS CONSTANT FOR EACH EPOCH OF STAR FORMATION: 
            # 'COEVAL' BINARIES
            #########################################################################
            tphysf = maxTime 
            
            # ALL PDFs TAKEN FROM Yu&Jeffery (2015) EXCEPT THE IMF TAKEN 
            # FROM Robin et al (2003)
            # IMF ~ m^(-2.35) m>0.7
            #########################################################################
            a_0 = np.random.uniform(0.0, 1.0)
            normFac = 1.19742                 
            mass_draw = ((0.7)**(-27.0/20.0)-(a_0*normFac/0.74071))**(-20.0/27.0)
            
            # COMPUTE BINARY FRACTION ACCORDING TO van Haaften et al.(2009) in appdx 
            #########################################################################
            binary_fraction = 1/2.0 + 1/4.0*math.log(mass_draw, 10)         
         
            binary_choose =  np.random.uniform(0, 1.0)
            if binary_choose < binary_fraction:
                nBinTot = nBinTot+1
                mass1 = mass_draw
            
                # SAMPLE MASS RATIO: uniform from 0.001 to 1
                ######################################################################
                # Generate a random number drawn uniformly from 0.001 to 1 
                mass_ratio = np.random.uniform(0.001, 1)
    
                # Compute the secondary from the mass ratio: ratio = M2/M1
                mass2 = mass_ratio*mass1
                 
                    
                # SAMPLE SEPARATION:
                ######################################################################
                a_0_sep = np.random.uniform(0, 1)
                if a_0_sep < 0.0583333:
                    separation = math.pow( (a_0_sep/0.00368058), 5/6.0)
                else:
                    separation = math.exp(a_0_sep/0.07+math.log(10.0))
        
                sep_m = separation*Rsun
        
                # Compute the orbital period in seconds using NEWTON III:
                porb_sec = pow( 4*pow(math.pi, 2)/(G*(mass1+mass2)*Msun)*\
                                 pow(sep_m, 3.0) , 0.5)
                # Convert the orbital period from seconds to days for the BSE input:
                tb = porb_sec/sec_in_day
        
        
                # SAMPLE ECCENTRICITY:
                ######################################################################
                a_0_ecc = np.random.uniform(0.0, 1.0)
                ecc = a_0_ecc**0.5
                
                # APPEND THE SAMPLED VALUES ACCORDING TO POPULATION FLAGS
                ######################################################################                   
            	mass_drawn = mass1 + mass2 
            	mass_bin_fixed = mass_bin_fixed + mass_drawn
                 
                if mass_drawn > minMass:
                    Nbinsave = Nbinsave+1
                    epochStart_pass.append(EpochStart)
                    epochMax_pass.append(maxTime)
                    mass1_pass.append(mass1)
                    mass2_pass.append(mass2)
                    tb_pass.append(tb)
                    ecc_pass.append(ecc)
                    met_pass.append(met)
                    tphysf_pass.append(tphysf)
                         
            # next iterations!
            else:
                mass_drawn = mass_draw 
                Nstars = Nstars + 1    
    
            mass_total_fixed = mass_total_fixed + mass_drawn
            
        fixedPopLog.append((epochNumber, nBin, mass_total_fixed, mass_bin_fixed))
            
    
    # GATHER ALL INITIAL PARAMETERS
    ##############################################################################
    initialParams = np.vstack((epochStart_pass, mass1_pass, mass2_pass, tb_pass, ecc_pass, met_pass, tphysf_pass, epochMax_pass)).T
                               
    rowSize = len(initialParams)
                             
    # SAVE THE PARAMETERS TO EVOLVE
    ##################################################################################
    binariesInFile = "bulgeTmp/lambda10/binaries_"+binTYPE+"_"+str(nCore)+".in"
    outputfile = open(binariesInFile,"w",0)
    np.savetxt(outputfile, initialParams, delimiter = ",")
    outputfile.flush()
    outputfile.close()

    # SAVE THE FIXED POP LOG
    ##################################################################################
    fixedPopLogSave = np.vstack(fixedPopLog)
    np.savetxt("bulgeTmp/lambda10/fixedPopLog_"+binTYPE+"_"+str(nCore)+".dat", fixedPopLogSave, delimiter = ',')
    
    # EVOLVE THE BINARIES WITH POPBIN
    ##################################################################################
    nCoreWrite = str(nCore).zfill(2)
    _popbin.runpopbin(binariesInFile, binID, binTYPE, nCore, nCoreWrite)
    
    
    
    #print "The initial binary sample time from core "+str(nCore)+ \
          #" and type "+binTYPE+" is: ", time.time()-start_time, "seconds"
    
######################################################################################
######################################################################################


