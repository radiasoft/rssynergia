import sys, os, time, random
import synergia
import synergia_workflow
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
from scipy import constants

from semi_gaussian import SemiGaussianBeam
from base_diagnostics import read_bunch
from base_diagnostics import workflow
from base_diagnostics import latticework



def semigaussianbeam6D(opts):
    '''
    
    Generate a costing KV beam with fixed Hamiltonian and returns the particle array.
    
    Coordinates are chosen with fixed random number generator seed, so that they should always produce the same initial distribution
    for a given emittance.
    
    '''

    # Beam parameters
    gamma0 = opts.gamma
    speciesMass = constants.m_p
    dgammaOgamma = 0
    #We want dpop no dE/E
    dpop = opts.dpop
    #Assume Gaussian longitudinal profile - put bunch length in m
    sigmaz = opts.stdz
    numMacroParticles = opts.macro_particles
    

    xOffset = 0. #m
    yOffset = 0. #m
    fileName = opts.bunch_file

    # Do not modify below this line
    E0 = gamma0 * constants.m_p * constants.c**2 * 6.24150934e9 #GeV/J
    ESpread = E0 * dgammaOgamma

    #EArray = [0.]*numMacroParticles
    pArray = [0.]*numMacroParticles
    #tArray = [0.]*numMacroParticles
    cdtArray = [0.]*numMacroParticles

    #fix a random seed!
    random.seed(opts.seed)

    bunch = []
    #innerBunch = np.zeros(numMacroParticles)

    for index,emit in enumerate(opts.emits):
        
        innerBunch = np.zeros(numMacroParticles) #bunch at this emittance
        transverseEmittance = emit

        if opts.betae:
            myBunchGenerator = SemiGaussianBeam(opts.betae) #use betae
        else:
            myBunchGenerator = SemiGaussianBeam() #fixed beta=1, betaprime=0
        #coords is an array of 4-vectors containing coordinate space information
        coords = myBunchGenerator.generatefixedbunch(transverseEmittance, numMacroParticles, opts.seed)
        
        
        lc = coords.shape[0]
        print lc
        if dpop == 0:
            pArray = np.zeros(lc)
        else:
            pArray = np.random.standard_normal(lc)*dpop
            
        cdtArray = np.random.standard_normal(lc)*sigmaz
        ID = index*lc + np.arange(lc)
        
        #print pArray.shape
        #print cdtArray.shape
        #print ID.shape
        coords = np.vstack((coords[:,0],coords[:,1],coords[:,2],coords[:,3],cdtArray,pArray,ID))
        
        #coords[idx] = np.append(coords[idx],[cdtArray[idx],pArray[idx],ID])
        
        bunch.append(coords.T)
        

    toFile = None
    if toFile:

        if os.path.isfile(fileName):
            newFileName = fileName+str(int(time.time()))
            print ' !Warning -- '
            print 'File '+fileName+' already exists. Renaming the old file to '+newFileName
            os.rename('./'+fileName, './'+newFileName)
    
        bunchFile = open(fileName, 'w')
        for idx in range(0,numMacroParticles):
            ptclString = str(bunch[idx][0])+' '+str(bunch[idx][1])+' '+str(bunch[idx][2])+' '+str(bunch[idx][3])+' '+str(cdtArray[idx])+' '+str(pArray[idx])+'\n'
            bunchFile.write(ptclString)

        bunchFile.close()
        
    return np.asarray(bunch)