import sys, os, time, random
import synergia
import synergia_workflow
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
from scipy import constants

import rssynergia
from rssynergia.standard import standard_beam
from rssynergia.base_diagnostics import read_bunch
from rssynergia.base_diagnostics import workflow
from rssynergia.base_diagnostics import latticework

def toyKVbeam6D(opts):
    '''Generate a costing KV beam with fixed Hamiltonian and returns the particle array.

    Coordinates are chosen with fixed random number generator seed, so that they should always
    produce the same initial distribution for a given emittance.

    Args:
        opts (object): Instance of the Synergia options class containing beam and macroparticle information

    Returns:
        bunch (ndarray): NumPy array with bunch with values {x, px, y, py, cdt, z, ID} for each particle
    '''

    # Beam parameters
    gamma0 = opts.gamma
    speciesmass = constants.m_p
    dgammaOgamma = 0
    #We want dpop no dE/E
    dpop = opts.dpop
    #Assume Gaussian longitudinal profile - put bunch length in m
    sigmaz = opts.stdz
    num_macro_particles = opts.macro_particles


    xOffset = 0. #m
    yOffset = 0. #m
    filename = opts.bunch_file

    # Do not modify below this line
    E0 = gamma0 * constants.m_p * constants.c**2 * 6.24150934e9 #GeV/J
    espread = E0 * dgammaOgamma

    #EArray = [0.]*numMacroParticles
    pArray = [0.]*num_macro_particles
    #tArray = [0.]*numMacroParticles
    cdtArray = [0.]*num_macro_particles

    #fix a random seed!
    random.seed(opts.seed)

    bunch = []
    #innerBunch = np.zeros(numMacroParticles)

    for index,emit in enumerate(opts.emits):

        innerbunch = np.zeros(num_macro_particles) #bunch at this emittance
        transverse_emittance = emit

        if opts.betae:
            mybunchgenerator = standard_beam.StandardBeam(opts.betae) #use betae
        else:
            mybunchgenerator = standard_beam.StandardBeam() #fixed beta=1, betaprime=0
        #coords is an array of 4-vectors containing coordinate space information
        coords = mybunchgenerator.generatefixedbunch(transverse_emittance, num_macro_particles, opts.seed)


        for idx in range(len(coords)):
            pArray[idx] = random.gauss(0, dpop)
            cdtArray[idx] = random.gauss(0, sigmaz)

            #assign unique index value to each particle
            ID = index*len(coords) + idx

            coords[idx] = np.append(coords[idx],[cdtArray[idx],pArray[idx],ID])

        bunch.append(coords)


    tofile = None
    if tofile:

        if os.path.isfile(filename):
            newfilename = filename+str(int(time.time()))
            print ' !Warning -- '
            print 'File '+filename+' already exists. Renaming the old file to '+newfilename
            os.rename('./'+filename, './'+newfilename)

        bunchfile = open(filename, 'w')
        for idx in range(0,num_macro_particles):
            ptclstring = str(bunch[idx][0])+' '+str(bunch[idx][1])+' '+str(bunch[idx][2])+' '+str(bunch[idx][3])+' '+str(cdtArray[idx])+' '+str(pArray[idx])+'\n'
            bunchfile.write(ptclstring)

        bunchfile.close()

    return np.asarray(bunch)
