#!/usr/bin/python

from elliptic_beam import EllipticBeam
from scipy import constants as consts
import random, time
import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow

import numpy as np
import matplotlib.pyplot as plt

import rssynergia
from rssynergia.base_diagnostics import read_bunch
from rssynergia.base_diagnostics import workflow
from rssynergia.base_diagnostics import latticework


def toyellipticalbeam6D(opts):
    '''

    Generate a toy costing beam with fixed elliptical Hamiltonian and returns the particle array, while also saving the array to
    a text file.

    Coordinates are chosen with fixed random number generator seed, so that they should always produce the same initial distribution
    for a given emittance.

    Args:
        opts (object): Instance of the Synergia options class containing beam and macroparticle information

    Returns:
        bunch (ndarray): NumPy array with bunch with values {x, px, y, py, cdt, z, ID} for each particle

    '''

    # Beam parameters
    gamma0 = 2.
    speciesMass = consts.m_p
    dgammaOgamma = 1.e-3
    #We want dpop no dE/E
    dpop = opts.dpop
    #Assume Gaussian longitudinal profile - put bunch length in m
    sigmaz = opts.stdz
    numMacroParticles = opts.macro_particles

    t = opts.t
    c = opts.c

    #Calculate bunch coordinates at entrance to NL magnet section
    beta = opts.betae
    alpha = -1*opts.alphae
    betaPrime = 2 * alpha

    xOffset = 0. #m
    yOffset = 0. #m
    fileName = opts.bunch_file

    # Do not modify below this line
    E0 = gamma0 * consts.m_p * consts.c**2 * 6.24150934e9 #GeV/J
    ESpread = E0 * dgammaOgamma

    pArray = [0.]*numMacroParticles
    cdtArray = [0.]*numMacroParticles

    #fix a random seed!
    random.seed(opts.seed)

    bunch = []
    #innerBunch = np.zeros(numMacroParticles)

    for index,emit in enumerate(opts.emits):

        innerBunch = np.zeros(numMacroParticles) #bunch at this emittance
        transverseEmittance = emit

        myBunchGenerator = EllipticBeam(t, c, beta, betaPrime)
        #coords is an array of 4-vectors containing coordinate space information
        coords = myBunchGenerator.generatefixedbunch(transverseEmittance, numMacroParticles, opts.seed)


        for idx in range(len(coords)):
            pArray[idx] = random.gauss(0, dpop)
            cdtArray[idx] = random.gauss(0, sigmaz)

            #assign unique index value to each particle
            ID = index*len(coords) + idx

            coords[idx] = np.append(coords[idx],[cdtArray[idx],pArray[idx],ID])

        bunch.append(coords)


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


def fixedellipticalbeam6D(opts):
    """
    Generate a toy costing beam with fixed elliptical Hamiltonian for fixed lattice parameters
    and returns the particle array, while also saving the array to a text file.

    Args:
        opts (object): Instance of the Synergia options class containing beam and macroparticle information

    Returns:
        bunch (ndarray): NumPy array with bunch with values {x, px, y, py, cdt, z, ID} for each particle

    """

    # Beam parameters
    gamma0 = 2.
    speciesMass = consts.m_p
    dgammaOgamma = 1.e-3
    #We want dpop no dE/E
    dpop = opts.dpop
    #Assume Gaussian longitudinal profile - put bunch length in m
    sigmaz = opts.stdz
    transverseEmittance = opts.emit_transverse
    numMacroParticles = opts.macro_particles

    t = 0.4 #fixed for IOTA 6-6 for now
    c = 0.01 #fixed for IOTA 6-6 for now

    #calculate beta at injection - center of NL element
    beta = 0.6538
    betaPrime = 0.0

    xOffset = 0. #m
    yOffset = 0. #m
    fileName = opts.bunch_file

    # Do not modify below this line
    E0 = gamma0 * consts.m_p * consts.c**2 * 6.24150934e9 #GeV/J
    ESpread = E0 * dgammaOgamma

    pArray = [0.]*numMacroParticles
    cdtArray = [0.]*numMacroParticles

    for idx in range(0,numMacroParticles):
        pArray[idx] = random.gauss(0, dpop)
        cdtArray[idx] = random.gauss(0, sigmaz)

    myBunchGenerator = EllipticBeam(t, c, beta, betaPrime)
    bunch = myBunchGenerator.generatebunch(transverseEmittance, numMacroParticles)

    for index,emit in enumerate(opts.emits):

        innerBunch = np.zeros(numMacroParticles) #bunch at this emittance
        transverseEmittance = emit

        myBunchGenerator = EllipticBeam(t, c, beta, betaPrime)
        #coords is an array of 4-vectors containing coordinate space information
        coords = myBunchGenerator.generatefixedbunch(transverseEmittance, numMacroParticles, opts.seed)


        for idx in range(len(coords)):
            pArray[idx] = random.gauss(0, dpop)
            cdtArray[idx] = random.gauss(0, sigmaz)

            #assign unique index value to each particle
            ID = index*len(coords) + idx

            coords[idx] = np.append(coords[idx],[cdtArray[idx],pArray[idx],ID])

        bunch.append(coords)

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


def make6Dellipticalbeam(lattice_file,lattice_name):

    '''Constructs a coasting elliptic beam with Gaussian longitduinal distribution and stores it in a Synergia bunch.

    A lattice file is required, as Synergia constructs the bunch relative to the reference particle as specified in mad-x.

    Arguments:
        lattice_file - A .madx file which stores the lattice
        lattice_name - The sequence name within the madx file

     '''

    #import lattice and construct options object
    #load lattice
    lattice = synergia.lattice.MadX_reader().get_lattice(lattice_name,lattice_file)
    length = lattice.get_length()
    ref = lattice.get_reference_particle() #reference particle
    ds = 0.01
    #nsteps = int(length/ds) +1 #calculate # of steps to take per turn
    #nsteps_per_element = nsteps/len(lattice.get_elements()) #not this isn't using future division, so returns int
    nsteps_per_element = 10 #hardcoded for IOTA nl elements
    nsteps = len(lattice.get_elements())*nsteps_per_element


    name = 'lattice_1IO_NL_Center'
    order = 1
    outputdir = 'order_'+str(order)+'_'+name

    opts = workflow.make_opts(name, order, outputdir, nsteps, nsteps_per_element)
    opts.macro_particles=1000
    opts.emitx = 7.0e-6
    workflow.make_path(outputdir)

    #make lattice chef propagating
    latticework.make_chef(lattice)

    stepper = synergia.simulation.Independent_stepper_elements(lattice, opts.map_order, opts.steps_per_element)
    lattice_simulator = stepper.get_lattice_simulator()
    #lattice_simulator = synergia.simulation.Lattice_simulator(lattice,opts.map_order)


    #construct bunch array and write it to a textfile
    ellipticalbeam6D(opts)

    #read in the file and construct a synergia bunch
    particles_file = 'myBunch.txt'
    comm = synergia.utils.Commxx(True)
    bucket_length = lattice_simulator.get_bucket_length()
    bucket_length = 0.05 #potential workaround
    myBunch = read_bunch.read_bunch(particles_file, ref, opts.real_particles, bucket_length, comm, verbose=False)

    return opts, lattice, lattice_simulator, stepper, myBunch
