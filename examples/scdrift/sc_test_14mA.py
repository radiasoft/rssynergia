#!/usr/bin/env python
#basic imports
import sys, os, time, random
import numpy as np
import scipy
from scipy.optimize import newton
from scipy import constants
import tables
from mpi4py import MPI

#Synergia imports
import synergia
import synergia_workflow

#diagnostics class imports
import read_bunch

#load options for SC_test
from sc_test_options import opts


#================ Some helpful scripts added to avoid using my homegrown modules ==================

def make_path(dirname):
    '''Create a directory with the specified name - avoid race conditions if possible'''

    try:
        os.makedirs(dirname)
    except OSError:
        if os.path.exists(dirname):
            # We are probably safe
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise
            
def cleanup(dirname):
    '''Cleanup files after run and move diagnostic outputs to proper directory.
    
    Arguments:
        -dirname: This is the relative path - e.g. full path = pwd + dirname
    
    '''
    curdir = os.getcwd()
    newdir = os.getcwd() + dirname
    
    for filename in os.listdir(curdir):
        if filename.endswith('.h5') or filename=='log':
            try:
                oldfn = '/'.join((curdir,filename))
                newfn = '/'.join((curdir,dirname,filename))
                #don't worry about shutil here, keep simple with os module
                os.rename(oldfn,newfn)
            except OSError:
                if os.path.exists(newfn):
                    #file already exists so delete and re-try
                    os.remove(newfn)
                    os.rename(oldfn,newfn)
                else:
                    #perhaps trying to move to a new disk or something that os can't handle
                    raise
            
def get_geometric_emittance(norm_emit,beta,gamma):
    '''Return the geometric emittance given a normalized emittance value and a beta/gamma.'''
    
    return norm_emit/(beta*gamma)

class standardBeam:
    
    ''' Generic class for generating a bunch distribution with certain properties. Generates a numpy array for easy output/input into other codes.'''
    
    def __init__(self, _beta=1, _betaPrime=0.):
        """ Generate a matched bunch for a fixed emittance
        Args:
        beta (float) the beta function where the bunch is being matched, defaults to 1
        betaPrime (float) the derivative of the beta function, defaults to 0
        """
        #self.ellipticT = -1.*_t
        #self.ellipticC = _c
        self.beta      = _beta
        self.betaPrime = _betaPrime

    def computeHamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (CS invariant) for the integrable potential"""

        quadratic = 0.5*(pxHat**2 + pyHat**2) #+ 0.5 * (xHat**2 + yHat**2)

        hamiltonian = quadratic 
        return hamiltonian
        
    def computePotential(self, xHat, yHat):
        quadratic = 0.5*(xHat**2 + yHat**2)

        potential = quadratic
        return potential
        
    def whatsLeft(self, yHat):
        return self.emittance - self.computePotential(0, yHat)

        
    def generateFixedBunch(self, emittance, nParticles, seed):
        """ Generate a matched bunch with single emittance and number of particles
        Args:
        emittance (float) the value of fixed H
        nParticles(int)   the number of particles for the bunch
        seed (int)        the random number generator seed for fixing particle coordinates
        
        Returns:
        bunch (list)  a list of numpy arrays of 4D phase space, (x, px, y, py)
        """
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
        
        # Use the lemming method to find the maximum y
        y0 = np.sqrt(emittance)
        #dy = 0.01*self.ellipticC
        #while self.computePotential(0, y0) < emittance:
        #    print self.computePotential(0,y0)
        #    y0 += dy
        
        #yMax = y0
        self.emittance = emittance
        yMax = newton(self.whatsLeft, y0)        
        
        # x is harder to bound due to the peanut nature of the potential -- estimate using conventional elliptic bunch
        xMax = yMax
        
        #seed the random generator
        random.seed(seed)
        
        # Generate particles by creating trials and finding particles with potential less than emittance, then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        while ptclsMade < nParticles:
            xTrial = 2.*(0.5 - random.random())*xMax
            yTrial = 2.*(0.5 - random.random())*yMax
            trialValue = self.computePotential(xTrial, yTrial)
            if trialValue < emittance:
                pMag = np.sqrt(2*(emittance - trialValue))
                pDir = 2*np.pi * random.random()
                pxHat = pMag * np.cos(pDir)
                pyHat = pMag * np.sin(pDir)
                xReal = xTrial * np.sqrt(self.beta)
                yReal = yTrial * np.sqrt(self.beta)
                #Changing these to minus signs - alpha should be -betaprime, and its x + alpha
                pxReal = (pxHat + 0.5*self.betaPrime*xTrial)/np.sqrt(self.beta)
                pyReal = (pyHat + 0.5*self.betaPrime*yTrial)/np.sqrt(self.beta)
                ptclCoords = np.array([xReal, pxReal, yReal, pyReal])
                phaseSpaceList.append(ptclCoords)
                ptclsMade += 1        
        
        return phaseSpaceList



def toyKVBeam6D(opts):
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
            myBunchGenerator = standardBeam(opts.betae) #use betae
        else:
            myBunchGenerator = standardBeam() #fixed beta=1, betaprime=0
        #coords is an array of 4-vectors containing coordinate space information
        coords = myBunchGenerator.generateFixedBunch(transverseEmittance, numMacroParticles, opts.seed)
        
        
        for idx in range(len(coords)):
            #EArray[idx] = random.gauss(E0, ESpread)
            pArray[idx] = random.gauss(0, dpop)
            #tArray[idx] = random.gauss(0., bunchLength)
            cdtArray[idx] = random.gauss(0, sigmaz)
            
            #assign unique index value to each particle
            ID = index*len(coords) + idx
        
            coords[idx] = np.append(coords[idx],[cdtArray[idx],pArray[idx],ID])
        
        #coords.append(pArray[idx])
        #coords2 = np.asarray(coords)
        #coords3 = coords2.flatten()
        
        
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


#================== Setting up logger and MPI comunicator ============================
#try:
#if True:
# this is the communicator object that will be used for MPI operations
comm = synergia.utils.Commxx()
myrank = comm.get_rank()
mpisize = comm.get_size()
verbose = opts.verbosity>0

logger = synergia.utils.Logger(0)


if myrank == 0:
    print "my rank is 0"
else:
    print "not rank 0"
    
    
#================ Lattice stuff =====================================
#Construct the lattice
#Construct the lattice
ol = 0.02 #2cm drift
steps_per_element = 2 #2 steps per drift
o = synergia.lattice.Lattice_element("drift", "o")
o.set_double_attribute("l", ol)

lattice = synergia.lattice.Lattice("test", synergia.lattice.Mad8_adaptor_map())
# Add copies of the lattice elements to the fodo lattice
lattice.append(o)

# Define reference particle to be a proton at 2.5 MeV
total_energy = synergia.foundation.pconstants.proton_mass + 2.5e-3 # 2.5 MeV protons
four_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.proton_mass, total_energy)
reference_particle = synergia.foundation.Reference_particle(synergia.foundation.pconstants.proton_charge,
                                        four_momentum)
opts.gamma = reference_particle.get_gamma()
opts.beta = reference_particle.get_beta()
lattice.set_reference_particle(reference_particle)

#64x64 2D grid with 1 longitudinal cell
gridx = 64
gridy = 64
gridz = 1
grid = [gridx, gridy, gridz]
opts.gridx = gridx
opts.gridy = gridy
opts.gridz = gridz

n_ppc = 100 #n_ppc particles per transverse cell
n_macro = n_ppc*opts.gridx*opts.gridy
opts.macro_particles = n_macro
outputdir = 'sc_drift_test'
opts.output_dir = outputdir
opts.relpath = opts.output_dir
make_path(outputdir)


opts.comm_divide = 8
if opts.comm_divide:
    sc_comm = synergia.utils.Commxx_divider(opts.comm_divide, False)
else:
    sc_comm = synergia.utils.Commxx(True)

#Use the 2D open Hockney solver
coll_operator = synergia.collective.Space_charge_2d_open_hockney(sc_comm, grid)

#generate stepper and lattice simulator
map_order = 1
nsteps_per_element = 2
opts.steps_per_element = nsteps_per_element
stepper = synergia.simulation.Split_operator_stepper_elements(lattice, map_order, coll_operator, opts.steps_per_element)
lattice_simulator = stepper.get_lattice_simulator()

opts.lattice = lattice
opts.lattice_simulator = lattice_simulator


#================================ Construct KV Beam =====================================

current = 14.e-3 #mA of current 
rp_perlength = current/(opts.beta*constants.c*constants.e)
bunch_length = 2*1.e-3 #effective bunch length 2 mm

opts.emit_n = 0.3*1.e-6 #We want 0.3 mm-mrad normalized emittance
opts.emits = [get_geometric_emittance(opts.emit_n,opts.beta,opts.gamma)] 
dpop = 0.0
opts.real_particles = rp_perlength*bunch_length
opts.betae = 0.5 #fix this
opts.alphae = 0.0
opts.macro_particles = n_macro

if myrank == 0:


    particles = toyKVBeam6D(opts)
    bunch = particles[0]
    bunch[:,4] = bunch_length*(np.random.random(len(bunch)) -0.5) #center at 0
    bunch[:,5] = opts.dpop*np.random.randn(1,len(bunch)) #set dp/p

    #Fix the weirdness with particle ID 4
    bunch[4] = bunch[100]
    bunch[4,6] = 4.0

    np.savetxt('myKVBunch.txt',bunch)         #write the bunch to a text file


bucket_length = bunch_length
particles_file = 'myKVBunch.txt'
myBunch = read_bunch.read_bunch(particles_file, reference_particle, opts.real_particles, bucket_length, comm)

# generated longitudinal coordinate is z position (beta*c*dt) but Synergia uses
# c*dt.  Divide by beta to get c*dt.
local_particles = myBunch.get_local_particles()
local_particles[:,4] /= opts.beta

if myrank ==0:
    print "Approximating a current of {}A in a bunch of length {} using {} particles".format(current,bunch_length, opts.real_particles)

#================================ Bunch Simulator =====================================
bunch_simulator = synergia.simulation.Bunch_simulator(myBunch)

#basic diagnostics - PER STEP
basicdiag = synergia.bunch.Diagnostics_basic("basic.h5", opts.output_dir)
bunch_simulator.add_per_step(basicdiag)

#include full diagnostics
fulldiag = synergia.bunch.Diagnostics_full2("full.h5", opts.output_dir)
bunch_simulator.add_per_turn(fulldiag)

#particle diagnostics - PER TURN
opts.turnsPerDiag = 1
particlediag = synergia.bunch.Diagnostics_particles("particles.h5",0,0,opts.output_dir)
bunch_simulator.add_per_turn(particlediag, opts.turnsPerDiag)

opts.turns = 200
opts.checkpointperiod = 50
opts.maxturns = opts.turns+1

print "setting up propagator for rank {}".format(myrank)

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(opts.checkpointperiod)
propagator.set_concurrent_io(opts.concurrent_io)

print "starting simulation for rank {}".format(myrank)
if myrank == 0:
    t_start = time.time()

propagator.propagate(bunch_simulator,opts.turns, opts.maxturns,opts.verbosity)

if myrank == 0:
    t_end = time.time()
    elapsed = t_end - t_start
    print "{} turns took {} seconds, or {} seconds per turn, with a grid size of {}x{}x{}.".format(opts.turns,elapsed, elapsed*1.0/opts.turns, gridx,gridy,gridz)

if myrank == 0:
    #clean up files
    cleanup(opts.output_dir)
