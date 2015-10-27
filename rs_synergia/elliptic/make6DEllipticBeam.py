import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow

import numpy as np
import matplotlib.pyplot as plt

from base_diagnostics import read_bunch
from base_diagnostics import workflow
from elliptic import EllipticBeam6D


#import lattice and construct options object
#load lattice
lattice = synergia.lattice.MadX_reader().get_lattice("iota", "/Users/ncook/Synergia_Tests/lattices/Iota6-6/lattice_2IO_nll.madx")
length = lattice.get_length()
ref = lattice.get_reference_particle() #reference particle
ds = 0.01
nsteps = int(length/ds) +1 #calculate # of steps to take per turn
nsteps_per_element = nsteps/len(lattice.get_elements()) #not this isn't using future division, so returns int

name = 'lattice_2IO_NL'
order = 1
outputdir = 'order_'+str(order)+'_'+name

opts = workflow.make_opts(name, order, outputdir, nsteps, nsteps_per_element)
opts.macro_particles=1000
opts.emitx = 1.0e-5
#workflow.make_path(outputdir)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice,opts.map_order)

#construct bunch array and write it to a textfile
EllipticBeam6D.EllipticalBeam6D(opts)

#read in the file and construct a synergia bunch
particles_file = 'myBunch.txt'
comm = synergia.utils.Commxx(True)
bucket_length = lattice_simulator.get_bucket_length()
myBunch = read_bunch.read_bunch(particles_file, ref, opts.real_particles, bucket_length, comm, verbose=False)
