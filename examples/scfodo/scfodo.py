# !/usr/bin/env synergia

import synergia

from synergia.foundation import Reference_particle
from synergia.lattice import Mad8_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_basic
from synergia.simulation import Bunch_simulator, Propagator
from mpi4py import MPI
import numpy as np
import sys, math

from scfodo_options import opts

# We wrap the entire simulation in a try..except block in order to allow
# for graceful failures under MPI.
try:
                                    
    # logger object only writes on one processor
    logger = synergia.utils.Logger(0)

    # a FODO lattice is used to define the initial bunch
    # read the lattice named "fodo" from the Mad8 file "fodo.lat"
    fodo_lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

    # Define a set of simulation steps for the FODO
    map_order = opts.map_order
    num_steps = opts.num_steps
    fodo_simulator = synergia.simulation.Lattice_simulator(fodo_lattice, map_order)

    # Define grid parameters
    gridx = 16
    gridy = 16
    gridz = 16
    grid = [gridx, gridy, gridz]
    ptcls_per_cell = opts.partpercell

    # Instantiate the bunch
    real_particles = opts.real_particles  # should be calculated from protons/m
    macro_particles = gridx * gridy * gridz * ptcls_per_cell

    x_emit = opts.x_emit  # m-rad, RMS
    y_emit = opts.y_emit  # m-rad, RMS
    slice_length = 0.01
    beta_x = 1.7
    beta_y = 1.7
    x_rms = math.sqrt(beta_x * x_emit)
    y_rms = math.sqrt(beta_y * y_emit)
    z_rms = 1.e6  # crazy large...?
    z_std = opts.z_std    # ad hoc choice...
    dp_std = 0.001

    periodic_bunch = opts.periodic
    if periodic_bunch == True:
        # next, populate
        bunch = synergia.optics.generate_matched_bunch(
                     fodo_simulator, x_rms, y_rms, z_rms,
                     real_particles, macro_particles,
                     rms_index=[0,2,4], seed=0,
                     bunch_index=0, comm=synergia.utils.Commxx(),
                     periodic=True)
    else:
        bunch = synergia.optics.generate_matched_bunch_transverse(
                     fodo_simulator, x_emit, y_emit, z_std, dp_std,
                     real_particles, macro_particles, seed=0)

    print >>logger, ' charge      = ', bunch.get_particle_charge()
    print >>logger, ' mass        = ', bunch.get_mass()
    print >>logger, ' num protons = ', bunch.get_real_num()
    print >>logger, ' num macro-p = ', bunch.get_total_num()

    # Create the space charge object
    use_simple = True
    if use_simple == False:
        # Trying for a meshed 3D slice, vey thin, with periodic longitudinal BCs
        space_charge = synergia.collective.Space_charge_2d_open_hockney(bunch.get_comm(), grid,
                                                                    True, False,
                                                                    slice_length,
                                                                    True, 0.)
    else:
        space_charge = synergia.collective.Space_charge_2d_open_hockney(bunch.get_comm(), grid
                                                                        )
    # mapping of the FODO lattice
    print_fodo_map = True
    if print_fodo_map == True:
        # Calculate the linear one-turn map for the FODO
        fodo_map = synergia.optics.one_turn_map.linear_one_turn_map(fodo_simulator)
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "Linear one turn map"
            print np.array2string(fodo_map,precision=6,suppress_small=True,max_line_width=200)

        l,v = np.linalg.eig(fodo_map)
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "eigenvalues of one turn map: ", l
            print "absolute values of eigenvalues (should all be 1): ", abs(l)
            print "fractional tunes from eigenvalues: ", np.log(l).imag/(2.0*np.pi)

    # Define a bunch simulator
    bunch_simulator = Bunch_simulator(bunch)

    # Create the split-operator stepping object
    fodo_stepper = synergia.simulation.Split_operator_stepper(fodo_simulator, space_charge, num_steps)

    # Number of turns
    turns = opts.turns  # a single pass through the line, since this isn't a ring
    max_turns = opts.max_turns  # Number of turns to run before writing checkpoint and stopping
                  # When max_turns is 0, the simulation continues until the end.
    verbosity = opts.verbosity  # Display information about each simulation step

    #
    # Define a set of bunch diagnostics
    #

    # diagnostics per step
    diagnostics = Diagnostics_basic("diag_basic_per_step.h5")
    bunch_simulator.add_per_step(diagnostics)

    # diagnostics per turn
    # bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("diag_full_per_turn.h5"))
    # bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("diag_ptcl_per_turn.h5"))
                    
    # Perform the simulation
    fodo_propagator = Propagator(fodo_stepper)
    fodo_propagator.propagate(bunch_simulator, turns, max_turns, verbosity)

except Exception, e:
    print ' '
    print 'An exception was thrown...'
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
