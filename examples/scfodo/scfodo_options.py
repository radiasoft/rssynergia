#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("map_order", 1, "map order")
opts.add("num_steps", 20, "number of steps through the lattice")
opts.add("gridx", 16, "x grid size")
opts.add("gridy", 16, "y grid size")
opts.add("gridz", 1,  "z grid size")
opts.add("partpercell", 8, "number of macro particles per cell in the bunch")
opts.add("real_particles", 1.0e8, "total charge (in units of e) in the bunch")

opts.add("periodic", True, "generate periodic bunches")

opts.add("x_emit", 2.0e-6, "x RMS emittance [m-rad]")
opts.add("y_emit", 2.0e-6, "y RMS emittance [m-rad]")
opts.add("z_std", 1.0, "z RMS length [m]")
opts.add("dpop", 1.0e-3, "(delta p)/p")
opts.add("seed", 1415926, "random number seed; 0 for automatic calculation (GSL)")
opts.add("turns", 4, "number of times to track through fodo lattice")
opts.add("max_turns", 0, "maximum number of turns to run before checkpointing and stopping; 0 to not stop")
opts.add("verbosity", 2, "verbosity level of simulation")

# Create the job manager for the simulation fodo.py, including the above options.
# When creating job directories, include the file fodo.lat.
job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat"])
