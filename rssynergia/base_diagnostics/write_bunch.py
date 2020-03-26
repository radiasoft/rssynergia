# -*- coding: utf-8 -*-
"""?

:copyright: Copyright (c) 2020 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function
import synergia
import numpy as np
import h5py

def write_bunch(bunch,  file_name = 'particles.h5' , reference_particle = None, txt = False):
    '''
    Write a Synergia bunch to file, mimicing the style of the bunch diagnostics. Defaults to .h5 output.

    Assumes that the main processor has access to all particles.

    Unlike the standard bunch_diagnostics, this method outputs the beam attributes as h5 attributes rather than as datasets.

    Arguments:
        - bunch (synergia.bunch.bunch.Bunch): A Synergia bunch object to be written to file
        - fn (Optional[String]): File name for saving the bunch - defaults to 'particles.h5'
        - txt (Optional[Bool]): Whether to write to a .txt file - defaults to False

    '''

    particles = bunch.get_local_particles()
    num_particles = particles.shape[0]

    #Define attributes for the bunch - this requires a reference particle or other context
    #check to see if reference particle is passed
    if reference_particle is not None:
        charge = reference_particle.get_charge()
        mass = reference_particle.get_total_energy()/reference_particle.get_gamma()
        pz = reference_particle.get_momentum()

    else:
        #pass defaults
        charge = 1
        mass = 0.9382723128
        pz = 10.

    #specify these as 0 because particle distribution is uncoupled from simulation
    rep = 0
    s_n = 0
    tlen = 0


    if txt:
        #txt file write is straightforward numpy save
        np.savetxt(file_name,particles)

    else:
        #write to HDF-5 by default

        #Create h5 file
        dump_file = h5py.File(file_name, 'w')

        #specify attributes
        dump_file.attrs['charge'] = charge #Particle charge in units of proton charge
        dump_file.attrs['mass'] = mass #Particle equivalent mass in GeV
        dump_file.attrs['pz'] = pz #Reference particle momentum in GeV/c

        dump_file.attrs['s_n'] = s_n #Current s value (between 0 and lattice length), defaults to 0
        dump_file.attrs['tlen'] = tlen #Total length traversed during simulation, defaults to 0
        dump_file.attrs['rep'] = rep #Current repetition, defaults to 0

        #write particle phase space array as a single dataset
        ptcl_data = dump_file.create_dataset('particles', (num_particles,7))
        ptcl_data[:] = particles

        dump_file.close()
