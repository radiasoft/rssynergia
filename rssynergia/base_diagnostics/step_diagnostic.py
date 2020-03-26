# -*- coding: utf-8 -*-
"""?

:copyright: Copyright (c) 2020 RadiaSoft LLC.  All Rights Reserved.
:license: http://www.apache.org/licenses/LICENSE-2.0.html
"""
from __future__ import absolute_import, division, print_function
import synergia
import os
import numpy as np
import h5py as h5
from mpi4py import MPI

comm_world = MPI.COMM_WORLD

try:
    import __builtin__
except ImportError:
    # Python 3
    import builtins as __builtin__

my_rank = comm_world.rank

def print(*args, **kwargs):
    """Overload print to prevent all ranks from printing"""
    if my_rank == 0:
        return __builtin__.print(*args, **kwargs)

class CustomDiagnostic(synergia.simulation.Propagate_actions):
    def __init__(self, stepper, element_names=[], step_numbers=[], positions=[]):
        """
        Create a step_end_action that will called by Synergia to output bunch data at arbitrary points in the lattice.
        These points can be chosen by element name, positions in the lattice, or direct input of step number.

        Each chosen diagnostic point should be accompanied by a callable function that will operate on the bunch object
        and return diagnostic data.

        Data is accumulated in Datum objects with information on where and when it was collected.

        :param stepper (synergia stepper object): Stepper being used in propagator.
        :param element_names (list of (str, func)): Optional list of tuples comprising element names and corresponding
            diagnostic function to call
        :param step_numbers (list of (int, func)): Optional list of tuples comprising step number and corresponding
            diagnostic function to call
        :param positions (list of (float, func)): Optional list of tuples comprising position (in meters) and
            corresponding diagnostic function to call
        """
        synergia.simulation.Propagate_actions.__init__(self)

        self.stepper = stepper
        self.steps = []
        self.diagnostics = []
        self.data = []

        for elem in element_names:
            elem_steps = self.find_element_steps(self.stepper, elem[0])
            if len(elem_steps) > 0:
                self.steps.append(elem_steps[len(elem_steps) // 2])  # output in middle of element
                self.diagnostics.append(elem[1])
            elif len(elem_steps) == 1:
                self.steps.append(elem_steps[0])  # output in middle of element
                self.diagnostics.append(elem[1])
            else:
                print("Could not find element: {}".format(elem[0]))
        for step in step_numbers:
            if step[0] not in self.steps:
                self.steps.append(step[0])
                self.diagnostics.append(step[1])
        for pos in positions:
            pos_step = self.find_step_position(self.stepper, pos[0])
            print("For position: {}, the closest step is {} m away".format(pos[0], pos_step[1]))
            if pos_step[0] not in self.steps:
                self.steps.append(pos_step[0])
                self.diagnostics.append(pos[1])

        for step, diag in zip(self.steps, self.diagnostics):
            assert callable(diag), "Diagnostic {} is not callable".format(diag)

            x = Datum(*self.find_step_information(self.stepper, step))
            x.diagnostic_name = diag.__name__
            self.data.append(x)


    @staticmethod
    def find_step_information(stepper, step_number):
        for i, step in enumerate(stepper.get_steps()):
            if i == step_number:
                oper = step.get_operators()[-1]
                slc = oper.get_slices()[-1]
                slice_element = slc.get_lattice_element()
                position = stepper.get_lattice_simulator().get_lattice_functions(slc).arc_length

                return i, slice_element.get_name(), position

    @staticmethod
    def find_element_steps(stepper, element_name):
        # TODO: Need to check element memory address to split out multiple element instances
        steps = []
        for step_number, step in enumerate(stepper.get_steps()):
            oper = step.get_operators()[-1]
            slc = oper.get_slices()[-1]
            slice_element = slc.get_lattice_element()
            if slice_element.get_name() == element_name:
                position = stepper.get_lattice_simulator().get_lattice_functions(slc).arc_length
                steps.append(step_number)
                print("Step: {}: Slice {}, Element {}, position {}".format(step_number, slc,
                                                                           slice_element.get_name(),
                                                                           position))

        return steps

    @staticmethod
    def find_step_position(stepper, target_position):
        closest_step = (0, 1e130)
        for step_number, step in enumerate(stepper.get_steps()):
            oper = step.get_operators()[-1]
            slc = oper.get_slices()[-1]
            position = stepper.get_lattice_simulator().get_lattice_functions(slc).arc_length
            if abs(target_position - position) < closest_step[1]:
                closest_step = (step_number, abs(target_position - position))

        return closest_step

    def write_datafiles(self, directory, filesnames=None):
        if filesnames:
            assert len(filesnames) == len(self.data), 'Number of supplied filenames not equal to number of datasets'
        else:
            filesnames = [None] * len(self.data)

        for data, name in zip(self.data, filesnames):
            data.write_data(directory, filename=name)

    def step_end_action(self, stepper, step, bunch, turn_num, step_num):
        """
        Overloads Synergia's default synergia.simulation.Propagate_actions.step_end_action method.
        Must maintain parameter input order (stepper, step, bunch, turn_num, step_num) to function.
        :param stepper: Synergia stepper object
        :param step: Individual step object
        :param bunch: Bunch object
        :param turn_num: (int) Current turn number
        :param step_num: (int) Current step on this turn
        :return:
        """
        # TODO: Add option for particular turn intervals
        # TODO: Collect particle data from other processors automatically if statistical data?
        if step_num in self.steps:
            indices = np.where(np.array(self.steps) == step_num)[0]
            for index in indices:
                datum = self.diagnostics[index](bunch)
                self.data[index].turn_num.append(turn_num)
                self.data[index].data.append(datum)


class Datum:
    def __init__(self, step_num, element_name, position):
        self.data = []
        self.turn_num = []
        self.step_num = step_num
        self.element_name = element_name
        self.position = position
        self.diagnostic_name = None

    def write_data(self, directory, filename=None):
        # TODO: allow for data caching on an interval
        if my_rank == 0:
            if not os.path.exists(directory):
                os.makedirs(directory)
            if not filename:
                filename = 'data_{}_step_{}.h5'.format(self.diagnostic_name, self.step_num)

            filepath = os.path.join(directory, filename)

            datafile = h5.File(filepath, 'w')
            for name, attr in [('step_num', self.step_num), ('element_name', self.element_name),
                               ('position', self.position)]:
                datafile.attrs[name] = attr

            datafile.create_dataset('turns', data=self.turn_num)
            if type(self.data[0]) == np.ndarray:
                for i, d in enumerate(self.data):
                    datafile.create_dataset('data{}'.format(i), data=d)
            else:
                datafile.create_dataset('data', data=np.array(self.data))

            datafile.close()

def xdistribution(bunch):
    all_particles = comm_world.gather(bunch.get_local_particles(), root=0)
    if comm_world.rank == 0:
        all_particles = np.vstack(all_particles)[:, 0]

        minx, maxx = all_particles.min(), all_particles.max()
        hist, bins = np.histogram(all_particles, range=(minx, maxx), bins='fd')
        centers = []
        for i in range(1, bins.size):
            centers.append((bins[i] + bins[i - 1]) / 2.)
        return np.array([hist, centers])

def xydistribution(bunch):
    all_particles = comm_world.gather(bunch.get_local_particles(), root=0)
    if comm_world.rank == 0:
        all_particles = np.vstack(all_particles)

        minx, maxx = all_particles[:, 0].min(), all_particles[:, 0].max()
        miny, maxy = all_particles[:, 1].min(), all_particles[:, 1].max()
        hist, binsx, binsy = np.histogram2d(all_particles[:, 0],
                                    all_particles[:, 1],
                                    range=[[minx, maxx], [miny, maxy]], bins=64)
        hist = np.append(hist, [-1])
        return np.array([hist, binsx, binsy])
