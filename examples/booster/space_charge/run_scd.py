# booster_cell

import synergia
import numpy as np

lattice_file = """
beam,
  particle=proton,
  gamma=1.4425916541;

L03A: DRIFT,l=0.31505;
L03B: DRIFT,l=0.300581433056;
L03C: DRIFT,l=0.300581433056;
L03D: DRIFT,l=0.43067;
L03E: DRIFT,l=0.54551;
L03F: DRIFT,l=0.76926992944;
L03G: DRIFT,l=0.2;
L03H: DRIFT,l=0.064165;
L03I: DRIFT,l=0.064165;
MINS: DRIFT,l=0.5;
SA: DRIFT,l=0.176;
SB: DRIFT,l=0.256;
SC: DRIFT,l=0.6;
SEPT03: DRIFT,l=1.524;
HL03: HKICKER,l=0.0240152947613;
HS02: HKICKER,l=0.024;
MARK2D: MARKER;
MARK3: MARKER;
START_3: MARKER;
BPML03: MONITOR;
BPMLU3: MONITOR,l=0.0240152947613;
BPMS02: MONITOR,l=0.024;
QSLERR03: MULTIPOLE,ksl={0.0,0.0};
QSSERR02: MULTIPOLE,ksl={0.0,0.0};
QL03: QUADRUPOLE,l=0.0240152947613, k1=-0.005660797;
QS02: QUADRUPOLE,l=0.024,k1=0.00483088;
QSL03: QUADRUPOLE,l=0.0240152947613, k1s=0.000741168;
QSS02: QUADRUPOLE,l=0.024, k1s=0.024701692;
DMAGD03: SBEND,angle=0.0601575118477,k1=-0.05772919636,k2=-0.0404062293898,l=2.889176299;
DMAGU03: SBEND,angle=0.0601575118477,k1=-0.05772919636,k2=-0.0404062293898,l=2.889176299;
DOG03_1A_U: SBEND,angle=0.0178457763,e2=0.0356915526,l=0.24727249619,tilt=1.57079632679;
DOG03_2A_U: SBEND,angle=-0.0178457763,e1=0.0356915526,l=0.24727249619,tilt=1.57079632679;
DOG03_3A_U: SBEND,angle=-0.0178457763,e2=0.0356915526,l=0.24727249619,tilt=1.57079632679;
DOG03_4A_U: SBEND,angle=0.0178457763,e1=0.0356915526,l=0.24727249619,tilt=1.57079632679;
FMAGD03: SBEND,angle=0.0707421821963,k1=0.05442533685,k2=0.000681204960849,l=2.889009499;
FMAGU03: SBEND,angle=0.0707421821963,k1=0.05442533685,k2=0.000681204960849,l=2.889009499;
IBEX03: SBEND,l=0.33,tilt=1.57079632679;
SSL03: SEXTUPOLE,l=0.0240152947613, k2s=-0.03271327;
SSS02: SEXTUPOLE,l=0.024, k2s=0.032763832;
SXL03: SEXTUPOLE,l=0.0240152947613, k2=0.114899751;
SXS02: SEXTUPOLE,l=0.024, k2=0.114899751;
VL03: VKICKER,l=0.0240152947613;
VS02: VKICKER,l=0.024;
BOOSTER_CELL: LINE=(START_3,SA,HS02,VS02,QS02,BPMS02,QSS02,QSSERR02,SXS02,SSS02,SB,FMAGU03,MINS,DMAGU03,MARK2D,L03A,DOG03_1A_U,L03B,HL03,VL03,QL03,BPMLU3,QSL03,QSLERR03,SXL03,SSL03,L03C,DOG03_2A_U,L03D,SEPT03,L03E,DOG03_3A_U,L03F,DOG03_4A_U,L03G,IBEX03,L03H,BPML03,L03I,DMAGD03,MINS,FMAGD03,SC,MARK3);
BOOSTER: LINE=(BOOSTER_CELL);
"""

# Pickle helper is not necessary but is retained for this example
class Pickle_helper:
    __getstate_manages_dict__ = 1
    def __init__(self, *args):
        self.args = args
    def __getinitargs__(self):
        return self.args
    def __getstate__(self):
        return self.__dict__
    def __setstate__(self, state):
        self.__dict__ = state
        
class Ramp_actions(synergia.simulation.Propagate_actions, Pickle_helper):
    # The arguments to __init__ are what the Ramp_actions instance is
    # initialized with
    def __init__(self, multiplier, interval, dither_f):
        synergia.simulation.Propagate_actions.__init__(self)
        # pickling the arguments to the initializer allows the
        # module to resume after checkpointing.  They should be in the same
        # order as the arguments to __init__.
        Pickle_helper.__init__(self, multiplier)
        self.multiplier = multiplier
        self.interval = interval
        self.centroid = []
        self.element_list = ["ql03", "qs02"] #, "dmagu03", "dmagd03" ]
        self.old_angle_list = []
        self.old_k1_list = []
        self.old_k2_list = []
        self.injection_f = 37.77e6
        self.dither_f = dither_f
        self.turn_p = 1. / self.injection_f
        self.dither_p = 1. / self.dither_f
        
    def turn_end_action(self, stepper, bunch, turn_num):
        self.centroid.append(np.average(bunch.get_local_particles()[:, 0]))
        
        if turn_num % self.interval != 0:
            pass
        else:
            k = 0
            for element in stepper.get_lattice_simulator().get_lattice().get_elements():        
                if element.get_name() in self.element_list:
                    if turn_num == 0:
                        #self.old_angle_list.append(element.get_double_attribute("angle"))
                        self.old_k1_list.append(element.get_double_attribute("k1"))
                        #self.old_k2_list.append(element.get_double_attribute("k2"))
                        

                    else:
                        modifier = 1.5 * np.sin(2.0 * np.pi * turn_num / 1000.)
                        #element.set_double_attribute("angle", self.old_angle_list[k] * modifier)
                        element.set_double_attribute("k1", self.old_k1_list[k] * modifier)        
                        #element.set_double_attribute("k2", self.old_k2_list[k] * modifier)
                        k += 1

            stepper.get_lattice_simulator().update()
        if turn_num % 25 == 0:
            print("Turn:", turn_num )
        
    # other possible actions which could be present.
    # actions that are not present will default to internal stubs
    def first_action(self, stepper, bunch):
        pass
    def step_end_action(self, stepper, step, bunch, turn_num, step_num):
        pass
    def before_resume_action(self, stepper, bunch):
        pass

    

LATTICE_FILE = 'lattice.madx'
with open(LATTICE_FILE, 'w') as f:
    f.write(lattice_file)


comm = synergia.utils.Commxx(True)
#collective = synergia.simulation.Dummy_collective_operator('stub')
collective = synergia.collective.Space_charge_2d_open_hockney(
    comm,
    [16, 16, 1] # grid
)

lattice = synergia.lattice.MadX_reader().get_lattice('booster', LATTICE_FILE)


stepper = synergia.simulation.Split_operator_stepper_elements(
    lattice,
    1, # map_order
    collective,
    1 # num_steps
)

intensity = 4.5e12 / 84. / 100.
dither_f = 1440.0 * 10 * 10.

bunch = synergia.optics.generate_matched_bunch(
    stepper.get_lattice_simulator(), 
    arms=1.0e-3, 
    brms=1.0e-3, 
    crms=1.0e-4, 
    num_real_particles = int(intensity), 
    num_macro_particles=5000, 
    rms_index=[0,2,4], 
    seed=12345)



ramp_actions = Ramp_actions(1e-1, 1, dither_f)

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
bunch_simulator.add_per_step(synergia.bunch.Diagnostics_basic("basic_diagnostics_scd.h5"))


propagator = synergia.simulation.Propagator(stepper, ramp_actions)
propagator.set_checkpoint_period(1000)
propagator.set_checkpoint_with_xml(True)
propagator.propagate(bunch_simulator,
                     5000, 5000, 1)

#synergia.simulation.Propagator(stepper).propagate(
#    bunch_simulator, #dramp_actions,
#    100, # number of turns
#    0, # max_turns, Number of turns to run before writing checkpoint and stopping
#       # When max_turns is 0, the simulation continues until the end.
#    1, # verbosity
#)
