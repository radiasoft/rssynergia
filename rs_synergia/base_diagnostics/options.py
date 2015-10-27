import sys

#Options class definition - All other diagnostics will make use of this class or subclasses

#To Do: Make some subclasses with more specific options for different plots

class Options:
    def __init__(self):
        self.hist = False
        self.plots = 1
        self.inputfile = None
        self.outputfile = None
        self.show = True
        self.hcoord = None
        self.vcoord = None
        self.bins = 20
        self.minh = -sys.float_info.max
        self.maxh = sys.float_info.max
        self.minv = -sys.float_info.max
        self.maxv = sys.float_info.max
        self.contour = None
        self.num_contour = None
        self.lattice_name = None
        self.save = False #default don't save
        self.ID = None
        self.path = None
        self.norm = False
        self.variance = None
        self.relpath = None
        self.lattice_simulator = None
        self.num = None #number of turns to plot
        self.scale = None #figure size scaling
        self.turns = None #specific turns to plot for certain types of plots
        self.plot_lost = False #specify to plot lost particles
        self.elliptic = False #specifies flag for plotting invariants for elliptic potential - defaults to false
        
    def set_new(self,name,val):
        '''A class method which wraps setattr()'''
        return setattr(self,name,val)