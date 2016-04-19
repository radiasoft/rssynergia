import sys

class Options:
    """
    A class designed to handle a fixed set of arguments shared across Synergia simulation and analysis tools.
    
    This iteration of the Options class draws on the Synergia Options class, and is meant to be broadly
    applicable for running Synergia simulations, and for providing relevant keyword arguments for processing
    diagnostic outputs and analysis tools for these simulations.
    
    The initializaiton of this class requires no arguments, and provides fixed initial values for a host
    of parameters.
    
    Initial parameter list:
        - hist, plots, inputfile, outputfile, show, hcoord, vcoord, bins. minh, maxh, minv, maxv, contour,
        num_countour, lattice_name, save, ID, path, norm, variance, relpath, lattice_simulator, num, scale
        turns, plot_lost, elliptic
        
    The class currently supports the creation of new attributes through the class method set_new(), which wraps
    Python's own setattr to provide fixed applicability towards class construction.
    
    
    """
    
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