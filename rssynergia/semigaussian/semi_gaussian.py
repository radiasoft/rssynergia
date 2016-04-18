import numpy as np
import random
from scipy.optimize import newton

class SemiGaussianBeam:
    
    """ Generic class for generating a bunch distribution with certain properties. 
    Generates a numpy array for easy output/input into other codes.
    
    Args:
        beta (float) the beta function where the bunch is being matched, defaults to 1
        betaPrime (float) the derivative of the beta function, defaults to 0
    
    """
    
    def __init__(self, _beta=1, _betaPrime=0.):
        
        self.beta      = _beta
        self.betaPrime = _betaPrime

    def computeHamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (CS invariant) for the potential"""

        quadratic = 0.5*(pxHat**2 + pyHat**2) #+ 0.5 * (xHat**2 + yHat**2)

        hamiltonian = quadratic 
        return hamiltonian
        
    def computepotential(self, xHat, yHat):
        """Compute the general potential"""
        quadratic = 0.5*(xHat**2 + yHat**2)

        potential = quadratic
        return potential
        
    def whatsleft(self, yHat):
        """Return the difference btween the emittance and potential"""
        return self.emittance - self.computepotential(0, yHat)
    
        
    def generatefixedbunch(self, emittance, nParticles, seed):
        """ Generate a matched bunch with RMS emittance and number of particles
        Args:
        emittance (float): the RMS emittance of the bunch
        nParticles(int): the number of particles for the bunch
        seed (int): the random number generator seed for fixing particle coordinates
        
        Returns:
        bunch (list): a list of numpy arrays of 4D phase space, (x, px, y, py)
        """
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
         
        
        xMax = 2*np.sqrt(self.beta*emittance)
        #seed the random generator
        random.seed(seed)
        
        
        ##Uniformly distribute particles in coordinate space
        rx = 2.*(0.5-np.random.rand(nParticles*2))*xMax #the majority of these values are <1
        ry = 2.*(0.5-np.random.rand(nParticles*2))*xMax
        rad = np.square(rx/xMax) + np.square(ry/xMax)
        indices = np.where(rad <= 1)[0][:nParticles]
        x_vals = np.asarray(rx[indices])
        y_vals = np.asarray(ry[indices])
        
        
        ##Distribute particles as a Gaussian in momentum space - throw out those more than 4 sigma away from mean
        std = 2*emittance/xMax
        #Create more particles than needed to select from
        px_full = np.random.standard_normal(nParticles*2)*std
        py_full = np.random.standard_normal(nParticles*2)*std
        indices_px = np.where(np.abs(px_full) <= 3.*std)[0][:nParticles]
        indices_py = np.where(np.abs(py_full) <= 3.*std)[0][:nParticles]

        px = np.asarray(px_full[indices_px])
        py = np.asarray(py_full[indices_py])
        
        
        p_Array = np.transpose(np.vstack((x_vals,px,y_vals,py)))
        
        return p_Array