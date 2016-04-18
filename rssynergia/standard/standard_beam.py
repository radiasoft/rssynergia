import numpy as np
import random
from scipy.optimize import newton

class StandardBeam:
    
    ''' Generic class for generating a symmetric bunch distribution. 
    Generates a numpy array for easy output/input into other codes.
    
    Args:
        beta (float): the beta function where the bunch is being matched, defaults to 1
        betaPrime (float): the derivative of the beta function, defaults to 0
    
    '''
    
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
        """ Generate a matched KV bunch with single emittance and number of particles
        Args:
            emittance (float): the value of fixed H
            nParticles(int): the number of particles for the bunch
            seed (int): the random number generator seed for fixing particle coordinates
        
        Returns:
            bunch (list)  a list of numpy arrays of 4D phase space, (x, px, y, py)
        """
        
        # Generate some bounds on the transverse size to reduce waste in generating the bunch
        # Use the lemming method to find the maximum y
        y0 = np.sqrt(emittance)
        self.emittance = emittance
        
        yMax = newton(self.whatsleft, y0)        
        xMax = yMax
        
        #seed the random generator
        random.seed(seed)
        
        # Generate particles by creating trials and finding particles with potential less than emittance, then assign the rest to momentum
        ptclsMade = 0
        phaseSpaceList = []
        while ptclsMade < nParticles:
            xTrial = 2.*(0.5 - random.random())*xMax
            yTrial = 2.*(0.5 - random.random())*yMax
            trialValue = self.computepotential(xTrial, yTrial)
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