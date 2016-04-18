import numpy as np
import random

class ellipticPostProcessing:

    def __init__(self, _t, _c, _beta, _betaPrime=0.):
        
        self.ellipticT = -1.*_t
        self.ellipticC = _c
        self.rtbeta    = np.sqrt(_beta)
        self.betaPrime = _betaPrime
        
    def computeHamiltonian(self, xHat, pxHat, yHat, pyHat):
        """Compute the Hamiltonian (1st invariant) for the integrable elliptic potential"""

        quadratic = 0.5 * (pxHat**2 + pyHat**2) + 0.5 * (xHat**2 + yHat**2)

        xN = xHat / self.ellipticC
        yN = yHat / self.ellipticC

        # Elliptic coordinates
        u = ( np.sqrt((xN + 1.)**2 + yN**2) +
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.
        v = ( np.sqrt((xN + 1.)**2 + yN**2) -
                  np.sqrt((xN - 1.)**2 + yN**2) )/2.

        f2u = u * np.sqrt(u**2 - 1.) * np.arccosh(u)
        g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))

        kfac = self.ellipticT * self.ellipticC**2
        elliptic = (f2u + g2v) / (u**2 - v**2)
        elliptic *= kfac

        hamiltonian = quadratic + elliptic
        return hamiltonian
        
        
    def computeinvariant(self, xHat, pxHat, yHat, pyHat):
        """Compute the 2nd invariant for the integrable elliptic potential"""

        xN = xHat / self.ellipticC
        yN = yHat / self.ellipticC

        # Elliptic coordinates
        u = (np.sqrt((xN + 1)**2 + yN**2) +\
                 np.sqrt((xN - 1)**2 + yN**2))/2.
        v = (np.sqrt((xN + 1)**2 + yN**2) -\
                 np.sqrt((xN - 1)**2 + yN**2))/2.

        # Elaborate x^2 + y^2
        f1u = self.ellipticC**2 * u**2 * (u**2 - 1)
        g1v = self.ellipticC**2 * v**2 * (1 - v**2)

        f2u = u * np.sqrt(u**2 - 1.) * np.arccosh(u)
        g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))

        kfac = self.ellipticT * self.ellipticC**2

        invariant = (x*py - y*px)**2 + self.ellipticC**2*(px**2) +\
            (2*self.ellipticC**2 ) * ((0.5*f1u - kfac*f2u)*v**2 \
            + (0.5 * g1v + kfac * g2v) * u**2 ) /(u**2 - v**2)

        return invariant
        
        
    def processparticles(self, x, px, y, py):
        """
        Given arrays of particles, compute the two invariants of the elliptic potential
        Args:
        x    canonical x position
        px   canonical x momentum
        y    canonical y position
        py   canonical y momentum
        
        Returns:
        H, I  The two invariants, with units of m-rad
        """
        
        xHat  = x/self.rtbeta
        pxHat = px*self.rtbeta - 0.5*self.betaPrime*xHat
        yHat  = y/self.rtbeta
        pyHat = py*self.rtbeta - 0.5*self.betaPrime*yHat
        
        H = self.computeHamiltonian(xHat, pxHat, yHat, pyHat)
        I = self.computeinvariant(xHat, pxHat, yHat, pyHat)
        
        return H, I