# Jaeyoon Kim
# University of Michigan
# jaeykim@umich.edu

# This program take an initial wavefunction as input
# Then simulates the time evolution of probability density according to Schordinger's Equation
# -hBar^2 / 2m gradient(psi) = E * psi = i partial (psi) / partial (time)

import matplotlib.pyplot as plt
import numpy as np


# Example functions
# Should be in strings with x as the variable
# Function f should be infinitely differentiable
# and have values close to 0 at x = 0 and x = L
f1 = "( np.sin(3 * np.pi * x / self.L))"
f2 = "( np.e ** (-(x - 3) ** 2) )"
f3 = "( x * np.sin(5 * np.pi * x / self.L ) )"
f4 = "( np.sin( np.pi * x ** 3 / self.L ** 3 ) )"


# Global constants
h = 1                           # 6.626*10**-34 planck's constant
hBar = h / (2 * np.pi)

# Returns n^th element of orthonormal basis of solutions of the Schrodinger's Equation
def particle_in_box(n, L):
    return "(np.sqrt(2.0/" + str(L) + ") * np.sin(" + str(n) + "*x*np.pi/" + str(L) + "))"
# Can simulate all other potentials as long as there is a orthonormal basis function defined as
# basis(n,L) that returns a string of n^th orthonormal basis


class Wavefunction:
    def __init__(self, psi_in, basis_in, L_in, numPoints_in, maxN_in, mass_in):
        self.psi = psi_in
        self.L = float(L_in)
        self.numPoints = numPoints_in
        self.basis = basis_in
        self.maxN = maxN_in
        self.mass = mass_in

        self.xs = [self.L / self.numPoints * i for i in range(self.numPoints + 1)]
        self.normalize()
        self.components = self.getComponents()
        self.energyCoef = self.getEnergy()

    # Input function f (in strings with x as a variable) integrate from a to b
    # Uses Reimann sum limits
    def integrate(self, f):
        total=0.0
        dx = self.L / self.numPoints
        for i in range(self.numPoints):
            x = dx*i
            total = total + eval(f)*dx

        return total

    # Normalizes Psi
    def normalize(self):
        self.psi = str(1.0 / np.sqrt(self.integrate( self.psi + "*" + self.psi)) ) + "*" + self.psi

    # Returns the array of coefficients in the fourier series transform of psi
    def fourierTransform(self):
        coef = []
        for n in range(1,self.maxN+1):
                f = self.psi + "*"+ self.basis(n, self.L)
                coef.append(self.integrate(f))
        return coef

    # returns the array of energy of each quantum level
    def getEnergy(self):
            energy = []
            for n in range(1, self.maxN + 1):
                    energy.append(n ** 2 * h ** 2 / (8 * self.mass * self.L ** 2))
            return energy

    def getComponents(self):
        coef = self.fourierTransform()
        components = [coef[n - 1] * np.array([eval(self.basis(n, self.L)) for x in self.xs]) for n in range(1,self.maxN+1)]
        return components

        

    # Returns the probability density (psi^*psi) at time t
    def prob(self,t):
            prob = []
            realComp = np.zeros(self.numPoints + 1)
            imgComp = np.zeros(self.numPoints + 1)
            for n in range(self.maxN):
                    realComp = realComp + np.cos(hBar * self.energyCoef[n] * t) * np.array(self.components[n])
                    imgComp = imgComp + np.sin(hBar * self.energyCoef[n] * t) * np.array(self.components[n])
            prob = np.sqrt(np.square(realComp) + np.square(imgComp))
            return prob

    def simulate(self, t_start, t_end, scaleFactor):
        plt.ion()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        line, = ax.plot(psi.xs, self.prob(t_start), 'r-')
        for t in range(t_start, t_end):
            scaleFactor = 10
            line.set_ydata(psi.prob(t * scaleFactor))
            fig.canvas.draw()
        plt.close()        

psi = Wavefunction(f2, particle_in_box, 10, 500, 40, 1)
psi.simulate(0,500,10)


