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
f1 = "( np.sin(3 * np.pi * x / L))"
f2 = "( np.e ** (-(x - 3) ** 2) )"
f3 = "( x * np.sin(5 * np.pi * x / L ) )"
f4 = "( np.sin( np.pi * x ** 3 / L ** 3 ) )"

# Actual function that is to be evaluated
psi = f2

# Constant Declaration
L=10.0                          # length of the box
numBoxes = 1000                 # Number of points to plot
dx = L / numBoxes               # X distance between each plot
maxN = 70                       # Maximum quantum number to check
h = 1                           # 6.626*10**-34 planck's constant
hBar = h / (2 * np.pi)
m = 1                           # mass. for an electron, 9.109*10**-31


# Input function f (in strings with x as a variable) integrate from a to b
# Uses Reimann sum limits
def integrate(f):
    total=0.0
    for i in range(numBoxes):
        x = dx*i
        total = total + eval(f)*dx

    return total

# Returns n^th element of orthonormal basis of solutions of the Schrodinger's Equation
def wave(n):
    return "(np.sqrt(2.0/L) * np.sin(" + str(n) + "*x*np.pi/L))"

# Returns the normalization constant of psi
def normalConst(psi):
    return str(1.0 / np.sqrt(integrate( psi + "*" + psi)) )
        
# Normalize psi
psi = normalConst(psi) + "*"+ psi

# Returns the array of coefficients in the fourier series transform of psi
def fourierTransform(psi):
    coef = []
    for n in range(1,maxN+1):
            f = psi + "*"+ wave(n)
            coef.append(integrate(f))
    return coef

# returns the array of energy of each quantum level
def getEnergy():
        energy = []
        for n in range(1, maxN + 1):
                energy.append(n**2 * h ** 2 / (8*m*L**2))
        return energy

coef = fourierTransform(psi)
energy = getEnergy()

xs = [dx*i for i in range(numBoxes + 1)]

components = []
for n in range(1, maxN + 1):
        components.append( coef[n - 1] * np.sin( np.array(xs) * n * np.pi / L )  )

# Returns the probability density (psi^*psi) at time t
def prob(t):
        prob = []
        realComp = np.zeros(numBoxes + 1)
        imgComp = np.zeros(numBoxes + 1)
        for n in range(maxN):
                realComp = realComp + np.cos(hBar * energy[n] * t) * components[n]
                imgComp = imgComp + np.sin(hBar * energy[n] * t) * components[n]
        prob = np.sqrt(np.square(realComp) + np.square(imgComp) )
        return prob

psiSqrd = prob(0)

# Display
plt.ion()

fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(xs, prob(0), 'r-')
for t in range(5000):
    scaleFactor = 10
    line.set_ydata(prob(t * scaleFactor))
    fig.canvas.draw()
plt.close()
