from __future__ import division # needed for 1/2 = 0.5
import math
from scipy import constants
from sympy import *

# assumes a large energy
def betaFromE(m_0,E):
    gammaSquared = E**2/m_0**2/constants.c**4+2*E/m_0/constants.c**2+1 # This will give the same matrices as TraceWin!!
    return sqrt(1-1/gammaSquared) # ^^
    #gammaSquared = (E/m_0/constants.c**2)**2 # wikipedia Energy-momentum_relation, which might be more correct
    #return sqrt(1-1/gammaSquared)

# assumes a large energy and not 100 % sure but I have calcs of this starting from the betaFromE and gammaFromBeta functions
def EFromBeta(m_0,beta):
    gamma = 1/sqrt(1-beta**2)
    return m_0*constants.c**2*(gamma-1)

def gammaFromBeta(beta):
    return 1/sqrt(1-beta**2)

def betaFromGamma(gamma):
    return sqrt(1-1/gamma**2)