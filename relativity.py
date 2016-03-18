from __future__ import division # needed for 1/2 = 0.5
import math
from scipy import constants
from sympy import *

# assumes a large energy
def betaFromE(m_0,E):
    gammaSquared = E**2/m_0**2/constants.c**4+2*E/m_0/constants.c**2+1
    return sqrt(1-1/gammaSquared)

def gammaFromBeta(beta):
    return 1/sqrt(1-beta**2)