from __future__ import division # needed for 1/2 = 0.5
import math
from scipy import constants
from sympy import *

# assumes a large energy
def betaFromE(m_0,E):
    return sqrt(1-m_0**2*constants.c**4/E**2)