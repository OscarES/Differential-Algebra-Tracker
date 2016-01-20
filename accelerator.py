from __future__ import division # needed for 1/2 = 0.5
import numpy as np
from scipy.misc import *
#from scipy.linalg import *
from scipy import linalg
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import math
import random
import matplotlib.pyplot as plt
import abc

class Lattice:
    def __init__(self,name):
        self.name = name
        self.lattice = list()

    def appendElement(self, element):
        self.lattice.append(element)

    def evaluate(self, multipart,envelope):
        for elem in self.lattice:
            #print "ett varv till"
            multipart,envelope = elem.evaluate(multipart,envelope)

        return multipart,envelope

class Element:
    def __init__(self,name, linear):
        self.name = name
        self.linear = linear

    # abstract evaluate here so that the lattice containing all kinds of elements can be evaluated from the same function
    @abc.abstractmethod
    def evaluate(self,multipart,envelope):
        """Calculate the new particle data as the old data goes through"""
        return

    @abc.abstractmethod
    def printInfo(self):
        """print the name and maybe some info"""
        return

    def isLinear(self):
        return self.linear

class LinearElement(Element):
    def __init__(self, name):
        Element.__init__(self, name, 1)

    def printInfo(self):
        return self.name

### DRIFT !!!!
class Drift(LinearElement):
    def __init__(self, name, L):
        LinearElement.__init__(self, name)
        self.L = L
        self.M = self.createMatrixM(L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(L) # M should be a 9x9 matrix

        # disunite the matrices, CLEAN UP THIS CODE!!!!!! (should disunite be removed or kept?)
        self.n = 5
        self.Lsp = self.L/self.n
        self.Msp, self.Tsp = self.disunite(self.M,self.T,self.n)

        # and the space charge matrices for the disunited matrices (SC always implies a split matrix)
        self.Msc, self.Tsc = self.spaceChargeMatrix()

    def createMatrixM(self,L):
        return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,L],
                [0,0,0,0,0,1]
                ])
    def createMatrixT(self, L):
        return 0

    def disunite(self,M,T,n):
        Msp = self.createMatrixM(self.Lsp)
        return Msp, 0

    #Maybe the space charge should be its own class????????
    def spaceChargeMatrix(self):
        deltas = self.Lsp

        # beam data is needed as input to calculate the following variables...
        beta = 0.9 # beam data needed
        gamma = 1/sqrt(1-beta**2)
        m = 1.0 # should be proton mass but would rather have mass and q as input so that other particles can be used
        q = 1.0 # -II-
        c = 299792458.0 # in metric (metric for everything perhaps?) 

        # <.> is called "norm_of_."
        # <x^2> (norm_of_xx), <xx'> (norm_of_xx') and <x'^2> (norm_of_xpxp) come from (119) in ref 1.
        norm_of_xx = 1.0
        # <x> (norm_of_x) come from !!!
        norm_of_x = 0.9
        # <xE_x> (norm_of_xE_x) come from !!!
        norm_of_xE_x = 1.0

        f_scx = gamma**3*beta**2*m*c**2/q*(norm_of_xx - norm_of_x**2)/norm_of_xE_x # shoud really be eqn (152) in ref 1.
        f_scy = 1.0 # shoud really be eqn (152) (with x->y) in ref 1.
        f_scz = 1.0 # shoud really be eqn (152) (with x->z) in ref 1.
        xbar = 1.0
        ybar = 1.0
        zbar = 1.0

        Msc = np.array([
                [1,0,0,0,0,0,0],
                [deltas/f_scx,1,0,0,0,0,-xbar*deltas/f_scx],
                [0,0,1,0,0,0,0],
                [0,0,deltas/f_scy,1,0,0,-ybar*deltas/f_scy],
                [0,0,0,0,1,0,0],
                [0,0,0,0,deltas/f_scz,1,-zbar*deltas/f_scz],
                [0,0,0,0,0,0,1]
            ])
        return Msc, 0

    def evaluate(self,multipart,envelope):
        # some for loop that goes through all of the disunited parts
        #newmultipart, newenvelope = multipart, envelope
        for i in range(0,self.n):

        #print "hej fran drift"
        #print "evaluate multipart: " + str(multipart)

            # The should be a check if Msc and Tsc need to be updated if the beam properties have changed a lot!!!!!!!!
            #if beamChanged(beam):
                #self.Msc, self.Tsc = spaceChargeMatrix()

            multipart, envelope = self.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
        return multipart, envelope

    # Evaluate for Space Charges, se ref A.
    def evaluateSC(self,multipart,envelope):
        # evaluate the SC
        # for n
        #     newdata = Msp*olddata
        #     olddata = newdata

        SC = 1
        if not SC:
            return multipart, envelope # just returns the same values as of now!
        else:
            for j in range(0,len(np.atleast_1d(multipart))):
                extendedphasespace = np.vstack((multipart[j][0][0:6], 1)) # adds the dispersion 1 term
                extendedphasespace = np.dot(self.Msc, extendedphasespace) # here calculations are made
                reducedphasespace = extendedphasespace[0:6] # throws away the dispersion 1 term
                multipart[j] = np.array([reducedphasespace, multipart[j][1]]) # s remains the same because the particles don't go anywhere. They "go" in evaluateM()
            return multipart,envelope

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        return multipart, 0


### QUAD
class Quad(LinearElement):
    def __init__(self, name, K, L):
        LinearElement.__init__(self, name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix
        
        # disunite matrices
        self.n = 5
        self.Lsp = self.L/self.n
        self.Msp, self.Tsp = self.disunite(self.M,self.T,self.n)

    """ taken from lie code where K = sqrt(k)"""
    def createMatrixM(self,K,L):
        if K > 0: # defocus in x
            return np.array([
                [cosh(K*L),sinh(K*L)/K,0,0,0,0],
                [K*sinh(K*L),cosh(K*L),0,0,0,0],
                [0,0,cos(K*L),sin(K*L)/K,0,0],
                [0,0,-K*sin(K*L),cos(K*L),0,0],
                [0,0,0,0,1,L],
                [0,0,0,0,0,1]
                ])
        elif K < 0: # focus in x
            return np.array([
                [cos(K*L),sin(K*L)/K,0,0,0,0],
                [-K*sin(K*L),cos(K*L),0,0,0,0],
                [0,0,cosh(K*L),sinh(K*L)/K,0,0],
                [0,0,K*sinh(K*L),cosh(K*L),0,0],
                [0,0,0,0,1,L],
                [0,0,0,0,0,1]
                ])
        else:
            return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,L],
                [0,0,0,0,0,1]
                ])

    def createMatrixT(self, K, L):
        return 0

    def disunite(self,M,T,n):
        Msp = self.createMatrixM(self.K, self.Lsp)
        return Msp, 0

    def printInfo(self):
        return self.name + "\n" + str(self.M)

    def evaluate(self,multipart,envelope):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            multipart, envelope = self.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
        return multipart, envelope

    # Evaluate for Space Charges
    def evaluateSC(self,multipart,envelope):
        # evaluate the SC

        SC = 0
        if not SC:
            return multipart, envelope # just returns the same values as of now!
        else:
            return 0,0

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        return multipart, 0

### SPACE CHARGE!!!!!
class SpaceCharge(LinearElement):
    def __init__(self, name, deltas,beam):
        LinearElement.__init__(self, name)
        self.deltas = deltas
        self.beam = beam

        self.Msc = spaceChargeMatrix()

    def updateBeam(self, beam):
        self.beam = beam
        return

    def beamChanged(self, newbeam):
        threshold = 0.1
        if diffBtwBeams(self.beam,newbeam) > threshold:
            return 1
        else:
            return 0

    def diffBtwBeams(self, beam1,beam2):
        diff = 0
#        for bla
#            diff += diff_each_variable
        return diff

    def spaceChargeMatrix(self, beam):
        deltas = self.Lsp

        # beam data is needed as input to calculate the following variables...
        beta = 0.9 # beam data needed
        gamma = 1/sqrt(1-beta**2)
        m = 1.0 # should be proton mass but would rather have mass and q as input so that other particles can be used
        q = 1.0 # -II-
        c = 299792458.0 # in metric (metric for everything perhaps?) 

        # <.> is called "norm_of_."
        # <x^2> (norm_of_xx), <xx'> (norm_of_xx') and <x'^2> (norm_of_xpxp) come from (119) in ref 1.
        norm_of_xx = 1.0
        # <x> (norm_of_x) come from !!!
        norm_of_x = 0.9
        # <xE_x> (norm_of_xE_x) come from !!!
        norm_of_xE_x = 1.0

        f_scx = gamma**3*beta**2*m*c**2/q*(norm_of_xx - norm_of_x**2)/norm_of_xE_x # shoud really be eqn (152) in ref 1.
        f_scy = 1.0 # shoud really be eqn (152) (with x->y) in ref 1.
        f_scz = 1.0 # shoud really be eqn (152) (with x->z) in ref 1.
        xbar = 1.0
        ybar = 1.0
        zbar = 1.0

        Msc = np.array([
                [1,0,0,0,0,0,0],
                [deltas/f_scx,1,0,0,0,0,-xbar*deltas/f_scx],
                [0,0,1,0,0,0,0],
                [0,0,deltas/f_scy,1,0,0,-ybar*deltas/f_scy],
                [0,0,0,0,1,0,0],
                [0,0,0,0,deltas/f_scz,1,-zbar*deltas/f_scz],
                [0,0,0,0,0,0,1]
            ])
        return Msc

    def evaluateSC(self,multipart,envelope):
        for j in range(0,len(np.atleast_1d(multipart))):
            # The should be a check if Msc and Tsc need to be updated if the beam properties have changed a lot!!!!!!!!
            #if beamChanged(envelope):
                #self.Msc, self.Tsc = spaceChargeMatrix(envlope)

            extendedphasespace = np.vstack((multipart[j][0][0:6], 1)) # adds the dispersion 1 term
            extendedphasespace = np.dot(self.Msc, extendedphasespace) # here calculations are made
            reducedphasespace = extendedphasespace[0:6] # throws away the dispersion 1 term
            multipart[j] = np.array([reducedphasespace, multipart[j][1]]) # s remains the same because the particles don't go anywhere. They "go" in evaluateM()
        return multipart,envelope


























class NonLinearElement(Element):
    def __init__(self, name):
        Element.__init__(self, name, 0)

    def printInfo(self):
        return self.name

# class NonLinearElement(Element):
#     def __init__(self, name, ham, kval, lval, order):
#         Element.__init__(self, name, 1)
#         super(name, self).__init__()
#         self.name = name
#         self.ham = ham
#         self.L = lval

#         # Same this done 4 times, time for a new function?
#         self.xfun = funFromHam(ham, order, qx)
#         self.xprimefun = funFromHam(ham, order, px)
#         self.yfun = funFromHam(ham, order, qy)
#         self.yprimefun = funFromHam(ham, order, py)

#         self.xf = self.xfun.subs([(k,kval),(l,lval)])
#         self.xpf = self.xprimefun.subs([(k,kval), (l,lval)])
#         self.yf = self.yfun.subs([(k,kval),(l,lval)])
#         self.ypf = self.yprimefun.subs([(k,kval), (l,lval)])

#         self.xf = lambdify((qx,px,qy,py),self.xf, "numpy")
#         self.xpf = lambdify((qx,px,qy,py),self.xpf, "numpy")
#         self.yf = lambdify((qx,px,qy,py),self.yf, "numpy")
#         self.ypf = lambdify((qx,px,qy,py),self.ypf, "numpy")

#     def printInfo(self):
#         return self.name

#     def evaluate(self, (mulxin, mulxpin, mulyin, mulypin)): # sending in xin and xpin as a vector (same for return) allows "recursive" calls and a lattice can be constructed
#         xout = self.xf(mulxin, mulxpin, mulyin, mulypin)
#         xpout = self.xpf(mulxin, mulxpin, mulyin, mulypin)
#         yout = self.yf(mulxin, mulxpin, mulyin, mulypin)
#         ypout = self.ypf(mulxin, mulxpin, mulyin, mulypin)

#         return (xout,xpout,yout,ypout)

# just a copy of Quad so far. TODO: use the scraps from the lie code!
class Sextu(NonLinearElement):
    def __init__(self, name, K, L, M):
        NonLinearElement.__init__(self, name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

    """ taken from lie code where K = sqrt(k)"""


    def printInfo(self):
        return self.name + "\n" + str(0)

    def evaluate(self,multipart,envelope):
        newmultipart, newenvelope = evaluateSC(multipart,envelope) # evaluate the SC
        newmultipart, newenvelope = evaluateM(newmultipart,newenvelope) # use the new data for "normal" evaluation
        return newmultipart, newenvelope

    # Evaluate for Space Charges
    def evaluateSC(self,multipart,envelope):
        # disunite matrices
        Msp, Tsp = disunite(M,T,n)
        # evaluate the SC

        return 0, 0

    def disunite(self,M,T,n):
        L_n = self.L/n
        return 0, 0

    def evaluateM(self,multipart,envelope):
        # something M*multipart
        return 0, 0

# just a copy of Quad so far. TODO: use the scraps from the lie code!
class Multi(NonLinearElement):
    def __init__(self, name, K, L, M):
        NonLinearElement.__init__(self, name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

# just a copy of Quad so far. TODO: use the scraps from the lie code!
class Cavity(NonLinearElement):
    def __init__(self, name, K, L, M):
        NonLinearElement.__init__(self, name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

# references
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition
# A. 7.2. Space Charge Impulses in simulatingbeamswithellipsoidalsymmetry-secondedition