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
            print "ett varv till"
            multipart,envelope = elem.evaluate(multipart,envelope)

        return multipart,envelope

class Element:
    def __init__(self,name):
        self.name = name

    # abstract evaluate here so that the lattice containing all kinds of elements can be evaluated from the same function
    @abc.abstractmethod
    def evaluate(self,multipart,envelope):
        print "hej fran element"
        """Calculate the new particle data as the old data goes through"""
        return

class LinearElement(Element):
    def __init__(self, name):
        Element.__init__(self, name)

    def printInfo(self):
        return self.name

class Drift(LinearElement):
    def __init__(self, name, L):
        LinearElement.__init__(self, name)
        self.L = L
        self.M = self.createMatrixM(L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(L) # M should be a 9x9 matrix

    def createMatrixM(self,L):
        return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]
                ])
    def createMatrixT(self, L):
        return 0

    def evaluate(self,multipart,envelope):
        print "hej fran drift"
        newmultipart, newenvelope = self.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
        newmultipart, newenvelope = self.evaluateM(newmultipart,newenvelope) # use the new data for "normal" evaluation
        return newmultipart, newenvelope

    # Evaluate for Space Charges
    def evaluateSC(self,multipart,envelope):
        # disunite matrices
        n = 5
        Msp, Tsp = self.disunite(self.M,self.T,n)
        # evaluate the SC

        return 0, 0

    def disunite(self,M,T,n):
        L_n = self.L/n
        return 0, 0


    def evaluateM(self,multipart,envelope):
        # something M*multipart
        return 0, 0

class Quad(LinearElement):
    def __init__(self, name, K, L):
        LinearElement.__init__(self, name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

    """ taken from lie code where K = sqrt(k)"""
    def createMatrixM(self,K,L):
        if K > 0: # defocus in x
            return np.array([
                [cosh(K*L),sinh(K*L)/K,0,0,0,0],
                [K*sinh(K*L),cosh(K*L),0,0,0,0],
                [0,0,cos(K*L),sin(K*L)/K,0,0],
                [0,0,-K*sin(K*L),cos(K*L),0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]
                ])
        elif K < 0: # focus in x
            return np.array([
                [cos(K*L/2),sin(K*L/2)/K,0,0,0,0],
                [-K*sin(K*L/2),cos(K*L/2),0,0,0,0],
                [0,0,cosh(K*L/2),sinh(K*L/2)/K,0,0],
                [0,0,K*sinh(K*L/2),cosh(K*L/2),0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]
                ])
        else:
            return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]
                ])

    def createMatrixT(self, K, L):
        return 0

    def printInfo(self):
        return self.name + "\n" + str(self.M)

    def evaluate(self,multipart,envelope):
        print "hej fran quad"
        newmultipart, newenvelope = self.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
        newmultipart, newenvelope = self.evaluateM(newmultipart,newenvelope) # use the new data for "normal" evaluation
        return newmultipart, newenvelope

    # Evaluate for Space Charges
    def evaluateSC(self,multipart,envelope):
        # disunite matrices
        n = 5
        Msp, Tsp = self.disunite(self.M,self.T,n)
        # evaluate the SC

        return 0, 0

    def disunite(self,M,T,n):
        L_n = self.L/n
        return 0, 0


    def evaluateM(self,multipart,envelope):
        # something M*multipart
        return 0, 0































class NonLinearElement(Element):
    def __init__(self, name, ham, kval, lval, order):
        super(name, self).__init__()
        self.name = name
        self.ham = ham
        self.L = lval

        # Same this done 4 times, time for a new function?
        self.xfun = funFromHam(ham, order, qx)
        self.xprimefun = funFromHam(ham, order, px)
        self.yfun = funFromHam(ham, order, qy)
        self.yprimefun = funFromHam(ham, order, py)

        self.xf = self.xfun.subs([(k,kval),(l,lval)])
        self.xpf = self.xprimefun.subs([(k,kval), (l,lval)])
        self.yf = self.yfun.subs([(k,kval),(l,lval)])
        self.ypf = self.yprimefun.subs([(k,kval), (l,lval)])

        self.xf = lambdify((qx,px,qy,py),self.xf, "numpy")
        self.xpf = lambdify((qx,px,qy,py),self.xpf, "numpy")
        self.yf = lambdify((qx,px,qy,py),self.yf, "numpy")
        self.ypf = lambdify((qx,px,qy,py),self.ypf, "numpy")

    def printInfo(self):
        return self.name

    def evaluate(self, (mulxin, mulxpin, mulyin, mulypin)): # sending in xin and xpin as a vector (same for return) allows "recursive" calls and a lattice can be constructed
        xout = self.xf(mulxin, mulxpin, mulyin, mulypin)
        xpout = self.xpf(mulxin, mulxpin, mulyin, mulypin)
        yout = self.yf(mulxin, mulxpin, mulyin, mulypin)
        ypout = self.ypf(mulxin, mulxpin, mulyin, mulypin)

        return (xout,xpout,yout,ypout)

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