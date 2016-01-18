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

class LinearElement:
    def __init__(self, name):
        self.name = name

    def printInfo(self):
        return self.name

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

    def evaluatue(self,multipart,envelope):
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































class NonLinearElement:
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