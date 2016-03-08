from __future__ import division # needed for 1/2 = 0.5
import numpy as np
from numpy import inf
from scipy.misc import *
#from scipy.linalg import *
from scipy import linalg
from scipy.integrate import quad
from scipy import constants, linspace
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import math
import random
import matplotlib.pyplot as plt
import abc
from particleFactory import envelopeFromMultipart
from relativity import betaFromE

class Lattice:
    def __init__(self,name,beamdata,twiss,multipart):
        self.name = name
        self.lattice = list()

        ## Beam properties
        self.beamdata = beamdata
        self.twiss = twiss
        self.multipart = multipart

    def appendElement(self, element):
        self.lattice.append(element)
        return 1

    def printLattice(self):
        text = ""
        for elem in self.lattice:
            text = text + elem.printInfo() + "\n"
        return text

    def evaluate(self, multipart,envelope,twiss):
        for elem in self.lattice:
            multipart,envelope, twiss = elem.evaluate(multipart,envelope,twiss)
            #print "twiss: " + str(twiss)
        return multipart,envelope, twiss

    def relativityAtTheEnd(self, multipart,envelope):
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
    def __init__(self, name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        LinearElement.__init__(self, "drift " + name)
        self.L = L
        self.M = self.createMatrixM(L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix

        # disunite the matrices, CLEAN UP THIS CODE!!!!!! (should disunite be removed or kept?)
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        self.Msp, self.Tsp = self.disunite(self.M,self.T,self.n)
        #print "self.Tsp: \n" + str(self.Tsp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('drift_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, multipart, twiss, beamdata)

    def printInfo(self):
        return self.name + "\t L: " + str(self.L)

    def createMatrixM(self,L):
        return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,L],
                [0,0,0,0,0,1]
                ])
    def createMatrixT(self, M):
        # T is size 9x9 and defined from eqn (1.56) in ref C
        return np.array([
                [M[0,0]**2, 2*M[0,0]*M[0,1], M[0,1]**2, 0, 0, 0, 0, 0, 0],
                [M[0,0]*M[1,0], M[0,0]*M[1,1]+M[0,1]*M[1,0], M[0,1]*M[1,1], 0, 0, 0, 0, 0, 0],
                [M[1,0]**2, 2*M[1,0]*M[1,1], M[1,1]**2, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, M[2,2]**2, 2*M[2,2]*M[2,3], M[2,3]**2, 0, 0, 0],
                [0, 0, 0, M[2,2]*M[3,2], M[2,2]*M[3,3]+M[2,3]*M[3,2], M[2,3]*M[3,3], 0, 0, 0],
                [0, 0, 0, M[3,2]**2, 2*M[3,2]*M[3,3], M[3,3]**2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, M[4,4]**2, 2*M[4,4]*M[4,5], M[4,5]**2],
                [0, 0, 0, 0, 0, 0, M[4,4]*M[5,4], M[4,4]*M[5,5]+M[4,5]*M[5,4], M[4,5]*M[5,5]],
                [0, 0, 0, 0, 0, 0, M[5,4]**2, 2*M[5,4]*M[5,5], M[5,5]**2]
                ])

    def disunite(self,M,T,n):
        Msp = self.createMatrixM(self.Lsp)
        Tsp = self.createMatrixT(Msp)
        return Msp, Tsp

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            if self.spaceChargeOn:
                self.sc.updateMatrix(multipart,twiss)
                multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC
                #print "envelope[0]: " + str(envelope[0])
                twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
                twiss[4] = envelope[3] / twiss[5]
                twiss[7] = envelope[6] / twiss[8]
            multipart, envelope = self.evaluateMT(multipart,envelope) # use the new data for "normal" evaluation
            
        return multipart, envelope, twiss

    def evaluateMT(self,multipart,envelope):
        # should just go through a disunited part
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        #envelope = np.dot(self.Tsp, envelope)
        #print "Envelope from multipart instead!"
        envelope = envelopeFromMultipart(multipart)
        return multipart, envelope

### DIPOLE
class Dipole(LinearElement):
    def __init__(self, name, rho, alpha, n, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        LinearElement.__init__(self, "dipole " + name)
        self.rho = rho
        self.alpha = alpha
        self.nparam = n
        self.beta = beamdata[0]
        self.gamma = 1/sqrt(1-self.beta**2)
        self.L = rho*alpha
        self.K_x = sqrt(1-self.nparam)/rho
        self.K_y = self.nparam/rho
        self.M = self.createMatrixM(self.rho, self.L, self.K_x, self.K_y, self.beta, self.gamma) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix
        
        # disunite matrices
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        self.Msp, self.Tsp = self.disunite(self.M,self.T,self.n)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('quad_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, multipart, twiss, beamdata)

    """ taken from lie code where K = sqrt(k)"""
    def createMatrixM(self, rho, L, K_x, K_y, beta, gamma):
        return np.array([
            [cos(K_x*L),sin(K_x*L)/K_x,0,0,0,(1-cos(K_x*L))/rho/K_x**2],
            [-K_x*sin(K_x*L),cos(K_x*L),0,0,0,sin(K_x*L)/rho/K_x],
            [0,0,cos(K_y*L),sin(K_y*L)/K_y,0,0],
            [0,0,-K_y*sin(K_y*L),cos(K_y*L),0,0],
            [-sin(K_x*L)/rho/K_x,-(1-cos(K_x*L))/rho/K_x/L**2,0,0,1,-(K_x*L*beta**2-sin(K_x*L))/(rho**2*K_x**3)+(L/gamma**2)*(1-1/rho**2/K_x**2)],
            [0,0,0,0,0,1]
            ])

    def createMatrixT(self, M):
        # T is size 9x9 and defined from eqn (1.56) in ref C
        return np.array([
                [M[0,0]**2, 2*M[0,0]*M[0,1], M[0,1]**2, 0, 0, 0, 0, 0, 0],
                [M[0,0]*M[1,0], M[0,0]*M[1,1]+M[0,1]*M[1,0], M[0,1]*M[1,1], 0, 0, 0, 0, 0, 0],
                [M[1,0]**2, 2*M[1,0]*M[1,1], M[1,1]**2, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, M[2,2]**2, 2*M[2,2]*M[2,3], M[2,3]**2, 0, 0, 0],
                [0, 0, 0, M[2,2]*M[3,2], M[2,2]*M[3,3]+M[2,3]*M[3,2], M[2,3]*M[3,3], 0, 0, 0],
                [0, 0, 0, M[3,2]**2, 2*M[3,2]*M[3,3], M[3,3]**2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, M[4,4]**2, 2*M[4,4]*M[4,5], M[4,5]**2],
                [0, 0, 0, 0, 0, 0, M[4,4]*M[5,4], M[4,4]*M[5,5]+M[4,5]*M[5,4], M[4,5]*M[5,5]],
                [0, 0, 0, 0, 0, 0, M[5,4]**2, 2*M[5,4]*M[5,5], M[5,5]**2]
                ])

    def disunite(self,M,T,n):
        Msp = self.createMatrixM(self.rho, self.Lsp, self.K_x, self.K_y, self.beta, self.gamma)
        Tsp = self.createMatrixT(Msp)
        return Msp, Tsp

    def printInfo(self):
        return self.name + "\t rho: " + str(self.rho) + "\t L: " + str(self.L) + "\t K_x: " + str(self.K_x) + "\t K_y: " + str(self.K_y) + "\t beta: " + str(self.beta) + "\t Alpha: " + str(self.alpha) + "\t nparam: " + str(self.nparam)

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            if self.spaceChargeOn:
                self.sc.updateMatrix(multipart,twiss)
                multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        return multipart, envelope

### QUAD
class Quad(LinearElement):
    def __init__(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        LinearElement.__init__(self, "quad " + name)
        #self.name = name
        self.K = K
        self.L = L
        self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix
        
        # disunite matrices
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        self.Msp, self.Tsp = self.disunite(self.M,self.T,self.n)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('quad_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, multipart, twiss, beamdata)

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

    def createMatrixT(self, M):
        # T is size 9x9 and defined from eqn (1.56) in ref C
        return np.array([
                [M[0,0]**2, 2*M[0,0]*M[0,1], M[0,1]**2, 0, 0, 0, 0, 0, 0],
                [M[0,0]*M[1,0], M[0,0]*M[1,1]+M[0,1]*M[1,0], M[0,1]*M[1,1], 0, 0, 0, 0, 0, 0],
                [M[1,0]**2, 2*M[1,0]*M[1,1], M[1,1]**2, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, M[2,2]**2, 2*M[2,2]*M[2,3], M[2,3]**2, 0, 0, 0],
                [0, 0, 0, M[2,2]*M[3,2], M[2,2]*M[3,3]+M[2,3]*M[3,2], M[2,3]*M[3,3], 0, 0, 0],
                [0, 0, 0, M[3,2]**2, 2*M[3,2]*M[3,3], M[3,3]**2, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, M[4,4]**2, 2*M[4,4]*M[4,5], M[4,5]**2],
                [0, 0, 0, 0, 0, 0, M[4,4]*M[5,4], M[4,4]*M[5,5]+M[4,5]*M[5,4], M[4,5]*M[5,5]],
                [0, 0, 0, 0, 0, 0, M[5,4]**2, 2*M[5,4]*M[5,5], M[5,5]**2]
                ])

    def disunite(self,M,T,n):
        Msp = self.createMatrixM(self.K, self.Lsp)
        Tsp = self.createMatrixT(Msp)
        return Msp, Tsp

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K)

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            if self.spaceChargeOn:
                self.sc.updateMatrix(multipart,twiss)
                multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        return multipart, envelope

### SPACE CHARGE!!!!! C. Allen's approach
class SpaceCharge(LinearElement):
    def __init__(self, name, deltas, multipart, twiss, beamdata):
        LinearElement.__init__(self, name)
        self.deltas = deltas
        self.multipart = multipart
        self.twiss = twiss
        self.beamdata = beamdata

        self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)

    #def updateBeam(self, twiss):
    #    self.twiss = twiss
    #    return

    def beamChanged(self, newtwiss):
        threshold = 0.1
        if diffBtwBeams(self.twiss, newtwiss) > threshold:
            return 1
        else:
            return 0

    def diffBtwBeams(self, twiss1,twiss2):
        diff = 0
#        for bla
#            diff += diff_each_variable
        return diff

    def updateMatrix(self,multipart,twiss):
        self.twiss = twiss
        self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)
        return 1

    def R_D(self, x, y, z):
        # from (110) in ref 1.
        #print "hello1"
        #print "x: " + str(x)
        #print "y: " + str(y)
        #print "z: " + str(z)
        result = quad(lambda t : 3/2*1/(sqrt(t+x) * sqrt(t+y) * (t+z)**(3/2)), 0, inf)
        #print "result: " + str(result)
        #print "hello2"
        # result[0] is the result and result[1] is the error
        return result[0]

    def spaceChargeMatrix(self, multipart, twiss, beamdata):
        # beamdata: beta (speed), mass, charge, lambda (RF-wavelength)
        beta = beamdata[0]
        rf_lambda = beamdata[1]
        m = beamdata[2]
        q = beamdata[3]
        
        # beam data is needed as input to calculate the following variables...
        N = len(multipart) # this info should come from the multipart (len(multipart))
        Q = q*N # from (19) in ref 1.
        gamma = 1/sqrt(1-beta**2)
        c = constants.c # in metric (metric for everything perhaps?)
        vac_perm = constants.epsilon_0

        I = N*q*c/rf_lambda # from ref. E

        ## Courant-Snyder or Twiss params
        # envelope comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
        alpha_x = twiss[0]
        beta_x = twiss[1]
        epsilon_rms_x = twiss[2]
        alpha_y = twiss[3]
        beta_y = twiss[4]
        epsilon_rms_y = twiss[5]
        alpha_z = twiss[6]
        beta_z = twiss[7]
        epsilon_rms_z = twiss[8]

        ## envelope X, Xp, Y, Yp, Z and Zp
        X = sqrt(5*beta_x*epsilon_rms_x)
        Xp = -alpha_x*sqrt(5*epsilon_rms_x/beta_x)

        Y = sqrt(5*beta_y*epsilon_rms_y)
        Yp = -alpha_y*sqrt(5*epsilon_rms_y/beta_y)

        Z = sqrt(5*beta_z*epsilon_rms_z)
        #print "beta_z: " + str(beta_z)
        Zp = -alpha_z*sqrt(5*epsilon_rms_z/beta_z)

        # <.> is called "norm_of_."
        # <x^2> (norm_of_xx), <xx'> (norm_of_xx') and <x'^2> (norm_of_xpxp) come from (119) in ref 1.
        norm_of_xx = 1/5*X**2
        # <x> (norm_of_x) come from (81) in ref 1.
        norm_of_x = 0.0
        # <xx'> (norm_of_xxp) come from (119) in ref 1. NOT used here but here as a reminder
        norm_of_xxp = 1/5*X*Xp
        # <x'x'> (norm_of_xpxp) come from (119) in ref 1. NOT used here but here as a reminder
        norm_of_xpxp = 1/5*Xp**2+5*epsilon_rms_x**2/X**2

        norm_of_yy = 1/5*Y**2
        norm_of_y = 0.0
        norm_of_yp = 1/5*Y*Yp
        norm_of_ypyp = 1/5*Yp**2+5*epsilon_rms_y**2/Y**2

        norm_of_zz = 1/5*Z**2
        norm_of_z = 0.0
        norm_of_zp = 1/5*Z*Zp
        norm_of_zpzp = 1/5*Zp**2+5*epsilon_rms_z**2/Z**2

        # <xE_x> (norm_of_xE_x) come from (115) in ref 1.
        norm_of_xE_x = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_xx*self.R_D(norm_of_yy, norm_of_zz, norm_of_xx)
        norm_of_yE_y = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_yy*self.R_D(norm_of_xx, norm_of_zz, norm_of_yy)
        norm_of_zE_z = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_zz*self.R_D(norm_of_xx, norm_of_yy, norm_of_zz)
        # R_D_bracket() = R_D_() from (111) and (112) in ref 1.

        # eqn (152) in ref 1. But maybe it should be (157) in ref 1. instead??????????
        f_scx = gamma**3*beta**2*m*c**2/q*(norm_of_xx - norm_of_x**2)/norm_of_xE_x # eqn (152) in ref 1. But maybe it should be (157) in ref 1. instead??????????
        f_scy = gamma**3*beta**2*m*c**2/q*(norm_of_yy - norm_of_y**2)/norm_of_yE_y # (152) (with x->y) in ref 1.
        f_scz = gamma**3*beta**2*m*c**2/q*(norm_of_zz - norm_of_z**2)/norm_of_zE_z # (152) (with x->z) in ref 1.        

        # Mean of x,y and z from all the particles
        xbar = sum([multipart[i][0][0] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting x out of each particle and then taking the mean
        ybar = sum([multipart[i][0][2] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting y out of each particle and then taking the mean
        zbar = sum([multipart[i][0][4] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting z out of each particle and then taking the mean

        # Matrix eqn (154) in ref 1.
        Msc = np.array([
                [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
                [self.deltas/f_scx,1.0,0.0,0.0,0.0,0.0,-xbar*self.deltas/f_scx],
                [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,self.deltas/f_scy,1.0,0.0,0.0,-ybar*self.deltas/f_scy],
                [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,self.deltas/f_scz,1.0,-zbar*self.deltas/f_scz],
                [0.0,0.0,0.0,0.0,0.0,0.0,1.0]
            ])
        return Msc

    def evaluateSC(self,multipart,envelope):
        for j in range(0,len(np.atleast_1d(multipart))):
            # The should be a check if Msc and Tsc need to be updated if the beam properties have changed a lot!!!!!!!!
            #if beamChanged(envelope):
                #self.Msc, self.Tsc = spaceChargeMatrix(envlope)

            extendedphasespace = np.append(multipart[j][0][0:6], 1) # adds the dispersion 1 term
            extendedphasespace = np.dot(self.Msc, extendedphasespace) # here calculations are made
            reducedphasespace = extendedphasespace[0:6] # throws away the dispersion 1 term
            multipart[j] = np.array([reducedphasespace, multipart[j][1]]) # s remains the same because the particles don't go anywhere. They "go" in evaluateM()
        envelope = envelopeFromMultipart(multipart) # the envelope is just calculated from the particles (NOT ON ITS OWN)
        return multipart,envelope






### SPACE CHARGE!!!!! Elliptical integral approach
class SpaceChargeEllipticalIntegral(LinearElement):
    def __init__(self, name, deltas, multipart, twiss, beamdata):
        LinearElement.__init__(self, name)
        self.deltas = deltas
        self.multipart = multipart
        self.twiss = twiss
        self.beamdata = beamdata

        self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)

    #def updateBeam(self, twiss):
    #    self.twiss = twiss
    #    return

    def beamChanged(self, newtwiss):
        threshold = 0.1
        if diffBtwBeams(self.twiss, newtwiss) > threshold:
            return 1
        else:
            return 0

    def diffBtwBeams(self, twiss1,twiss2):
        diff = 0
#        for bla
#            diff += diff_each_variable
        return diff

    def updateMatrix(self,multipart,twiss):
        self.twiss = twiss
        self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)
        return 1

    def R_D(self, x, y, z):
        # from (110) in ref 1.
        result = quad(lambda t : 3/2*1/(sqrt(t+x) * sqrt(t+y) * (t+z)**(3/2)), 0, inf)
        # result[0] is the result and result[1] is the error
        return result[0]

    def spaceChargeMatrix(self, multipart, twiss, beamdata):        
        #print "\n\n\nNew element!"
        #print "deltas: " + str(self.deltas)
        # beamdata: beta (speed), lambda (RF-wavelength), mass, charge
        beta = beamdata[0]
        #print "beta: " + str(beta)
        rf_lambda = beamdata[1]
        #print "rf_lambda: " + str(rf_lambda)
        m = beamdata[2]
        q = beamdata[3]
        
        # beam data is needed as input to calculate the following variables...
        N = len(multipart) # this info should come from the multipart (len(multipart))
        Q = q*N # from (19) in ref 1.
        gamma = 1/sqrt(1-beta**2)
        #print "gamma: " + str(gamma)
        c = constants.c # in metric (metric for everything perhaps?)
        vac_perm = constants.epsilon_0

        I = N*q*c/rf_lambda # from ref. E #I = 0.065 # beam data from ref F
        #print "I: " + str(I)

        ## Courant-Snyder or Twiss params
        # envelope comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
        alpha_x = twiss[0]
        beta_x = twiss[1]
        #print "beta_x: " + str(beta_x)
        epsilon_rms_x = twiss[2]
        #print "epsilon_rms_x: " + str(epsilon_rms_x)
        alpha_y = twiss[3]
        beta_y = twiss[4]
        epsilon_rms_y = twiss[5]
        alpha_z = twiss[6]
        beta_z = twiss[7]
        epsilon_rms_z = twiss[8]

        ## envelope X, Xp, Y, Yp, Z and Zp
        X = sqrt(5*beta_x*epsilon_rms_x)
        #print "X: " + str(X)
        Xp = -alpha_x*sqrt(5*epsilon_rms_x/beta_x)
        #print "Xp: " + str(Xp)

        Y = sqrt(5*beta_y*epsilon_rms_y)
        #print "Y: " + str(Y)
        Yp = -alpha_y*sqrt(5*epsilon_rms_y/beta_y)
        #print "Yp: " + str(Yp)

        Z = sqrt(5*beta_z*epsilon_rms_z)
        #print "Z: " + str(Z)
        Zp = -alpha_z*sqrt(5*epsilon_rms_z/beta_z)
        #print "Zp: " + str(Zp)

        ## Deviations from multipart (x,y,z). So there should really be one matrix per particle. Now I just use one particle's data
        #r_x = multipart[0][0][0]
        r_x = X
        #print "r_x: " + str(r_x)
        #r_y = multipart[0][0][2]
        r_y = Y
        #print "r_y: " + str(r_y)
        #r_z = multipart[0][0][4]
        r_z = Z
        #print "r_z: " + str(r_z)
        
        #r_x = 1 # semi-axes of the ellipsiod. What are they really? putting them to one works
        #r_y = 1
        #r_z = 1

        #x = multipart[0][0][0]
        #print "x: " + str(x)
        #y = multipart[0][0][2]
        #print "y: " + str(y)
        #z = multipart[0][0][4]
        #print "z: " + str(z)

        ## Eliptical integral?
        #s = Z/sqrt(X*Y)
        #epsilon_of_s_integral = s*quad(lambda t : 1/(sqrt(t+1)* (t+s**2)**(3/2)), 0, inf)[0] # eqn 5 in ref B

        # Eliptical integral
        g = gamma*r_z/sqrt(r_x*r_y) # eqn 40 from ref E.
        #print "g: " + str(g)
        f_of_g_integral = g/2*quad(lambda t : 1/((t+1)*(t+g**2)**(3/2)), 0, inf)[0] # eqn 41 from ref E.
        #print "f_of_g_integral: " + str(f_of_g_integral)

        #G_x = 3*(1-f_of_g_integral)*x/(r_x*(r_x+r_y)*r_z) # eqn 36 from ref E.
        #print "G_x: " + str(G_x)
        #G_y = 3*(1-f_of_g_integral)*y/(r_y*(r_x+r_y)*r_z) # eqn 37 from ref E.
        #print "G_y: " + str(G_y)
        #G_z = 3*f_of_g_integral*z/(r_x*r_y*r_z) # eqn 38 from ref E.
        #print "G_z: " + str(G_z)

        G_x_without_x = 3*(1-f_of_g_integral)/(r_x*(r_x+r_y)*r_z) # saving x for eval (x will come from multipart)
        G_y_without_y = 3*(1-f_of_g_integral)/(r_y*(r_x+r_y)*r_z)
        G_z_without_z = 3*f_of_g_integral/(r_x*r_y*r_z)

        #U_scx = I*rf_lambda*G_x/(4*math.pi*constants.epsilon_0*c*gamma**2) # eqn 33 from ref E.
        #print "U_scx: " + str(U_scx)
        #U_scy = I*rf_lambda*G_y/(4*math.pi*constants.epsilon_0*c*gamma**2) # eqn 34 from ref E.
        #print "U_scy: " + str(U_scy)
        #U_scz = I*rf_lambda*G_z/(4*math.pi*constants.epsilon_0*c) # eqn 35 from ref E.
        #print "U_scz: " + str(U_scz)

        U_scx_without_x = I*rf_lambda*G_x_without_x/(4*math.pi*constants.epsilon_0*c*gamma**2)
        U_scy_without_y = I*rf_lambda*G_y_without_y/(4*math.pi*constants.epsilon_0*c*gamma**2)
        U_scz_without_z = I*rf_lambda*G_z_without_z/(4*math.pi*constants.epsilon_0*c)

        #delta_P_x = q*U_scx*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_x: " + str(delta_P_x)
        #delta_P_y = q*U_scy*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_y: " + str(delta_P_y)
        #delta_P_z = q*U_scz*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_z: " + str(delta_P_z)

        delta_P_x_without_x = q*U_scx_without_x*self.deltas/(m*c**2*beta)
        delta_P_y_without_y = q*U_scy_without_y*self.deltas/(m*c**2*beta)
        delta_P_z_without_z = q*U_scz_without_z*self.deltas/(m*c**2*beta)

        # Converting from delta_P to delta_xp
        #v = beta*c
        #p = gamma*m*v
        ##print "p: " + str(p)
        #delta_xp = delta_P_x/p # xp = p_x/p. eqn 150 and 151 from ref 1.
        #print "delta_xp: " + str(delta_xp)
        #delta_yp = delta_P_y/p # yp = p_y/p. eqn 150 and 151 from ref 1.
        #print "delta_yp: " + str(delta_yp)
        #delta_zp = delta_P_z/p # zp = p_z/p. eqn 150 and 151 from ref 1.
        #print "delta_zp: " + str(delta_zp)

        delta_xp_without_x = delta_P_x_without_x#/p # #/p means that delta_P_i is actually delta_ip
        delta_yp_without_y = delta_P_y_without_y#/p
        delta_zp_without_z = delta_P_z_without_z#/p

        #Msc = np.array([
        #        [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
        #        [0.0,1.0,0.0,0.0,0.0,0.0,delta_xp],
        #        [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
        #        [0.0,0.0,0.0,1.0,0.0,0.0,delta_yp],
        #        [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
        #        [0.0,0.0,0.0,0.0,0.0,1.0,delta_zp],
        #        [0.0,0.0,0.0,0.0,0.0,0.0,1.0]
        #    ])

        Msc_new_from_without = np.array([
                [1.0,0.0,0.0,0.0,0.0,0.0],
                [delta_xp_without_x,1.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,1.0,0.0,0.0,0.0],
                [0.0,0.0,delta_yp_without_y,1.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,1.0,0.0],
                [0.0,0.0,0.0,0.0,delta_zp_without_z,1.0]
            ])

        return Msc_new_from_without
        #return Msc

    def evaluateSC(self,multipart,envelope):
        for j in range(0,len(np.atleast_1d(multipart))):
            # The should be a check if Msc and Tsc need to be updated if the beam properties have changed a lot!!!!!!!!
            #if beamChanged(envelope):
                #self.Msc, self.Tsc = spaceChargeMatrix(envlope)

            #zz
            #extendedphasespace = np.append(multipart[j][0][0:6], 1) # adds the dispersion 1 term
            #extendedphasespace = np.dot(self.Msc, extendedphasespace) # here calculations are made
            #reducedphasespace = extendedphasespace[0:6] # throws away the dispersion 1 term
            #multipart[j] = np.array([reducedphasespace, multipart[j][1]]) # s remains the same because the particles don't go anywhere. They "go" in evaluateM()
            #zz

            ### New way of calculating (from the without params)
            multipart[j][0][0:6] = np.dot(self.Msc, multipart[j][0][0:6]) # self.Msc was set to the new matrix

        envelope = envelopeFromMultipart(multipart) # the envelope is just calculated from the particles (NOT ON ITS OWN)
        return multipart,envelope























## Here all the heavy Lie calculations will be made
class LieAlgebra():
    def __init__(self):
        self.qx = Symbol('qx')
        self.qy = Symbol('qy')
        self.qz = Symbol('qz')
        self.px = Symbol('px')
        self.py = Symbol('py')
        self.pz = Symbol('pz')

        self.l = Symbol('l') # arbitrary length
        self.k = Symbol('k')        

    def lieop(self, f,g):
        dfdqx = f.diff(self.qx) 
        dfdpx = f.diff(self.px)
        dgdqx = g.diff(self.qx)
        dgdpx = g.diff(self.px)
        sumtermx = dfdqx*dgdpx-dfdpx*dgdqx
    
        dfdqy = f.diff(self.qy) 
        dfdpy = f.diff(self.py)
        dgdqy = g.diff(self.qy)
        dgdpy = g.diff(self.py)
        sumtermy = dfdqy*dgdpy-dfdpy*dgdqy
    
        dfdqz = f.diff(self.qz) 
        dfdpz = f.diff(self.pz)
        dgdqz = g.diff(self.qz)
        dgdpz = g.diff(self.pz)
        sumtermz = dfdqz*dgdpz-dfdpz*dgdqz
    
        colfcolg = sumtermx + sumtermy + sumtermz
    
        return colfcolg

    def lietransform(self, ham, vof0, order):#,t):, #zzz
        voft = vof0
        for i in range(1,order+1):
            lieterm = simplify(self.lieop(ham,vof0))
            for j in range(0,i-1):
                lieterm = simplify(self.lieop(ham,lieterm))

            #voft = voft + t**i / factorial(i) * lieterm # for my formalism, #zzz
            voft = voft + lieterm / factorial(i) # for Ems formalism, shouldn't each term also 

        return voft

    def funFromHam(self, ham, order, vof0):
        transresult = self.lietransform(ham, vof0, order)
        fun = simplify(transresult)
        return fun

    def hamToNumFuns(self, ham, kval, lval,order):
        xfun = self.funFromHam(ham, order, self.qx)
        xprimefun = self.funFromHam(ham, order, self.px)
        yfun = self.funFromHam(ham, order, self.qy)
        yprimefun = self.funFromHam(ham, order, self.py)
        zfun = self.funFromHam(ham, order, self.qz)
        zprimefun = self.funFromHam(ham, order, self.pz)

        xf = xfun.subs([(self.k,kval),(self.l,lval)])
        xpf = xprimefun.subs([(self.k,kval), (self.l,lval)])
        yf = yfun.subs([(self.k,kval),(self.l,lval)])
        ypf = yprimefun.subs([(self.k,kval), (self.l,lval)])
        zf = zfun.subs([(self.k,kval),(self.l,lval)])
        zpf = zprimefun.subs([(self.k,kval), (self.l,lval)])

        protoFunctions = (xf, xpf, yf, ypf, zf, zpf)
        symplectic = self.isSymplectic(protoFunctions)
        if not symplectic:
            print "Not symplectic!!!!!"

        xNumFun = lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),xf, "numpy")
        xpNumFun= lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),xpf, "numpy")
        yNumFun = lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),yf, "numpy")
        ypNumFun= lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),ypf, "numpy")
        zNumFun = lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),zf, "numpy")
        zpNumFun= lambdify((self.qx,self.px,self.qy,self.py,self.qz,self.pz),zpf, "numpy")

        numFuns = (xNumFun, xpNumFun, yNumFun, ypNumFun, zNumFun, zpNumFun)
        return numFuns

    def isSymplectic(self, protoFuns):
        size = len(protoFuns)/2 # numFun must be of even length
        S = np.append(np.append(np.zeros((size,size)),np.identity(size), axis=1),np.append(-np.identity(size),np.zeros((size,size)), axis=1), axis=0)

        J = np.array([
            [diff(protoFuns[0],self.qx), diff(protoFuns[0],self.qy), diff(protoFuns[0],self.qz), diff(protoFuns[0],self.px), diff(protoFuns[0],self.py), diff(protoFuns[0],self.pz)],
            [diff(protoFuns[2],self.qx), diff(protoFuns[2],self.qy), diff(protoFuns[2],self.qz), diff(protoFuns[2],self.px), diff(protoFuns[2],self.py), diff(protoFuns[2],self.pz)],
            [diff(protoFuns[4],self.qx), diff(protoFuns[4],self.qy), diff(protoFuns[4],self.qz), diff(protoFuns[4],self.px), diff(protoFuns[4],self.py), diff(protoFuns[4],self.pz)],
            [diff(protoFuns[1],self.qx), diff(protoFuns[1],self.qy), diff(protoFuns[1],self.qz), diff(protoFuns[1],self.px), diff(protoFuns[1],self.py), diff(protoFuns[1],self.pz)],
            [diff(protoFuns[3],self.qx), diff(protoFuns[3],self.qy), diff(protoFuns[3],self.qz), diff(protoFuns[3],self.px), diff(protoFuns[3],self.py), diff(protoFuns[3],self.pz)],
            [diff(protoFuns[5],self.qx), diff(protoFuns[5],self.qy), diff(protoFuns[5],self.qz), diff(protoFuns[5],self.px), diff(protoFuns[5],self.py), diff(protoFuns[5],self.pz)],
            ]) # from 1.80 in ref C

        approxOffsetQ = 0.01
        approxOffsetP = 0.00001
        for i in range(0,6):
            for j in range(0,6):
                J[i,j] = J[i,j].subs([(self.qx,approxOffsetQ),(self.qy,approxOffsetQ),(self.qz,approxOffsetQ),(self.px,approxOffsetP),(self.py,approxOffsetP),(self.pz,approxOffsetP)])

        J_transpose = np.transpose(J)
        LHS = np.dot(J_transpose, np.dot(S, J)) # from 1.81 in ref C
        LHSminusS = LHS - S
        detOfJ = np.linalg.det(J)

        ### Pick one of these
        ## symplecticity check with det of J
        if abs(detOfJ-1) <0.000001:
            return 1

        ## symplecticity check with 1.81 in ref C (checks if J^T*S*J = S)
        #if sum(sum(np.isclose(LHSminusS.astype(np.float64), np.zeros((2*size,2*size)).astype(np.float64)))) == LHSminusS.size: # if all elements are close the sum will be a sum of ones
        #    return 1

        return 0

# General class for elements from Hamiltonians, can be linear but since all is based on differential algebra "linear" is set to 0
class LieAlgElement(Element):
    def __init__(self, name, hamToUse, K, L, order, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        Element.__init__(self, "liealgelem " + name, 0)

        self.L = L
        self.K = K
        self.n = nbrOfSplits # matches matrix approach well
        #self.n = 5 # matches matrix approach not as well but is needed for split
        self.Lsp = L/self.n

        self.qx = Symbol('qx')
        self.qy = Symbol('qy')
        self.qz = Symbol('qz')
        self.px = Symbol('px')
        self.py = Symbol('py')
        self.pz = Symbol('pz')

        self.l = Symbol('l') # arbitrary length
        self.k = Symbol('k')

        self.LA = LieAlgebra() # Lie algebra object

        ## Hamiltonians
        self.driftham = -self.l/2*(self.px**2 + self.py**2 + self.pz**2)
        self.quadham = -self.l/2*(self.k**2*(self.qx**2-self.qy**2)+self.px**2+self.py**2+self.pz**2) # replace k with -k for defocus. Without quad term in z dir
        self.quaddefocusham = -self.l/2*(-self.k**2*(self.qx**2-self.qy**2)+self.px**2+self.py**2+self.pz**2) # replace k with -k for defocus. Without quad term in z dir
        self.sextupoleham = -self.l/2*(2/3*self.k**2*(self.qx**3-3*self.qx*self.qy**2)+(self.px**2+self.py**2)) # should the ps' perhaps be divided by 2 as in nonlinear2013_3.pdf? That division is assumed to be the l/2 in the beginning, . k is actually k**2
        self.octupoleham = -self.l/2*(2/4*self.k*(self.qx**4-6*self.qx**2*self.qy**2+self.qy**4)+(self.px**2+self.py**2)) # same decision as above

        if hamToUse == "drift":
            self.numFuns = self.LA.hamToNumFuns(self.driftham, K, self.Lsp, order)
        elif hamToUse == "quad":
            self.numFuns = self.LA.hamToNumFuns(self.quadham, K, self.Lsp, order)
        elif hamToUse == "quaddefocusham":
            self.numFuns = self.LA.hamToNumFuns(self.quaddefocusham, K, self.Lsp, order)
        elif hamToUse == "sextupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.sextupoleham, K, self.Lsp, order)
        elif hamToUse == "octupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.octupoleham, K, self.Lsp, order)

        self.hamUsed = hamToUse
        self.order = order

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('quad_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, multipart, twiss, beamdata)

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K) + "\t HamUsed: " +  str(self.hamUsed) + "\t Order: " +  str(self.order)

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            if self.spaceChargeOn:
                self.sc.updateMatrix(multipart,twiss)
                multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss

    def evaluateNumFun(self,multipart,envelope):
        for particle in multipart:
            #print "hej"
            #print "particle: " + str(particle)
            x = particle[0][0]
            xp = particle[0][1]
            y = particle[0][2]
            yp = particle[0][3]
            z = particle[0][4]
            zp = particle[0][5]
            #s = particle[1]

            particle[0][0] = self.numFuns[0](x, xp, y, yp, z, zp)
            particle[0][1] = self.numFuns[1](x, xp, y, yp, z, zp)
            particle[0][2] = self.numFuns[2](x, xp, y, yp, z, zp)
            particle[0][3] = self.numFuns[3](x, xp, y, yp, z, zp)
            particle[0][4] = self.numFuns[4](x, xp, y, yp, z, zp)
            particle[0][5] = self.numFuns[5](x, xp, y, yp, z, zp)

            particle[1] += self.Lsp
            # Doesn't change envelope at all
        envelope = envelopeFromMultipart(multipart)
        return multipart, envelope





## Leapfrog algorithm mostly from ref D.
def leapfrog(x_0, v_0, F, h, n):
    """ The leapfrog algorithm: x_0 is the single particle's x,y and z. v_0 is the single particle's xp,yp and zp. h is the step length. F() is a function (of just x,y,z or all (x,xp,y,yp,z, zp)?). n is the amount of steps """

    x_of_i = list()
    x_of_i.append(x_0) # first value in x_of_i
    
    v_of_i_plus_half = list()
    # the first value of v_of_i_plus_half is given in the first iteration of the velocity verlet

    v_of_i = list()
    v_of_i.append(v_0)

    # velocity verlet (better version of leapfrog). To get the wholes of v (v_of_i) from v_of_i_plus_half
    for i in range(0,n): # +1 to get a i to n as well
        
        v_of_i_plus_half.append(v_of_i[i] + 1/2*h*F(x_of_i[i])) # v_n_plus_half = v_n + 1/2*h*F(x_n)
        
        x_of_i.append(x_of_i[i] + h*v_of_i_plus_half[i]) # x_n_plus_one = x_n + h*v_n_plus_half
        
        v_of_i.append(v_of_i_plus_half[i] + 1/2*h*F(x_of_i[i+1])) # v_n_plus_one = v_n_plus_half + 1/2*h*F(x_n_plus_one)
        # send F(x_of_i[i+1]) to next lap to save some time

    return x_of_i, v_of_i






# Comes from ref E.
class Cavity(Element):
    def __init__(self, name, L, Ezofs, beamdata, nbrOfSplits):#, K, M):
        Element.__init__(self, "cavity " + name, 0) # linear set to zero
        #self.name = name
        #self.K = K
        #self.L = L
        #self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        #self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

        ## Input
        self.L = L
        self.n = nbrOfSplits
        # E-field
        self.oscillations = Ezofs[0]
        self.halfnbrofoscillations = self.oscillations/2
        self.amplitudeA = Ezofs[1]
        self.amplitudeB = Ezofs[2]
        self.E_0 = Ezofs[3]
        self.sigma = Ezofs[4]
        self.p = Ezofs[5]
        # beam
        self.beta_i = beamdata[0]
        self.rf_lambda = beamdata[1]
        self.m = beamdata[2]
        self.q = beamdata[3]
        E_i = beamdata[4]

        # calc params
        self.k = 2*constants.pi/self.beta_i/self.rf_lambda
        self.phi_s = self.findPhi_s(self.beta_i)
        print "self.phi_s: " + str(self.phi_s)

        ## Calculations
        T_of_beta = self.timeTransitFactor(self.beta_i)
        print "T_of_beta: " + str(T_of_beta)

        deltaW = self.q*self.E_0*T_of_beta*cos(self.phi_s)
        print "deltaW: " + str(deltaW)
        self.E_f = E_i + deltaW
        self.beta_f = betaFromE(self.m, self.E_f)
        beamdata[0] = self.beta_f

        gamma_i = 1/sqrt(1-self.beta_i**2)
        gamma_f = 1/sqrt(1-self.beta_f**2)

        beta_avg = (self.beta_i + self.beta_f)/2
        gamma_avg = (gamma_i + gamma_f)/2

        betaTimesGamma_i = self.beta_i*gamma_i
        betaTimesGamma_f = self.beta_f*gamma_f

        betaTimesGammaSquared_i = self.beta_i*gamma_i**2
        betaTimesGammaSquared_f = self.beta_f*gamma_f**2

        C = betaTimesGamma_i/betaTimesGamma_f

        k_11_x = self.E_0*T_of_beta*cos(self.phi_s)
        k_21_x = -self.q*self.k*self.E_0*T_of_beta*sin(self.phi_s)/(2*beta_avg*gamma_avg**2*self.m*constants.c**2)
        k_22_x = self.E_0*T_of_beta*cos(self.phi_s)

        G_x = np.array([[k_11_x*C, 0],
            [k_21_x/betaTimesGamma_f, k_22_x*C]])
        G_y = G_x # since in ref E. it says so after eqn 29

        k_21_z = self.q*self.k*self.E_0*T_of_beta*sin(self.phi_s)/(beta_avg**2*self.m*constants.c)
        G_z = np.array([[gamma_f/gamma_i, 0],
            [k_21_z/(gamma_i*betaTimesGammaSquared_f), betaTimesGammaSquared_i/betaTimesGammaSquared_f]])


        zeros_22 = np.zeros((2,2))
        row0 = np.hstack((np.hstack((G_x, zeros_22)),zeros_22))
        row1 = np.hstack((np.hstack((zeros_22, G_y)),zeros_22))
        row2 = np.hstack((np.hstack((zeros_22, zeros_22)),G_z))

        self.M = np.vstack((
            np.vstack(
            (row0, row1)
            ),
            row2))
        print str(self.M)

        self.Msp = self.M # splitting doesn't do anything
        self.Lsp = self.L

    # with E_z0 from ref C.
    def timeTransitFactor(self, beta):
        # E_z0_of_s
        z1z2 = -2*self.halfnbrofoscillations*self.L/3 # /3 since A = 0 and B =/= 0
        #print "z1z2: " + str(z1z2)
        z4z5 = 2*self.halfnbrofoscillations*self.L/3 # /3 since A = 0 and B =/= 0
        #print "z4z5: " + str(z4z5)
        ## Integral
        # -inf to z1|z2
        I1 = quad(lambda s: self.amplitudeB*exp(((s+z1z2)/self.sigma)**self.p)*cos(2*constants.pi/(beta*self.rf_lambda)*s - self.phi_s),-inf,z1z2)
        #print "I1: " + str(I1)
        # z2 to z4||z5
        I2 = quad(lambda s: (self.amplitudeA*cos(constants.pi*s/self.L)+self.amplitudeB*cos(3*constants.pi*s/self.L))*cos(2*constants.pi/(beta*self.rf_lambda)*s - self.phi_s),z1z2,z4z5)
        #print "I2: " + str(I2)
        # z5 to inf
        I3 = quad(lambda s: self.amplitudeB*exp((-(s-z4z5)/self.sigma)**self.p)*cos(2*constants.pi/(beta*self.rf_lambda)*s - self.phi_s),z4z5,inf)
        #print "I3: " + str(I3)

        # sum up
        res = I1[0]+I2[0]+I3[0]
        res = res/self.E_0
        return res

    def findPhi_s(self, beta):
        # E_z0_of_s
        z1z2 = -2*self.halfnbrofoscillations*self.L/3 # /3 since A = 0 and B =/= 0
        #print "z1z2: " + str(z1z2)
        z4z5 = 2*self.halfnbrofoscillations*self.L/3 # /3 since A = 0 and B =/= 0
        #print "z4z5: " + str(z4z5)
         
        Ithatgoeswithphi_s = inf
        for i in linspace(0,2*constants.pi,30):
            ## Integral
            # -inf to z1|z2
            I1 = quad(lambda s: self.amplitudeB*exp(((s+z1z2)/self.sigma)**self.p)*sin(2*constants.pi/(beta*self.rf_lambda)*s - i),-inf,z1z2)
            #print "I1: " + str(I1)
            # z2 to z4||z5
            I2 = quad(lambda s: (self.amplitudeA*cos(constants.pi*s/self.L)+self.amplitudeB*cos(3*constants.pi*s/self.L))*sin(2*constants.pi/(beta*self.rf_lambda)*s - i),z1z2,z4z5)
            #print "I2: " + str(I2)
            # z5 to inf
            I3 = quad(lambda s: self.amplitudeB*exp((-(s-z4z5)/self.sigma)**self.p)*sin(2*constants.pi/(beta*self.rf_lambda)*s - i),z4z5,inf)
            #print "I3: " + str(I3)
    
            # sum up
            res = I1[0]+I2[0]+I3[0]

            if abs(res) < Ithatgoeswithphi_s:
                phi_s = i
                Ithatgoeswithphi_s = abs(res)

        print "Ithatgoeswithphi_s: " + str(Ithatgoeswithphi_s)
        return phi_s

    def getNewE(self):
        return self.E_f

    def printInfo(self):
        return self.name + "\t L: " + str(self.L) + "\t Oscillations: " + str(self.oscillations) + "\t AmplitudeA: " + str(self.amplitudeA) + "\t AmplitudeB: " + str(self.amplitudeB) + "\t E_0: " + str(self.E_0) + "\t sigma: " + str(self.sigma) + "\t p: " + str(self.p)

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            #if self.spaceChargeOn:
                #self.sc.updateMatrix(multipart,twiss)
                #multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        return multipart, envelope


class FieldMapCavity(Element):
    def __init__(self, name, L, beamdata, nbrOfSplits):
        Element.__init__(self, "fieldmap " + name, 1)
        
        ## Input
        self.L = L
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n

        self.beta_s = beamdata[0]
        self.gamma_s = 1/sqrt(1-self.beta_s**2)
        self.m = beamdata[2]
        self.q = beamdata[3]

        #notdone#self.E_z =

        ## Calculation
        #dfdx = f.diff(x)

        #notdone#self.dE_xoverdx = 
        #notdone#self.dE_yoverdy = 
        #notdone#self.dE_zoverdz = self.E_z.diff(z)

        #notdone#self.dB_yoverdx =
        #notdone#self.dB_xoverdy = 

        self.prefactor = self.q/(self.gamma_s*self.beta_s**2*self.m*constants.c**2)

        self.Msp = np.array([
            [1,0,0,0,0,0],
            [self.Lsp*self.prefactor*(self.dE_xoverdx-self.beta_s*constants.c*self.dB_yoverdx),1-self.Lsp*self.prefactor*self.E_z,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,self.Lsp*self.prefactor*(self.dE_yoverdy+self.beta_s*constants.c*self.dB_xoverdy),1-self.Lsp*self.prefactor*self.E_z,0,0],
            [0,0,0,0,1,0],
            [0,0,0,0,self.Lsp*self.prefactor*1/self.gamma_s**2*self.dE_zoverdz,1-self.Lsp*self.prefactor*self.E_z]
            ])

# references
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition
# A. 7.2. Space Charge Impulses in simulatingbeamswithellipsoidalsymmetry-secondedition
# B. A MODIFIED QUADSCAN TECHNIQUE FOR EMITTANCE.pdf
# C. Accelerator-Recipies.pdf by E. Laface
# D. The leapfrog method and other symplectic algorithms for integrating Newtons laws of motion Peter Young Dated April 21 2014
# E. ESS Linac simulator
# F. WEPEA 040