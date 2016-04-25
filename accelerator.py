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
from relativity import betaFromE, gammaFromBeta, EFromBeta
import copy
import time

# beamdata comes as [beta, rf_lambda, m, q, E, nbrOfParticles, I]
# the units for beamdata: beta is unitless, rf_lambda is in m, m is in kg, q is in C, E is in J (should be in Mev later), nbrOfParticles is unitless, I is Ameperes
# twiss comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
# the units for twiss: alpha is unitless, beta is in m, epsilon is in m*rad
# multipart is an array of particles which comes as [[x, xp, y, yp, z, zp], s] with unit m for x, y, z and s. xp, yp and zp are unitless.
class Lattice:
    def __init__(self,name,beamdata,twiss,multipart):
        self.name = name
        self.lattice = list()

        ## Beam properties
        self.beamdata = beamdata
        self.twiss = twiss
        self.multipart = multipart

        ## space-charge, not sure how these shall be represented in the GUI. spaceChargeOn can be a radio button, nbrOfSplits can be a text field
        self.spaceChargeOn = 0
        self.nbrOfSplits = 1

        ## laps to be made
        self.laps = 1

    def createDrift(self, name, L):
        drift = Drift(name, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(drift)

    def createDipole(self, name, rho, alpha, n):
        dipole = Dipole(name, rho, alpha, n, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(dipole)

    def createQuadrupole(self, name, K, L):
        quad = Quad(name, K, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(quad)

    def createLieDrift(self, name, L, compOrder):
        hamToUse = "driftham"
        liedrift = LieAlgElement(name, hamToUse, 0, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(liedrift)

    def createSextupole(self, name, K, L, compOrder):
        hamToUse = "sextupoleham"
        sextu = LieAlgElement(name, hamToUse, K, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(sextu)

    def createOctupole(self, name, K, L, compOrder):
        hamToUse = "octupoleham"
        octu = LieAlgElement(name, hamToUse, K, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(octu)

    def createSextupolerel(self, name, K, L, compOrder):
        hamToUse = "sextupolehamrel"
        sextu = LieAlgElement(name, hamToUse, K, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(sextu)

    def createSextupolekin(self, name, K, L, compOrder):
        hamToUse = "sextupolehamkin"
        sextu = LieAlgElement(name, hamToUse, K, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(sextu)

    def createOctupolekin(self, name, K, L, compOrder):
        hamToUse = "octupolehamkin"
        octu = LieAlgElement(name, hamToUse, K, L, compOrder, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(octu)

    def createSextupolemat(self, name, K, L):
        sextu = SextupoleMat(name, K, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(sextu)

    def createSextupolematema(self, name, K, L):
        sextu = SextupoleMatEma(name, K, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(sextu)

    def createRotation(self, name, nu_x, nu_y):
        rot = Rotation(name, nu_x, nu_y)
        self.appendElement(rot)

    def createCavity(self, name, L, Ezofs):
        cavity = Cavity(name, L, Ezofs, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        self.appendElement(cavity)
        self.beamdata[4] = cavity.getNewE() # Update energy

    def createCavityMatrix(self, name, L, a, Efield_0, phi_0):
        cavity = CavityMatrix(name, L, a, Efield_0, phi_0, self.beamdata)
        self.appendElement(cavity)

    def appendElement(self, element):
        self.lattice.append(element)
        return 1

    ## Passing of varibles to and fro
    def getLattice(self):
        return self.lattice

    def getBeamdata(self):
        return self.beamdata

    def getTwiss(self):
        return self.twiss

    def getMultipart(self):
        return self.multipart

    def getSpaceChargeOn(self):
        return self.spaceChargeOn

    def getNbrOfSplits(self):
        return self.nbrOfSplits

    def setBeamdata(self,beamdata):
        self.beamdata = beamdata

    def setTwiss(self,twiss):
        self.twiss = twiss

    def setMultipart(self,multipart):
        self.multipart = multipart

    def setSpaceChargeOnAndSplits(self, spaceChargeOn, nbrOfSplits):
        self.spaceChargeOn = spaceChargeOn
        self.nbrOfSplits = nbrOfSplits
        # Now the lattice shall remake each element with the new spaceChargeOn and nbrOfSplits
        for elem in self.lattice:
            elem.updateSC(self.spaceChargeOn, self.nbrOfSplits, self.multipart, self.twiss, self.beamdata)

    def getLaps(self):
        return self.laps

    def setLaps(self, laps):
        self.laps = laps
    ## End Passing

    def printLattice(self):
        text = ""
        for elem in self.lattice:
            text = text + elem.printInfo() + "\n"
        return text

    def printMatrices(self):
        text = ""
        for elem in self.lattice:
            text = text + elem.name + ":\n"
            hasFun = getattr(elem, "printMatrix", None)
            if not callable(hasFun):
                continue
            text = text + elem.printMatrix() + "\n"
        return text

    def evaluate(self, multipartin,envelopein,twissin):
        print "Evaluation started..."
        t0 = time.clock()
        multipart = copy.deepcopy(multipartin)
        envelope = copy.deepcopy(envelopein)
        twiss = copy.deepcopy(twissin)
        envlist = list()
        envlist.append(np.array([envelope, 0]))
        multipartafterall = copy.deepcopy(multipartin)
        if self.spaceChargeOn:
            for elem in self.lattice:
                print elem.printInfo()
                multipart,envelope, twiss, env_with_s = elem.evaluateWithSC(multipart,envelope,twiss)
                env_with_s[1] = env_with_s[1] + envlist[-1][1]
                envlist.append(env_with_s)
        else:
            for i in range(0,self.laps):
                for elem in self.lattice:
                    #print elem.printInfo()
                    if elem.printInfo().startswith("cavmat"):
                        multipart,envelope, twiss, env_with_s, self.beamdata = elem.evaluateWithoutSC(multipart,envelope,twiss, self.beamdata)
                    else:
                        multipart,envelope, twiss, env_with_s = elem.evaluateWithoutSC(multipart,envelope,twiss)
                    env_with_s[1] = env_with_s[1] + envlist[-1][1]
                    envlist.append(env_with_s)
                    #for i in range(len(multipart)): # for saving the particles after each element
                    #    multipartafterall.append(copy.deepcopy(multipart[i]))
        t1 = time.clock()
        tdiff = t1-t0
        print "Evaluation finished, time = " + str(tdiff)
        return multipart,envelope, twiss, envlist, multipartafterall

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
        gamma = gammaFromBeta(beamdata[0])
        self.M = self.createMatrixM(L, beamdata[0], gamma) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix

        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        
        self.Msp = self.createMatrixM(self.Lsp, beamdata[0], gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('drift_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('drift_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, beamdata)

    def printInfo(self):
        return self.name + "\t L: " + str(self.L)

    def printMatrix(self):
        MspAsMatrix = np.asmatrix(self.Msp)
        return str(MspAsMatrix**self.n)

    def createMatrixM(self,L, beta, gamma):
        return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,L/(beta**2*gamma**2)],
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

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        
        gamma = gammaFromBeta(beamdata[0])
        self.Msp = self.createMatrixM(self.Lsp, beamdata[0], gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('drift_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('drift_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('drift_sc', self.Lsp, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):  
            multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
            multipart, envelope, env_with_s = self.evaluateMT(multipart,envelope) # use the new data for "normal" evaluation
            
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
            multipart, envelope, env_with_s = self.evaluateMT(multipart,envelope) # use the new data for "normal" evaluation
            
        return multipart, envelope, twiss, env_with_s

    def evaluateMT(self,multipart,envelope):
        # should just go through a disunited part
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j][0][0:6] = np.dot(self.Msp, multipart[j][0][0:6])
            multipart[j][1] = multipart[j][1] + self.Lsp
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

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

        self.Msp = self.createMatrixM(self.rho, self.Lsp, self.K_x, self.K_y, self.beta, self.gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('dipole_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('dipole_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('dipole_sc', self.Lsp, beamdata)

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

    def printInfo(self):
        return self.name + "\t rho: " + str(self.rho) + "\t L: " + str(self.L) + "\t Alpha: " + str(self.alpha) + "\t nparam: " + str(self.nparam)
        #return self.name + "\t rho: " + str(self.rho) + "\t L: " + str(self.L) + "\t K_x: " + str(self.K_x) + "\t K_y: " + str(self.K_y) + "\t beta: " + str(self.beta) + "\t Alpha: " + str(self.alpha) + "\t nparam: " + str(self.nparam)

    def printMatrix(self):
        MspAsMatrix = np.asmatrix(self.Msp)
        return str(MspAsMatrix**self.n)        

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        
        self.Msp = self.createMatrixM(self.rho, self.Lsp, self.K_x, self.K_y, self.beta, self.gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('dipole_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('dipole_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('dipole_sc', self.Lsp, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j][0][0:6] = np.dot(self.Msp, multipart[j][0][0:6])
            multipart[j][1] = multipart[j][1] + self.Lsp
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

### QUAD
class Quad(LinearElement):
    def __init__(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        LinearElement.__init__(self, "quad " + name)
        #self.name = name
        self.K = K
        self.L = L
        gamma = gammaFromBeta(beamdata[0])
        self.M = self.createMatrixM(K, L, beamdata[0], gamma) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix
        
        # disunite matrices
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        self.Msp = self.createMatrixM(self.K, self.Lsp, beamdata[0], gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('quad_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('quad_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('quad_sc', self.Lsp, beamdata)

    """ taken from lie code where K = sqrt(k)"""
    def createMatrixM(self,K,L,beta,gamma):
        if K > 0: # defocus in x
            return np.array([
                [cosh(K*L),sinh(K*L)/K,0,0,0,0],
                [K*sinh(K*L),cosh(K*L),0,0,0,0],
                [0,0,cos(K*L),sin(K*L)/K,0,0],
                [0,0,-K*sin(K*L),cos(K*L),0,0],
                [0,0,0,0,1,L/(beta**2*gamma**2)],
                [0,0,0,0,0,1]
                ])
        elif K < 0: # focus in x
            return np.array([
                [cos(K*L),sin(K*L)/K,0,0,0,0],
                [-K*sin(K*L),cos(K*L),0,0,0,0],
                [0,0,cosh(K*L),sinh(K*L)/K,0,0],
                [0,0,K*sinh(K*L),cosh(K*L),0,0],
                [0,0,0,0,1,L/(beta**2*gamma**2)],
                [0,0,0,0,0,1]
                ])
        else:
            return np.array([
                [1,L,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,L,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,L/(beta**2*gamma**2)],
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

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K)

    def printMatrix(self):
        MspAsMatrix = np.asmatrix(self.Msp)
        return str(MspAsMatrix**self.n)

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        
        gamma = gammaFromBeta(beamdata[0])
        self.Msp = self.createMatrixM(self.K, self.Lsp, beamdata[0], gamma)
        self.Tsp = self.createMatrixT(self.Msp)

        # space charge class
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            #self.sc = SpaceCharge('quad_sc', self.Lsp, multipart, twiss, beamdata) # OLD
            self.sc = SpaceCharge('quad_sc', self.Lsp, beamdata) # NEW
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('quad_sc', self.Lsp, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        #print "hej fran quad"
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j][0][0:6] = np.dot(self.Msp, multipart[j][0][0:6])
            multipart[j][1] = multipart[j][1] + self.Lsp
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

### SPACE CHARGE!!!!! C. Allen's approach
class SpaceCharge(LinearElement):
    #def __init__(self, name, deltas, multipart, twiss, beamdata): # old style
    def __init__(self, name, deltas, beamdata): # new style
        LinearElement.__init__(self, name)
        self.deltas = deltas
        self.beamdata = beamdata

        ## NEW STUFF
        self.beta = beamdata[0]
        self.m = beamdata[2]
        self.q = beamdata[3]
        self.N = beamdata[5]
        self.I = beamdata[6]
        self.gamma = gammaFromBeta(self.beta)
        self.bc2 = self.gamma*self.gamma-1 # wierd gamma**2 that CKA uses

        ## END NEW STUFF

        #########self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)

    #def updateBeam(self, twiss):
    #    self.twiss = twiss
    #    return

    #########def beamChanged(self, newtwiss):
    #########    threshold = 0.1
    #########    if diffBtwBeams(self.twiss, newtwiss) > threshold:
    #########        return 1
    #########    else:
    #########        return 0
#########
    #########def diffBtwBeams(self, twiss1,twiss2):
    #########    diff = 0
#   #########     for bla
#   #########         diff += diff_each_variable
    #########    return diff
#########
    #########def updateMatrix(self,multipart,twiss):
    #########    self.twiss = twiss
    #########    self.Msc = self.spaceChargeMatrix(multipart,twiss, self.beamdata)
    #########    return 1

    def R_D(self, x, y, z):
        # from (110) in ref 1.
        result = quad(lambda t : 3/2*1/(sqrt(t+x) * sqrt(t+y) * (t+z)**(3/2)), 0, inf)
        # result[0] is the result and result[1] is the error
        return result[0]

    #def spaceChargeMatrix(self, multipart, twiss, beamdata, envelope): # OLD
    def spaceChargeMatrix(self, multipart, envelope): # NEW
        ########## beamdata: beta (speed), mass, charge, lambda (RF-wavelength)
        #########beta = beamdata[0]
        #########rf_lambda = beamdata[1]
        #########m = beamdata[2]
        #########q = beamdata[3]
#########
        #########I = beamdata[6]
        #########
        ########## beam data is needed as input to calculate the following variables...
        #########N = len(multipart) # this info should come from the multipart (len(multipart))
        #########Q = q*N # from (19) in ref 1.
        #########gamma = 1/sqrt(1-beta**2)
        #########c = constants.c # in metric (metric for everything perhaps?)
        #########vac_perm = constants.epsilon_0
#########
        ##########I = N*q*c/rf_lambda # from ref. E
#########
        ########### Courant-Snyder or Twiss params
        ########## envelope comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
        #########alpha_x = twiss[0]
        #########beta_x = twiss[1]
        #########epsilon_rms_x = twiss[2]
        #########alpha_y = twiss[3]
        #########beta_y = twiss[4]
        #########epsilon_rms_y = twiss[5]
        #########alpha_z = twiss[6]
        #########beta_z = twiss[7]
        #########epsilon_rms_z = twiss[8]
#########
        ########### envelope X, Xp, Y, Yp, Z and Zp
        #########X = sqrt(5*beta_x*epsilon_rms_x)
        #########Xp = -alpha_x*sqrt(5*epsilon_rms_x/beta_x)
#########
        #########Y = sqrt(5*beta_y*epsilon_rms_y)
        #########Yp = -alpha_y*sqrt(5*epsilon_rms_y/beta_y)
#########
        #########Z = sqrt(5*beta_z*epsilon_rms_z)
        ##########print "beta_z: " + str(beta_z)
        #########Zp = -alpha_z*sqrt(5*epsilon_rms_z/beta_z)
#########
        ########## <.> is called "norm_of_."
        ########## <x^2> (norm_of_xx), <xx'> (norm_of_xx') and <x'^2> (norm_of_xpxp) come from (119) in ref 1.
        #########norm_of_xx = 1/5*X**2
        ########## <x> (norm_of_x) come from (81) in ref 1.
        #########norm_of_x = 0.0
        ########## <xx'> (norm_of_xxp) come from (119) in ref 1. NOT used here but here as a reminder
        #########norm_of_xxp = 1/5*X*Xp
        ########## <x'x'> (norm_of_xpxp) come from (119) in ref 1. NOT used here but here as a reminder
        #########norm_of_xpxp = 1/5*Xp**2+5*epsilon_rms_x**2/X**2
#########
        #########norm_of_yy = 1/5*Y**2
        #########norm_of_y = 0.0
        #########norm_of_yp = 1/5*Y*Yp
        #########norm_of_ypyp = 1/5*Yp**2+5*epsilon_rms_y**2/Y**2
#########
        #########norm_of_zz = 1/5*Z**2
        #########norm_of_z = 0.0
        #########norm_of_zp = 1/5*Z*Zp
        #########norm_of_zpzp = 1/5*Zp**2+5*epsilon_rms_z**2/Z**2
#########
        ########## <xE_x> (norm_of_xE_x) come from (115) in ref 1.
        #########norm_of_xE_x = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_xx*self.R_D(norm_of_yy, norm_of_zz, norm_of_xx)
        #########norm_of_yE_y = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_yy*self.R_D(norm_of_xx, norm_of_zz, norm_of_yy)
        #########norm_of_zE_z = 1/5**(3/2)*Q/(4*math.pi*vac_perm)*norm_of_zz*self.R_D(norm_of_xx, norm_of_yy, norm_of_zz)
        ########## R_D_bracket() = R_D_() from (111) and (112) in ref 1.
#########
        ########## eqn (152) in ref 1. But maybe it should be (157) in ref 1. instead??????????
        #########f_scx = gamma**3*beta**2*m*c**2/q*(norm_of_xx - norm_of_x**2)/norm_of_xE_x # eqn (152) in ref 1. But maybe it should be (157) in ref 1. instead??????????
        #########f_scy = gamma**3*beta**2*m*c**2/q*(norm_of_yy - norm_of_y**2)/norm_of_yE_y # (152) (with x->y) in ref 1.
        #########f_scz = gamma**3*beta**2*m*c**2/q*(norm_of_zz - norm_of_z**2)/norm_of_zE_z # (152) (with x->z) in ref 1.        
#########
        ########## Mean of x,y and z from all the particles
        #########xbar = sum([multipart[i][0][0] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting x out of each particle and then taking the mean
        #########ybar = sum([multipart[i][0][2] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting y out of each particle and then taking the mean
        #########zbar = sum([multipart[i][0][4] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting z out of each particle and then taking the mean
#########
        ########## Matrix eqn (154) in ref 1.
        #########Msc = np.array([
        #########        [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
        #########        [self.deltas/f_scx,1.0,0.0,0.0,0.0,0.0,-xbar*self.deltas/f_scx],
        #########        [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
        #########        [0.0,0.0,self.deltas/f_scy,1.0,0.0,0.0,-ybar*self.deltas/f_scy],
        #########        [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
        #########        [0.0,0.0,0.0,0.0,self.deltas/f_scz,1.0,-zbar*self.deltas/f_scz],
        #########        [0.0,0.0,0.0,0.0,0.0,0.0,1.0]
        #########    ])

        ## NEW STUFF
        #K = self.q**2*self.N/2/constants.pi/constants.epsilon_0/self.gamma/self.bc2/self.beta**2/self.m/constants.c**2 # (18) (bunched beam) with mods q*q -> q**2 and uses CKA's bc2 instead of a gamma**2
        #K = self.I*self.q/2/constants.pi/constants.epsilon_0/gamma**3/beta**3/self.m/constants.c**3 # (20) (continous beam)
        ## NEWNEW STUFF
        K = self.I*self.q/2/constants.pi/constants.epsilon_0/self.gamma/self.bc2/self.beta**2/self.m/constants.c**2 # I = q*N

        sigma_x = sqrt(envelope[0])
        sigma_y = sqrt(envelope[3])
        sigma_z = sqrt(envelope[6])

        # R_D[X**2/Z**2,Y**2/Z**2,1] = Z**3*R_D(X**2,Y**2,Z**2)
        one_over_f_scx = K/2*(5*sigma_x)**(-3.0/2)*sigma_x**3*self.R_D(sigma_y**2,sigma_z**2,sigma_x**2) # negative in exponent instead of 1/... . *sigma_x**3 since (112)
        one_over_f_scy = K/2*(5*sigma_y)**(-3.0/2)*sigma_y**3*self.R_D(sigma_z**2,sigma_x**2,sigma_y**2)
        one_over_f_scz = K/2*(5*sigma_z)**(-3.0/2)*sigma_z**3*self.R_D(sigma_x**2,sigma_y**2,sigma_z**2)

        # Mean of x,y and z from all the particles
        xbar = sum([multipart[i][0][0] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting x out of each particle and then taking the mean
        ybar = sum([multipart[i][0][2] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting y out of each particle and then taking the mean
        zbar = sum([multipart[i][0][4] for i in xrange(len(multipart))])/(len(multipart)) # Tedious way of getting z out of each particle and then taking the mean

        Msc = np.array([
                [1.0,0.0,0.0,0.0,0.0,0.0,0.0],
                [self.deltas*one_over_f_scx,1.0,0.0,0.0,0.0,0.0,-xbar*self.deltas*one_over_f_scx],
                [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,self.deltas*one_over_f_scy,1.0,0.0,0.0,-ybar*self.deltas*one_over_f_scy],
                [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,self.deltas*one_over_f_scz,1.0,-zbar*self.deltas*one_over_f_scz],
                [0.0,0.0,0.0,0.0,0.0,0.0,1.0]
            ])

        #R = ... # eqn (161)
        #Msc = np.dot(R.transpose(),np.dot(Msc, R))

        return Msc

    def evaluateSC(self,multipart,envelope):
        #self.updateMatrix(multipart,twiss)
        self.Msc = self.spaceChargeMatrix(multipart, envelope)
        for j in range(0,len(np.atleast_1d(multipart))):
            # The should be a check if Msc and Tsc need to be updated if the beam properties have changed a lot!!!!!!!!
            #if beamChanged(envelope):
                #self.Msc, self.Tsc = spaceChargeMatrix(envlope)

            extendedphasespace = np.append(multipart[j][0][0:6], 1) # adds the dispersion 1 term
            extendedphasespace = np.dot(self.Msc, extendedphasespace) # here calculations are made
            multipart[j][0][0:6] = extendedphasespace[0:6] # throws away the dispersion 1 term # s remains the same because the particles don't go anywhere. They "go" in evaluateM()
        envelope = envelopeFromMultipart(multipart) # the envelope is just calculated from the particles (NOT ON ITS OWN)
        #env_with_s = np.array([copy.deepcopy(envelope), self.Lsp]) # there is only an angle kick so the deviation envelope wont change
        return multipart, envelope






### SPACE CHARGE!!!!! Elliptical integral approach
class SpaceChargeEllipticalIntegral(LinearElement):
    def __init__(self, name, deltas, beamdata):
        LinearElement.__init__(self, name)
        self.deltas = deltas
        self.beamdata = beamdata

        # beamdata: beta (speed), lambda (RF-wavelength), mass, charge
        self.beta = beamdata[0]
        self.rf_lambda = beamdata[1]
        self.m = beamdata[2]
        self.q = beamdata[3]
        self.N = beamdata[5]
        self.I = beamdata[6]
        
        # beam data is needed as input to calculate the following variables...
        self.Q = self.q*self.N # from (19) in ref 1.
        self.gamma = 1/sqrt(1-self.beta**2)

#    def beamChanged(self, newtwiss):
#        threshold = 0.1
#        if diffBtwBeams(self.twiss, newtwiss) > threshold:
#            return 1
#        else:
#            return 0

#    def diffBtwBeams(self, twiss1,twiss2):
#        diff = 0
#        for bla
#            diff += diff_each_variable
#        return diff

    def spaceChargeMatrix(self, envelope):        
        ## The semi-axes are the to the sigmas
        r_x = sqrt(envelope[0])
        r_y = sqrt(envelope[3])
        r_z = sqrt(envelope[6])

        # Eliptical integral
        g = self.gamma*r_z/sqrt(r_x*r_y) # eqn 40 from ref E.
        f_of_g_integral = g/2*quad(lambda t : 1/((t+1)*(t+g**2)**(3/2)), 0, inf)[0] # eqn 41 from ref E.

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

        U_scx_without_x = self.I*self.rf_lambda*G_x_without_x/(4*math.pi*constants.epsilon_0*constants.c*self.gamma**2)
        U_scy_without_y = self.I*self.rf_lambda*G_y_without_y/(4*math.pi*constants.epsilon_0*constants.c*self.gamma**2)
        U_scz_without_z = self.I*self.rf_lambda*G_z_without_z/(4*math.pi*constants.epsilon_0*constants.c)

        #delta_P_x = q*U_scx*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_x: " + str(delta_P_x)
        #delta_P_y = q*U_scy*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_y: " + str(delta_P_y)
        #delta_P_z = q*U_scz*self.deltas/(m*c**2*beta) # eqn 42 from ref E.
        #print "delta_P_z: " + str(delta_P_z)

        delta_P_x_without_x = self.q*U_scx_without_x*self.deltas/(self.m*constants.c**2*self.beta)
        delta_P_y_without_y = self.q*U_scy_without_y*self.deltas/(self.m*constants.c**2*self.beta)
        delta_P_z_without_z = self.q*U_scz_without_z*self.deltas/(self.m*constants.c**2*self.beta)

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

        Msc = np.array([
                [1.0,0.0,0.0,0.0,0.0,0.0],
                [delta_xp_without_x,1.0,0.0,0.0,0.0,0.0],
                [0.0,0.0,1.0,0.0,0.0,0.0],
                [0.0,0.0,delta_yp_without_y,1.0,0.0,0.0],
                [0.0,0.0,0.0,0.0,1.0,0.0],
                [0.0,0.0,0.0,0.0,delta_zp_without_z,1.0]
            ])

        return Msc

    def evaluateSC(self,multipart,envelope):
        self.Msc = self.spaceChargeMatrix(envelope) # Always update Msc
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j][0][0:6] = np.dot(self.Msc, multipart[j][0][0:6])
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
            #voft = voft + lieterm / factorial(i) # for Ems formalism, (For all the old hamiltonians)
            voft = voft + (-self.l)**i * lieterm / factorial(i) # adds L**i (used for relativistic hamiltonians (those with sqrt:s))

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
        print "distance from symplecticity: " + str(abs(detOfJ-1))
        return 0

# General class for elements from Hamiltonians, can be linear but since all is based on differential algebra "linear" is set to 0
class LieAlgElement(Element):
    def __init__(self, name, hamToUse, K, L, order, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        Element.__init__(self, "liealgelem " + name, 0)

        self.hamUsed = hamToUse
        self.hamToUse = hamToUse
        self.L = L
        self.K = K
        self.order = order
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

        ## Hamiltonians, works with Lie transform which says "Ems formalism"
        #self.driftham = -self.l/2*(self.px**2 + self.py**2 + self.pz**2)
        self.quadham = -self.l/2*(self.k**2*(self.qx**2-self.qy**2)+self.px**2+self.py**2+self.pz**2) # replace k with -k for defocus. Without quad term in z dir
        self.quaddefocusham = -self.l/2*(-self.k**2*(self.qx**2-self.qy**2)+self.px**2+self.py**2+self.pz**2) # replace k with -k for defocus. Without quad term in z dir
        #self.sextupoleham = -self.l/2*(2/3*self.k**2*(self.qx**3-3*self.qx*self.qy**2)+(self.px**2+self.py**2)) # should the ps' perhaps be divided by 2 as in nonlinear2013_3.pdf? That division is assumed to be the l/2 in the beginning, . k is actually k**2
        #self.octupoleham = -self.l/2*(2/4*self.k*(self.qx**4-6*self.qx**2*self.qy**2+self.qy**4)+(self.px**2+self.py**2)) # same decision as above

        ## redefined Hamiltonians (works!!!!! for sextu, ) works with the "relativistic" transform
        self.driftham = 1/2*(self.px**2 + self.py**2 + self.pz**2)
        self.sextupoleham = 1/6*self.k*(self.qx**3-3*self.qx*self.qy**2)+1/2*(self.px**2+self.py**2) # works with lie trans for the rel (which is the better trans)
        self.octupoleham = 1/8*self.k*(self.qx**4-6*self.qx**2*self.qy**2+self.qy**4)+1/2*(self.px**2+self.py**2)

        ## Hamiltonians with only the kinetic parts (used for symplectic integrator described in Wolski ch 10)
        self.sextupolehamkin = 1/6*self.k*(self.qx**3-3*self.qx*self.qy**2)
        self.octupolehamkin = 1/8*self.k*(self.qx**4-6*self.qx**2*self.qy**2+self.qy**4)

        ## Relativistic Hamiltonian from Wolski
        #beta = beamdata[0]
        #gamma = gammaFromBeta(beta)
        #beta = 1 # ultra relativistic
        #gamma = 1e20 # ultra relativistic
        # delta =approx= gamma**2*z' from tracewins conversion section
        #self.sextupolehamrel = gamma**2*self.pz/beta - sqrt((1/beta + gamma**2*self.pz)**2 - self.px**2 - self.py**2 - 1/(beta**2*gamma**2)) + self.k/6*(self.qx**3 - 3*self.qx*self.qy**2) # (9.47)
        # delta =approx= z'
        #self.sextupolehamrel = self.pz/beta - sqrt((1/beta + self.pz)**2 - self.px**2 - self.py**2 - 1/(beta**2*gamma**2)) + self.k/6*(self.qx**3 - 3*self.qx*self.qy**2) # (9.47)
        # only x and xp (y, yp, z and delta = 0), with beta = 1 and gamma = inf
        self.sextupolehamrel = -sqrt(1 - self.px**2) + self.k/6*(self.qx**3) # (9.47) # only in x dimension!!!!
        #self.octupolehamrel = -sqrt(1 - self.px**2) + 2/8*self.k*(self.qx**4)

        if self.hamToUse == "driftham":
            self.numFuns = self.LA.hamToNumFuns(self.driftham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "quad":
            self.numFuns = self.LA.hamToNumFuns(self.quadham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "quaddefocusham":
            self.numFuns = self.LA.hamToNumFuns(self.quaddefocusham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.sextupoleham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "octupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.octupoleham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupolehamrel":
            self.numFuns = self.LA.hamToNumFuns(self.sextupolehamrel, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupolehamkin":
            self.numFuns = self.LA.hamToNumFuns(self.sextupolehamkin, self.K, self.Lsp, self.order)
        elif self.hamToUse == "octupolehamkin":
            self.numFuns = self.LA.hamToNumFuns(self.octupolehamkin, self.K, self.Lsp, self.order)

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K) + "\t HamUsed: " +  str(self.hamUsed) + "\t Order: " +  str(self.order)

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n
        
        if self.hamToUse == "driftham":
            self.numFuns = self.LA.hamToNumFuns(self.driftham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "quad":
            self.numFuns = self.LA.hamToNumFuns(self.quadham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "quaddefocusham":
            self.numFuns = self.LA.hamToNumFuns(self.quaddefocusham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.sextupoleham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "octupoleham":
            self.numFuns = self.LA.hamToNumFuns(self.octupoleham, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupolehamrel":
            self.numFuns = self.LA.hamToNumFuns(self.sextupolehamrel, self.K, self.Lsp, self.order)
        elif self.hamToUse == "sextupolehamkin":
            self.numFuns = self.LA.hamToNumFuns(self.sextupolehamkin, self.K, self.Lsp, self.order)
        elif self.hamToUse == "octupolehamkin":
            self.numFuns = self.LA.hamToNumFuns(self.octupolehamkin, self.K, self.Lsp, self.order)

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            self.sc.updateMatrix(multipart,twiss)
            multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateNumFun(self,multipart,envelope):
        for particle in multipart:
            x = particle[0][0]
            xp = particle[0][1]
            y = particle[0][2]
            yp = particle[0][3]
            z = particle[0][4]
            zp = particle[0][5]

            particle[0][0] = self.numFuns[0](x, xp, y, yp, z, zp)
            particle[0][1] = self.numFuns[1](x, xp, y, yp, z, zp)
            particle[0][2] = self.numFuns[2](x, xp, y, yp, z, zp)
            particle[0][3] = self.numFuns[3](x, xp, y, yp, z, zp)
            particle[0][4] = self.numFuns[4](x, xp, y, yp, z, zp)
            particle[0][5] = self.numFuns[5](x, xp, y, yp, z, zp)

            particle[1] += self.Lsp
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

class SextupoleMat(Element): # Based on wolski ch 10, mainly eqn (10.6)
    def __init__(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        Element.__init__(self, "liealgelem mat " + name, 0)

        self.L = L
        self.K = K
        self.n = nbrOfSplits # matches matrix approach well
        self.Lsp = L/self.n

        self.qx = Symbol('qx')
        self.px = Symbol('px')

        xf = self.qx + self.L*self.px/sqrt(1-self.px**2)
        xpf = self.px - self.K*self.L*self.qx**2/2 - self.K*self.L**2*self.qx*self.px/sqrt(1-self.px**2) - self.K*self.L**3*self.px**2/(2*(1-self.px**2))
        #xfun = -self.k*self.l**2*self.qx**2/(4*(-self.px**2 + 1)**(3/2)) + self.l*self.px/sqrt(-self.px**2 + 1) + self.qx # moved here from a temp place in LieAlg
        #xprimefun = -self.k*self.l**2*self.px*self.qx/(2*sqrt(-self.px**2 + 1)) + self.k*self.l*self.qx**2/2 + self.px # moved here from a temp place in LieAlg

        self.xNumFun = lambdify((self.qx,self.px),xf, "numpy")
        self.xpNumFun = lambdify((self.qx,self.px),xpf, "numpy")

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K)

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            self.sc.updateMatrix(multipart,twiss)
            multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluate(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateNumFun(self,multipart,envelope):
        for particle in multipart:

            x = particle[0][0]
            xp = particle[0][1]

            particle[0][0] = self.xNumFun(x, xp)
            particle[0][1] = self.xpNumFun(x, xp)

            particle[1] += self.Lsp
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

class SextupoleMatEma(Element): # from ema's lecture 4
    def __init__(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        Element.__init__(self, "liealgelem matema " + name, 0)

        self.L = L
        self.K = K
        self.n = nbrOfSplits # matches matrix approach well
        self.Lsp = L

        self.qx = Symbol('qx')
        self.px = Symbol('px')

        xf = self.K**6*self.L**4*self.qx**3/12 - self.K**3*self.L**4*self.px**2/12 - self.K**3*self.L**3*self.qx*self.px/3 - self.L**2*self.K**3*self.qx**2/2 + self.L*self.px + self.qx
        xpf = self.K**6*self.L**4*self.qx**2*self.px/6 + self.K**6*self.L**3*self.qx**3/3 - self.K**3*self.L**3*self.px**2/3 - self.K**3*self.L**2*self.qx*self.px + self.K**6*self.L**4*self.qx**2*self.px/6 + self.K**6*self.L**4*self.qx**2*self.px/12 - self.K**3*self.L*self.qx**2 + self.px

        self.xNumFun = lambdify((self.qx,self.px),xf, "numpy")
        self.xpNumFun = lambdify((self.qx,self.px),xpf, "numpy")

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def printInfo(self):
        return self.name + "\t L: " +  str(self.L) + "\t K: " +  str(self.K)

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n

        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        self.sc.updateMatrix(multipart,twiss)
        multipart, envelope = self.sc.evaluateSC(multipart,envelope) # evaluate the SC # not needed since it is linear
        multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
        #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
        #twiss[4] = envelope[3] / twiss[5]
        #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
        #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
        #twiss[4] = envelope[3] / twiss[5]
        #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluate(self,multipart,envelope,twiss):
        multipart, envelope, env_with_s = self.evaluateNumFun(multipart,envelope) # use the new data for "normal" evaluation
        #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
        #twiss[4] = envelope[3] / twiss[5]
        #twiss[7] = envelope[6] / twiss[8]
        return multipart, envelope, twiss, env_with_s

    def evaluateNumFun(self,multipart,envelope):
        for particle in multipart:

            particle[0][0] = self.xNumFun(particle[0][0], particle[0][1])
            particle[0][1] = self.xpNumFun(particle[0][0], particle[0][1])

            particle[1] += self.Lsp
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

### Rotation !!!!
class Rotation(LinearElement):
    def __init__(self, name, nu_x, nu_y):
        LinearElement.__init__(self, "rotation " + name)        
        self.L = 0.0
        self.M = self.createMatrixM(nu_x, nu_y) # M should be a 6x6 matrix
        self.T = self.createMatrixT(self.M) # M should be a 9x9 matrix
        
        self.Lsp = 0.0
        self.Msp = self.M
        self.Tsp = self.T

    def printInfo(self):
        return self.name + "\t L: " + str(self.L)

    def printMatrix(self):
        MspAsMatrix = np.asmatrix(self.Msp)
        return str(MspAsMatrix**self.n)

    def createMatrixM(self,nu_x, nu_y):
        return np.array([
                [cos(2*constants.pi*nu_x),sin(2*constants.pi*nu_x),0,0,0,0],
                [-sin(2*constants.pi*nu_x),cos(2*constants.pi*nu_x),0,0,0,0],
                [0,0,cos(2*constants.pi*nu_y),sin(2*constants.pi*nu_y),0,0],
                [0,0,-sin(2*constants.pi*nu_y),cos(2*constants.pi*nu_y),0,0],
                [0,0,0,0,1,0],
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

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        return

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
        #twiss[4] = envelope[3] / twiss[5]
        #twiss[7] = envelope[6] / twiss[8]
        multipart, envelope, env_with_s = self.evaluateMT(multipart,envelope) # use the new data for "normal" evaluation
            
        return multipart, envelope, twiss, env_with_s

    def evaluateMT(self,multipart,envelope):
        # should just go through a disunited part
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j][0][0:6] = np.dot(self.Msp, multipart[j][0][0:6])
            multipart[j][1] = multipart[j][1] + self.Lsp
            # remove particles far away
            ####if sqrt(multipart[j][0][0]**2 + multipart[j][0][1]**2 + multipart[j][0][2]**2 + multipart[j][0][3]**2 + multipart[j][0][4]**2 + multipart[j][0][5]**2) > 1:
            ####    print "Particle far away:"
            ####    print "x: " + str(multipart[j][0][0])
            ####    print "xp: " + str(multipart[j][0][1])
            ####    print "y: " + str(multipart[j][0][2])
            ####    print "yp: " + str(multipart[j][0][3])
            ####    print "z: " + str(multipart[j][0][4])
            ####    print "zp: " + str(multipart[j][0][5])
            ####    print "radius: " + str(sqrt(multipart[j][0][0]**2 + multipart[j][0][1]**2 + multipart[j][0][2]**2 + multipart[j][0][3]**2 + multipart[j][0][4]**2 + multipart[j][0][5]**2))
            ####    multipart[j][0][0:6] = np.array([0,0,0,0,0,0])
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s





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
    def __init__(self, name, L, Ezofs, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):#, K, M):
        Element.__init__(self, "cavity " + name, 0) # linear set to zero
        #self.name = name
        #self.K = K
        #self.L = L
        #self.M = self.createMatrixM(K, L) # M should be a 6x6 matrix
        #self.T = self.createMatrixT(K, L) # M should be a 9x9 matrix

        ## Input
        self.L = L
        self.n = nbrOfSplits
        self.spaceChargeOn = spaceChargeOn
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

        self.Lsp = self.L/self.n

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

        # Space Charge!
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('cavity_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('cavity_sc', self.Lsp, multipart, twiss, beamdata)

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

    def updateSC(self, spaceChargeOn, nbrOfSplits, multipart, twiss, beamdata):
        self.n = nbrOfSplits
        self.Lsp = self.L/self.n

        # WARNING THE ELEMENT IS NEVER SPLIT
        
        self.spaceChargeOn = spaceChargeOn
        if self.spaceChargeOn == 1:
            self.sc = SpaceCharge('liealg_sc', self.Lsp, multipart, twiss, beamdata)
        elif self.spaceChargeOn == 2:
            self.sc = SpaceChargeEllipticalIntegral('liealg_sc', self.Lsp, multipart, twiss, beamdata)

    def evaluateWithSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        for i in range(0,self.n):
            self.sc.updateMatrix(multipart,twiss)
            multipart, envelope = self.sc.evaluateSC(multipart,envelope)
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateWithoutSC(self,multipart,envelope,twiss):
        # some for loop that goes through all of the disunited parts
        # WARNING THE ELEMENT IS NEVER SPLIT SO THAT IT WILL BE n NORMAL SIZED ELEMENT IN THE PLACE OF THE WHOLE ELEMENT
        for i in range(0,self.n):
            multipart, envelope, env_with_s = self.evaluateM(multipart,envelope) # use the new data for "normal" evaluation
            #twiss[1] = envelope[0] / twiss[2] # updating beta: beta = sigma**2/epsilon (envelope[0] is sigma_x**2)
            #twiss[4] = envelope[3] / twiss[5]
            #print "twiss[7] before: " + str(twiss[7]) + " \t name: " + self.name
            #twiss[7] = envelope[6] / twiss[8]
            #print "twiss[7] after: " + str(twiss[7])
        return multipart, envelope, twiss, env_with_s

    def evaluateM(self,multipart,envelope):
        # should just go through a disunited part
        # each loop iteration is for a new particle
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Msp, multipart[j][0][0:6]), multipart[j][1] + self.Lsp])
        #envelope = np.dot(self.Tsp, envelope)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.Lsp])
        return multipart, envelope, env_with_s

# Comes from ref E.
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

# Cavity a la Wolski ch 3.6
class CavityMatrix(Element):
    def __init__(self, name, L, a, Efield_0, phi_0, beamdata):
        Element.__init__(self, "cavmat " + name, 1)

        #beta_0 = beamdata[0]
        #gamma_0 = gammaFromBeta(beta_0)
        self.m_0 = beamdata[2]
        self.q = beamdata[3]
        #E_0 = beamdata[4]
    
        #P_0 = self.momFromE(E_0,m_0)
        
        self.p_01 = 2.405
        self.a = a # radius of cavity
        self.Efield_0 = Efield_0 # Electric field amplitude

        self.L = L # Length of cavity
        self.k = self.p_01/self.a # wavenumber for the TM_010 mode

        self.phi_0 = phi_0
        # What can be set has been set the rest depends on the incoming beam

    def printInfo(self):
        return self.name + "\t L: " + str(self.L)

    def EFromMom(self, P, m_0): # Takes momentum and converts it into energy
        return sqrt((P*constants.c)**2+(m_0*constants.c**2)**2)

    def momFromE(self, E, m_0): # Takes energy and converts it into momentum
        return 1/constants.c*sqrt(E**2-(m_0*constants.c**2)**2)

    def elementFromBeamE(self, beta_0, E_0):
        #beta_0 = beamdata[0]
        gamma_0 = gammaFromBeta(beta_0)
        P_0 = self.momFromE(E_0,self.m_0)
        T = 2*constants.pi*beta_0/self.k**2/self.L**2*sin(self.k*self.L/2/beta_0) # Transit-time factor
        #update everything

        V_0 = self.L*self.Efield_0*T # Cavity voltage
        alpha = self.q*V_0/P_0/constants.c # just a definition
        print "alpha: " + str(alpha)
        
    
        wt = self.k*sqrt(alpha*cos(self.phi_0)/2/constants.pi)# omega transversal
        wp = self.k/(beta_0*gamma_0)*sqrt(alpha*cos(self.phi_0)/constants.pi)# omega parallel
    
        ct = cos(wt*self.L) # cos transversal
        cp = cos(wp*self.L) # cos parallel
        st = sin(wt*self.L)/wt # sin transversal
        sp = sin(wp*self.L)/wp # sin parallel
        # the transfer matrices
        self.Rrf = np.array([
            [ct, st, 0, 0, 0, 0],
            [-wt**2*st, ct, 0, 0, 0, 0],
            [0, 0, ct, st, 0, 0],
            [0, 0, -wt**2*st, ct, 0, 0],
            [0, 0, 0, 0, cp, 1/(beta_0**2*gamma_0**2)*sp],
            [0, 0, 0, 0, -beta_0**2*gamma_0**2*wp**2*sp, cp] # Deltadelta!!! (see eqn 2.54 and 2.55 in wolski)
            ])
        # same
        self.mrf = np.array([
            [0],
            [0],
            [0],
            [0],
            [(1-cos(wp*self.L))*tan(self.phi_0)/self.k],
            [beta_0**2*gamma_0**2*wp*sin(wp*self.L)*tan(self.phi_0)/self.k] # Deltadelta!!! (see eqn 2.54 and 2.55 in wolski)
            ])
        print "beta_0**2*gamma_0**2: " + str(beta_0**2*gamma_0**2)
        print "wp: " + str(wp)
        print "sin(wp*self.L): " + str(sin(wp*self.L))
        print "tan(self.phi_0): " + str(tan(self.phi_0))
        print "self.k: " + str(self.k)
    
        # change in reference momentum
        Deltadelta = self.q*V_0/P_0/constants.c*self.k*self.L/constants.pi*sin(self.phi_0) # eqn (3.152)
        E_1 = E_0 + Deltadelta
        # if k*L=pi: Deltadelta = q*V_0/P_0/constants.c*sin(phi_0)
        #z#P_1 = self.momFromE(self.EFromMom(P_0, self.m_0) + Deltadelta, self.m_0)
        #P_1 = self.momFromE(E_0 + Deltadelta, self.m_0)
        P_1 = self.momFromE(E_1, self.m_0)
        print "P_0:" + str(P_0)
        print "P_1:" + str(P_1)
        beta_1 = betaFromE(self.m_0,E_1)
        gamma_1 = gammaFromBeta(beta_1)
        print "gamma_0:" + str(gamma_0)
        print "gamma_1:" + str(gamma_1)
        self.R_deltaP = np.array([
            [1, 0, 0, 0, 0, 0],
            [0, P_0/P_1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, P_0/P_1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, P_0/P_1] # Deltadelta!!! (see eqn 2.54 and 2.55 in wolski)
            ])
        # same
        self.m_deltaP = np.array([
            [0],
            [0],
            [0],
            [0],
            [0],
            [1/(beta_1*constants.c)*(gamma_0/gamma_1-1)] # Deltadelta!!! (see eqn 2.54 and 2.55 in wolski) # The divide by c isn't in wolski and doesn't give right dimensions but it makes my results work...
            ])
        self.newbeta = beta_1

    def evaluateWithoutSC(self,multipart,envelope,twiss, beamdata): # beamdata can be made seperately since energy change is stored in z'
        self.elementFromBeamE(beamdata[0],beamdata[4]) # Update the T and all else that follows...
        print "self.Rrf: \n" + str(self.Rrf)
        print "self.mrf: \n" + str(self.mrf)
        print "self.R_deltaP: \n" + str(self.R_deltaP)
        print "self.m_deltaP: \n" + str(self.m_deltaP)
        print "Total matrix: \n" + str(np.dot(self.R_deltaP, self.Rrf))
        print "Total +m: \n" + str(np.dot(self.R_deltaP, self.mrf) + self.m_deltaP)
        for j in range(0,len(np.atleast_1d(multipart))):
            multipart[j] = np.array([np.dot(self.Rrf, multipart[j][0][0:6]) + self.mrf, multipart[j][1] + self.L])
            multipart[j] = np.array([np.dot(self.R_deltaP, multipart[j][0][0:6]) + self.m_deltaP, multipart[j][1]])
        beamdata[0] = self.newbeta # new beta 
        beamdata[4] = EFromBeta(beamdata[2],beamdata[0]) # new E: E(m_0,beta)
        envelope = envelopeFromMultipart(multipart)
        env_with_s = np.array([envelope, self.L])
        return multipart, envelope, twiss, env_with_s, beamdata

# references
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition
# A. 7.2. Space Charge Impulses in simulatingbeamswithellipsoidalsymmetry-secondedition
# B. A MODIFIED QUADSCAN TECHNIQUE FOR EMITTANCE.pdf
# C. Accelerator-Recipies.pdf by E. Laface
# D. The leapfrog method and other symplectic algorithms for integrating Newtons laws of motion Peter Young Dated April 21 2014
# E. ESS Linac simulator
# F. WEPEA 040