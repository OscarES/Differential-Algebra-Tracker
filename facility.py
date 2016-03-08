from accelerator import *
from scipy import constants
from relativity import betaFromE
from particleFactory import gaussianTwiss3D
import numpy as np




# code that will be the middle man between accelerator and qtinterface

class Facility():
    def __init__(self):
        

        

        ## Beam properties
        self.beamdata = self.getDefaultBeamdata()
        self.twiss = self.getDefaultTwiss()
        self.multipart = gaussianTwiss3D(self.beamdata[5], self.twiss)

        self.lattice = Lattice("Facility", self.beamdata, self.twiss, self.multipart)

        ## Start with just a drift
        self.getDefaultLatticeElement()

    def createDrift(self, name, L):
        self.lattice.createDrift(name, L)
        #drift = Drift(name, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        #self.lattice.appendElement(drift)

    def createDipole(self, name, rho, alpha, n):
        self.lattice.createDipole(name, rho, alpha, n)

    def createQuad(self, name, K, L):
        self.lattice.createQuad(name, K, L)
        #quad = Quad(name, K, L, self.spaceChargeOn, self.multipart, self.twiss, self.beamdata, self.nbrOfSplits)
        #self.lattice.appendElement(quad)

    def createSextupole(self,name, K, L, compOrder):
        self.lattice.createSextupole(name, K, L, compOrder)



    ## Setup
    def getDefaultLatticeElement(self):
        self.lattice.createDrift("Null", 1)

    def getDefaultBeamdata(self):
        E = 2e9*constants.e # 2GeV to joule from ref F.
        freq = 704.42e6 # (Hz) from ref. F
        
        rf_lambda = float(constants.c)/freq  # beam data needed
        m = float(constants.m_p)
        beta = betaFromE(m, E)
        q = float(constants.e)
        nbrOfParticles = 1000
        beamdata = [beta, rf_lambda, m, q, E, nbrOfParticles]
        return beamdata

    def getDefaultTwiss(self):
        twiss = np.array([0.0, 10.3338028723, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06, 0.0, 10.3338028723, 1e-06])
        return twiss
    ## End Setup

    ## Passing of varibles to and fro
    def getLattice(self):
        return self.lattice.getLattice()
        #return self.lattice

    def getBeamdata(self):
        return self.lattice.getBeamdata()
        #return self.beamdata

    def getTwiss(self):
        return self.lattice.getTwiss()
        #return self.twiss

    def getMultipart(self):
        return self.lattice.getMultipart()
        #return self.multipart

    def getSpaceChargeOn(self):
        return self.lattice.getSpaceChargeOn()
        #return self.spaceChargeOn

    def getNbrOfSplits(self):
        return self.lattice.getNbrOfSplits()
        #return self.nbrOfSplits

    def setLattice(self,lattice): # changes the entire lattice object!
        #self.lattice.setLattice(lattice)
        self.lattice = lattice

    def setBeamdata(self,beamdata):
        self.lattice.setBeamdata(beamdata)
        #self.beamdata = beamdata

    def setTwiss(self,twiss):
        self.lattice.setTwiss(twiss)
        #self.twiss = twiss

    def setMultipart(self,multipart):
        self.lattice.setMultipart(multipart)
        #self.multipart = multipart

    def setSpaceChargeOn(self, spaceChargeOn):
        self.lattice.setSpaceChargeOn(spaceChargeOn)
        #self.spaceChargeOn = spaceChargeOn
    ## End Passing

    def printLattice(self):
        return self.lattice.printLattice()

    def evaluate(self, multipart, envelope, twiss):
        resultmultipart, resultenvelope, resulttwiss = self.lattice.evaluate(multipart,envelope,twiss)
        return resultmultipart, resultenvelope, resulttwiss