from accelerator import *
from scipy import constants
from relativity import betaFromE
from particleFactory import gaussianTwiss3D, envelopeFromMultipart
import numpy as np
from plotting import plotEverything


# code that will be the middle man between accelerator and qtinterface

class Facility():
    def __init__(self):
        ## Beam properties
        self.beamdata = self.getDefaultBeamdata()
        self.twiss = self.getDefaultTwiss()
        self.multipart = gaussianTwiss3D(self.beamdata[5], self.twiss)
        self.envelope = envelopeFromMultipart(self.multipart)

        # empty lattice
        self.lattice = Lattice("Facility", self.beamdata, self.twiss, self.multipart)

    def createDrift(self, name, L):
        self.lattice.createDrift(name, L)

    def createDipole(self, name, rho, alpha, n):
        self.lattice.createDipole(name, rho, alpha, n)

    def createQuadrupole(self, name, K, L):
        self.lattice.createQuadrupole(name, K, L)

    def createSextupole(self,name, K, L, compOrder):
        self.lattice.createSextupole(name, K, L, compOrder)

    def createCavity(self, name, L, Ezofs):
        self.lattice.createCavity(name, L, Ezofs)

    ## Setup
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

    def getBeamdata(self):
        return self.lattice.getBeamdata()

    def getTwiss(self):
        return self.lattice.getTwiss()

    def getMultipart(self):
        return self.lattice.getMultipart()

    def getSpaceChargeOn(self):
        return self.lattice.getSpaceChargeOn()

    def getNbrOfSplits(self):
        return self.lattice.getNbrOfSplits()

    def setLattice(self,lattice): # changes the entire lattice object!
        self.lattice = lattice

    def setBeamdata(self,beamdata):
        self.lattice.setBeamdata(beamdata)
        self.beamdata = self.lattice.getBeamdata()

    def setTwiss(self,twiss):
        self.lattice.setTwiss(twiss)
        self.twiss = self.lattice.getTwiss()

    def generateMultipart(self, nbrOfParticles, twiss):
        multipart = gaussianTwiss3D(nbrOfParticles, twiss)
        self.setMultipart(multipart)

    def setMultipart(self,multipart):
        self.lattice.setMultipart(multipart)
        self.multipart = self.lattice.getMultipart()

    def setSpaceChargeOn(self, spaceChargeOn):
        self.lattice.setSpaceChargeOn(spaceChargeOn)
    ## End Passing

    def printLattice(self):
        return self.lattice.printLattice()

    def evaluate(self):
        self.resultmultipart, self.resultenvelope, self.resulttwiss, self.resultenvlist = self.lattice.evaluate(self.multipart,self.envelope,self.twiss)
        #return self.resultmultipart, self.resultenvelope, self.resulttwiss, self.resultenvlist

    def plotAfterEval(self):
        plotEverything(self.multipart, self.twiss, self.resultmultipart, self.resultenvlist)