from accelerator import *
from scipy import constants
from relativity import betaFromE
from particleFactory import gaussianTwiss3D
import numpy as np




# code that will be the middle man between accelerator and qtinterface

class Facility():
    def __init__(self):
        self.lattice = Lattice("Facility")

        self.getDefaultLatticeElement()

        ## Beam properties
        self.beamdata = self.getDefaultBeamdata()
        self.twiss = self.getDefaultTwiss()
        self.multipart = gaussianTwiss3D(self.beamdata[5], self.twiss)

    def createDrift(self, name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        drift = Drift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        self.lattice.appendElement(drift)

    def createQuad(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        quad = Quad(name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        self.lattice.appendElement(quad)

    def getDefaultLatticeElement(self):
        drift = Drift("Null", 1, 0, 0, 0, 0, 1)
        self.lattice.appendElement(drift)

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

    def getLattice(self):
        return self.lattice

    def getBeamdata(self):
        return self.beamdata

    def getTwiss(self):
        return self.twiss

    def getMultipart(self):
        return self.multipart

    def setLattice(self,lattice):
        self.lattice = lattice

    def setBeamdata(self,beamdata):
        self.beamdata = beamdata

    def setTwiss(self,twiss):
        self.twiss = twiss

    def setMultipart(self,multipart):
        self.multipart = multipart

def getBeamdata():
        E = 2e9*constants.e # 2GeV to joule from ref F.
        freq = 704.42e6 # (Hz) from ref. F
        
        rf_lambda = float(constants.c)/freq  # beam data needed
        m = float(constants.m_p)
        beta = betaFromE(m, E)
        q = float(constants.e)
        nbrOfParticles = 1000
        beamdata = [beta, rf_lambda, m, q, E, nbrOfParticles]
        return beamdata