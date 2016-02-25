from accelerator import *
from scipy import constants
from relativity import betaFromE




# code that will be the middle man between accelerator and qtinterface

class Facility():
    def __init__(self):
        self.lattice = Lattice("Facility")

        # really ugly solution too avoid an empty lattice
        drift = Drift("Null", 1, 0, 0, 0, 0, 1)
        self.lattice.appendElement(drift)

    def createDrift(self, name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        drift = Drift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        self.lattice.appendElement(drift)

    def createQuad(self, name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        quad = Quad(name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        self.lattice.appendElement(quad)

    def getLattice(self):
        return self.lattice

    def setLattice(self,lattice):
        self.lattice = lattice

def getBeamdata():
        E = 2e9*constants.e # 2GeV to joule from ref F.
        freq = 704.42e6 # (Hz) from ref. F
        
        rf_lambda = constants.c/freq  # beam data needed
        m = constants.m_p
        beta = betaFromE(m, E)
        q = constants.e
        beamdata = [beta, rf_lambda, m, q]
        return beamdata