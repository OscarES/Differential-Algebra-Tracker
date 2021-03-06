from accelerator import *
from scipy import constants
from relativity import betaFromE
from particleFactory import gaussianTwiss3D, envelopeFromMultipart, straightxxp, singleparticle
import numpy as np
from plotting import plotEverything, plot_x_before_and_after, plot_x_after, plot_z_before_and_after, plot_x_with_mulpart_and_twiss
from profilehooks import profile


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

    def createLieDrift(self, name, L, compOrder):
        self.lattice.createLieDrift(name, L, compOrder)

    def createSextupole(self, name, K, L, compOrder):
        self.lattice.createSextupole(name, K, L, compOrder)

    def createOctupole(self, name, K, L, compOrder):
        self.lattice.createOctupole(name, K, L, compOrder)

    def createSextupolerel(self, name, K, L, compOrder):
        self.lattice.createSextupolerel(name, K, L, compOrder)

    def createSextupolekin(self, name, K, L, compOrder):
        self.lattice.createSextupolekin(name, K, L, compOrder)

    def createOctupolekin(self, name, K, L, compOrder):
        self.lattice.createOctupolekin(name, K, L, compOrder)

    def createSextupolemat(self, name, K, L):
        self.lattice.createSextupolemat(name, K, L)

    def createSextupolematema(self, name, K, L):
        self.lattice.createSextupolematema(name, K, L)

    def createRotation(self, name, nu_x, nu_y):
        self.lattice.createRotation(name, nu_x, nu_y)

    def createCavity(self, name, L, Ezofs):
        self.lattice.createCavity(name, L, Ezofs)

    def createCavityMatrix(self, name, L, a, Efield_0, phi_0):
        self.lattice.createCavityMatrix(name, L, a, Efield_0, phi_0)

    ## Setup
    def getDefaultBeamdata(self):
        E = 2e9*constants.e # 2GeV to joule from ref F. # Why show this in Joules?
        freq = 704.42e6 # (Hz) from ref. F
        
        rf_lambda = float(constants.c)/freq  # beam data needed
        m = float(constants.m_p)
        beta = betaFromE(m, E)
        q = float(constants.e)
        nbrOfParticles = 1000
        I = 0.0625
        beamdata = [beta, rf_lambda, m, q, E, nbrOfParticles, I]
        return beamdata

    def getDefaultTwiss(self):
        twiss = np.array([0.0, 10.3338028723, 3.14159e-06, -3.331460652e-16, 8.85901414121, 3.14159e-06, 0.0, 10.3338028723, 1e-06])
        return twiss
    ## End Setup

    ## Passing of varibles to and fro
    def getLattice(self):
        return self.lattice

    def getLatticeObj(self):
        return self.lattice

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

    def deleteLattice(self):
        # empty lattice
        self.lattice = Lattice("Facility", self.beamdata, self.twiss, self.multipart)

    def setBeamdata(self,beamdata):
        self.lattice.setBeamdata(beamdata)
        self.beamdata = self.lattice.getBeamdata()

    def setTwiss(self,twiss):
        self.lattice.setTwiss(twiss)
        self.twiss = self.lattice.getTwiss()

    def generateMultipart(self, nbrOfParticles, twiss):
        multipart = gaussianTwiss3D(nbrOfParticles, twiss)
        self.setMultipart(multipart)

    def generateGridpart(self, nbrOfParticles):
        multipart = straightxxp(nbrOfParticles) #### !!!!! wolski
        self.setMultipart(multipart)

    def setMultipart(self,multipart): # Also updates the envelope
        self.lattice.setMultipart(multipart)
        self.multipart = self.lattice.getMultipart()
        self.envelope = envelopeFromMultipart(self.multipart)

    def setSpaceChargeOnAndSplits(self, spaceChargeOn, nbrOfSplits):
        self.lattice.setSpaceChargeOnAndSplits(spaceChargeOn, nbrOfSplits)

    def getLaps(self):
        return self.lattice.getLaps()

    def setLaps(self, laps):
        self.lattice.setLaps(laps)
    ## End Passing

    def printLattice(self):
        return self.lattice.printLattice()

    def printMatrices(self):
        return self.lattice.printMatrices()

    def evaluate(self):
        self.resultmultipart, self.resultenvelope, self.resulttwiss, self.resultenvlist, self.multipartafterall, self.resultrealtwiss = self.lattice.evaluate(self.multipart,self.envelope,self.twiss)

    @profile(immediate=true)
    def evaluateWithProfiling(self):
        self.evaluate()

    def getResults(self):
        #if 'self.resultmultipart' in locals():
        if hasattr(self,'resultmultipart'):
            return self.resultmultipart, self.resultenvelope, self.resulttwiss, self.resultenvlist
        else:
            return

    def plotAfterEval(self):
        plotEverything(self.multipart, self.twiss, self.resultmultipart, self.resultenvlist, self.multipartafterall)
        plot_x_before_and_after(self.multipart, self.resultmultipart)
        plot_x_after(self.resultmultipart)
        plot_z_before_and_after(self.multipart, self.resultmultipart)

        # plot mulpart and twiss
        realtwissin = np.array([self.twiss[1], self.twiss[0], (1+self.twiss[0]**2)/self.twiss[1], self.twiss[4], self.twiss[3], (1+self.twiss[3]**2)/self.twiss[4], self.twiss[7], self.twiss[6], (1+self.twiss[6]**2)/self.twiss[7]])
        emittance_x = self.twiss[2]
        plot_x_with_mulpart_and_twiss(self.multipart, realtwissin, self.resultmultipart, self.resultrealtwiss, emittance_x)