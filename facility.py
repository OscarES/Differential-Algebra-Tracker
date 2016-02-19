from accelerator import *




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