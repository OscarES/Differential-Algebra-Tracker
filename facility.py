from accelerator import *




# code that will be the middle man between accelerator and qtinterface

class Facility():
    def __init__(self):
        self.lattice = Lattice("Facility")

    def createDrift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        drift = Drift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        self.lattice.appendElement(drift)