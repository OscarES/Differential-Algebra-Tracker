import numpy as np
import pickle
import os

def saveMultipart(filename, multipart):
    np.save(filename, multipart)
    return 1

# OK

def loadMultipart(filename):
    multipart = np.load(filename)
    return multipart

# OK

def saveTwiss(filename, twiss):
    np.save(filename, twiss)
    return 1

def loadTwiss(filename):
    twiss = np.load(filename)
    return twiss

def saveEnvelope(filename, envelope):
    np.save(filename, envelope)
    return 1

def loadEnvelope(filename):
    envelope = np.load(filename)
    return envelope

### NEEDS fixin
def saveLattice(filename, lattice):
    #np.save(filename, lattice)
    if not os.path.isfile(filename):
        os.mknod(filename)
    pickle.dump(lattice, open(filename, 'rwb'))
    return 1

def loadLattice(filename):
    #lattice = np.load(filename)
    lattice = pickle.load(filename)
    return lattice