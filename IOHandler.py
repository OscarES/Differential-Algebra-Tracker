import numpy as np
import pickle
import os

def saveAll(filename, multipart, twiss, envelope, lattice):
    saveMultipart(filename + "multipart", multipart)
    saveTwiss(filename + "twiss", twiss)
    saveEnvelope(filename + "envelope", envelope)
    saveLattice(filename + "lattice" + ".npy", lattice)
    return 1

def loadAll(filename):
    multipart = loadMultipart(filename + "multipart" + ".npy")
    twiss = loadTwiss(filename + "twiss" + ".npy")
    envelope = loadEnvelope(filename + "envelope" + ".npy")
    lattice = loadLattice(filename + "lattice" + ".npy")
    return multipart, twiss, envelope, lattice

def saveMultipart(filename, multipart):
    np.save(filename, multipart)
    return 1

def loadMultipart(filename):
    multipart = np.load(filename)
    return multipart

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

def saveLattice(filename, lattice):
    #np.save(filename, lattice)
    if not os.path.isfile(filename):
        os.mknod(filename)
    pickle.dump(lattice, open(filename, 'wb'))
    return 1

def loadLattice(filename):
    #lattice = np.load(filename)
    lattice = pickle.load(open(filename, 'rb'))
    return lattice

def saveSummer2015Format():
    return 0

def loadSummer2015Format():
    return 0