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
    try:
        multipart = loadMultipart(filename + "multipart" + ".npy")
        twiss = loadTwiss(filename + "twiss" + ".npy")
        envelope = loadEnvelope(filename + "envelope" + ".npy")
        lattice = loadLattice(filename + "lattice" + ".npy")
        return multipart, twiss, envelope, lattice
    except:
        print 'Bad datafile!'
        quit()
        return 0

def saveMultipart(filename, multipart):
    np.save(filename, multipart)
    return 1

def loadMultipart(filename):
    try:
        multipart = np.load(filename)
        return multipart
    except:
        print 'Bad datafile!'
        quit()
        return 0

def saveTwiss(filename, twiss):
    np.save(filename, twiss)
    return 1

def loadTwiss(filename):
    try:
        twiss = np.load(filename)
        return twiss
    except:
        print 'Bad datafile!'
        quit()
        return 0

def saveEnvelope(filename, envelope):
    np.save(filename, envelope)
    return 1

def loadEnvelope(filename):
    try:
        envelope = np.load(filename)
        return envelope
    except:
        print 'Bad datafile!'
        quit()
        return 0

def saveLattice(filename, lattice):
    #np.save(filename, lattice)
    if not os.path.isfile(filename):
        os.mknod(filename)
    pickle.dump(lattice, open(filename, 'wb'))
    return 1

def loadLattice(filename):
    try:
        lattice = pickle.load(open(filename, 'rb'))
        return lattice
    except:
        print 'Bad datafile!'
        quit()
        return 0

def saveSummer2015Format():
    return 0

def loadSummer2015Format(datafilepart, datafiletwiss):
    try:
        x, xp, y, yp = np.loadtxt(datafilepart,unpack = True)
        alpha_x, beta_x, epsilon_x, alpha_y, beta_y, epsilon_y = np.loadtxt(datafiletwiss,unpack = True)

        nbrOfParticles = len(x)

        z = np.linspace(0,0,nbrOfParticles)
        zp = np.linspace(0,0,nbrOfParticles)

        s = np.linspace(0,0,nbrOfParticles)

        bigMatrix = np.array([x, xp, y, yp, z, zp])

        multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

        epsilon_rms_x = epsilon_x
        epsilon_rms_y = epsilon_y

        alpha_z = 0
        beta_z = 0
        epsilon_rms_z = 0

        twiss = np.array([alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z])

        return multipart, twiss
    except:
        print 'Bad datafile!'
        quit()
        return 0