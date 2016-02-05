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

# Note that this doesn't have z and zp, or envelope or lattice. Just x, xp, y, yp and the twiss
def saveSummer2015Format(filenamepart, filenametwiss, multipart, twiss):
    dta = np.dtype([('x', 'd'), ('xp', 'd'), ('y', 'd'), ('yp', 'd')])
    dtb = np.dtype([('alpha_x', 'd'), ('beta_x', 'd'), ('epsilon_x', 'd'), ('alpha_y', 'd'), ('beta_y', 'd'), ('epsilon_y', 'd')])    

    x = [multipart[i][0][0] for i in xrange(len(multipart))]
    xp = [multipart[i][0][1] for i in xrange(len(multipart))]
    y = [multipart[i][0][2] for i in xrange(len(multipart))]
    yp = [multipart[i][0][3] for i in xrange(len(multipart))]

    nbrOfParticles = len(x)

    alpha_x = twiss[0]
    beta_x = twiss[1]
    epsilon_rms_x = twiss[2]
    alpha_y = twiss[3]
    beta_y = twiss[4]
    epsilon_rms_y = twiss[5]
    alpha_z = twiss[6]

    epsilon_x = epsilon_rms_x
    epsilon_y = epsilon_rms_y

    a = np.zeros(nbrOfParticles, dta)
    b = np.zeros(1,dtb)

    a['x'] = x
    a['xp'] = xp
    a['y'] = y
    a['yp'] = yp
    b['alpha_x'] = alpha_x
    b['beta_x'] = beta_x
    b['epsilon_x'] = epsilon_x
    b['alpha_y'] = alpha_y
    b['beta_y'] = beta_y
    b['epsilon_y'] = epsilon_y

    np.savetxt(filenamepart, a, '%10s')
    np.savetxt(filenametwiss, b, '%10s')

    return 1

# Note that this don't have z and zp
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