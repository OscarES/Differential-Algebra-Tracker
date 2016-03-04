import numpy as np
import pickle
import os
import copy
from accelerator import Lattice, Element, LinearElement, Quad, Drift, LieAlgebra, LieAlgElement, leapfrog, Dipole, Cavity

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

def saveBeamdata(filename, beamdata):
    np.save(filename, beamdata)
    return 1

def loadBeamdata(filename):
    try:
        beamdata = np.load(filename)
        return beamdata
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

#Can't handle LieAlgElements!!!!!!
#def saveLattice(filename, lattice):
#    if not os.path.isfile(filename):
#        os.mknod(filename)
#    pickle.dump(lattice, open(filename, 'wb'))
#    return 1
#
#def loadLattice(filename):
#    try:
#        lattice = pickle.load(open(filename, 'rb'))
#        return lattice
#    except:
#        print 'Bad datafile!'
#        quit()
#        return 0

# Can handle sextupole elements!
# creates a lattice equal to that described by text (which is output from a printLattice call)
def parseLatticeString(text, facility):
    spaceChargeOn = facility.getSpaceChargeOn()
    multipart = facility.getMultipart()
    twiss = facility.getTwiss()
    beamdata = facility.getBeamdata()
    nbrOfSplits = facility.getNbrOfSplits()

    lattice = Lattice('ParsedLattice')
    for line in iter(text.splitlines()):
        words = line.split()
        typeOfElem = words[0]
        name = words[1]
        l = float(words[words.index("L:") + 1]) #what comes after "L:"
        if typeOfElem == "cavity":
            cavityOscillations = float(words[words.index("Oscillations:") + 1])
            cavityAmplitudeA = float(words[words.index("AmplitudeA:") + 1])
            cavityAmplitudeB = float(words[words.index("AmplitudeB:") + 1])
            cavityE_0 = float(words[words.index("E_0:") + 1])
            cavitySigma = float(words[words.index("sigma:") + 1])
            cavityP = float(words[words.index("p:") + 1])

            cavityEzofs = [cavityOscillations, cavityAmplitudeA, cavityAmplitudeB, cavityE_0, cavitySigma, cavityP]
            elem = Cavity(name, l, cavityEzofs, beamdata, nbrOfSplits)
        if typeOfElem == "dipole":
            rho = float(words[words.index("rho:") + 1]) #what comes after "rho:"
            #k_x = what comes after "K_x: "
            #k_y = what comes after "K_y: " # not needed for construction
            beta = beamdata[0]
            nparam = float(words[words.index("nparam:") + 1]) #what comes after "nparam:"
            alpha = float(words[words.index("Alpha:") + 1]) #what comes after "Alpha:"
            elem = Dipole(name, rho, alpha, nparam, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        elif typeOfElem != "drift" and typeOfElem != "cavity":
            k = float(words[words.index("K:") + 1]) #what comes after "K:"
        if typeOfElem == "liealgelem":
            hamToUse = words[words.index("HamUsed:") + 1] #what comes after "HamUsed:" and before next whitespace
            order = int(words[words.index("Order:") + 1]) #what comes after "Order:"
            elem = LieAlgElement(name, hamToUse, k, l, order, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        
        if typeOfElem == "quad":
            elem = Quad(name, k, l, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        if typeOfElem == "drift":
            elem = Drift(name, l, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        lattice.appendElement(elem)           
    return lattice

def saveLattice(filename, lattice):
    latticeString = lattice.printLattice()
    np.save(filename, latticeString)
    return 1

def loadLattice(filename, facility):
    try:
        latticeString = str(np.load(filename))
    except:
        print 'Bad datafile!'
        return 0
    try:
        lattice = parseLatticeString(latticeString, facility)
        return lattice
    except:
        print 'Error while parsing lattice string!'
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

def loadSummer2015Formatzasx(datafilepart, datafiletwiss):
    try:
        x, xp, y, yp = np.loadtxt(datafilepart,unpack = True)
        alpha_x, beta_x, epsilon_x, alpha_y, beta_y, epsilon_y = np.loadtxt(datafiletwiss,unpack = True)

        nbrOfParticles = len(x)

        z = copy.deepcopy(x)
        zp = copy.deepcopy(xp)

        s = np.linspace(0,0,nbrOfParticles)

        bigMatrix = np.array([x, xp, y, yp, z, zp])

        multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

        epsilon_rms_x = epsilon_x
        epsilon_rms_y = epsilon_y

        alpha_z = copy.deepcopy(alpha_x)
        beta_z = copy.deepcopy(beta_x)
        epsilon_rms_z = copy.deepcopy(epsilon_rms_x)

        twiss = np.array([alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z])

        return multipart, twiss
    except:
        print 'Bad datafile!'
        quit()
        return 0