from __future__ import division # needed for 1/2 = 0.5
import numpy as np
from accelerator import Lattice, Element, LinearElement, Quad, Drift, LieAlgebra, LieAlgElement, leapfrog, Dipole, Cavity
#from sympy.parsing.sympy_parser import parse_expr
#from sympy import *
from particleFactory import straight, scanned, randomed, gaussian, gaussianTwiss3D, envelopeFromMultipart
from plotting import plotEverything, plotEnvelope, plotPhaseSpace
from IOHandler import saveAll, loadAll, saveMultipart, loadMultipart, saveTwiss, loadTwiss, saveEnvelope, loadEnvelope, saveLattice, loadLattice, saveSummer2015Format, loadSummer2015Format, loadSummer2015Formatzasx, parseLatticeString, saveLatticeString, loadLatticeString
import copy
from scipy import constants
from relativity import betaFromE
import math


#sigma = 0.001 # the standard deviation that the user will enter
#epsilon = sqrt(sigmax**2*sigmaxp**2-sigmaxxp**2)

### Non-linear
#order = 5 # user sets this
#
#LA = LieAlgebra()
#
#quaddefocuskval = 0.1
#quaddefocuslval = 1
#quaddefocusNumFuns = LA.hamToNumFuns(quadhamdefocus, quaddefocuskval, quaddefocuslval, order)
#
##print "quaddefocusNumFuns:" + str(quaddefocusNumFuns)
#
##diffAlgRes = (quaddefocusNumFuns[0](particle1), quaddefocusNumFuns[1](particle1), quaddefocusNumFuns[2](particle1), quaddefocusNumFuns[3](particle1), quaddefocusNumFuns[4](particle1), quaddefocusNumFuns[5](particle1))
#
##sextukval = 
##sextupoleNumFuns = DA.hamToNumFuns(sextupoleham, kval, lval, order)
#
### Multiparticle
#x = 0.0
#xp =1.0
#y = 0.0
#yp =1.0
#z = 0.0
#zp =0.0
#s = 1.0 # from ref 1 in section 2.3
#zvector = np.array([x, xp, y, yp, z, zp])
#particle1 = np.array([zvector, s])
#particle2 = np.array([-zvector, s]) # careful, if this is the same expression as for particle1 there will only be one object!
#multipart = np.array([particle1, particle2])
#multipart2 = copy.copy(multipart)
##print "len(np.atleast_1d(multipart))" + str(len(np.atleast_1d(multipart)))
##print "multipart[0][0:6]" + str(multipart[0][0:6])
##print "multipart[1] (s)" + str(multipart[1])
#print "multipart" + str(multipart)
##print "len(multipart) " + str(len(multipart))
#
#
## twiss comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
#twiss = np.array([0.0, 10.3338028723, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06])
#
## this new def needs twiss params perhaps? YES, since otherwise SC won't work
## envelope comes as [sigma_x**2, sigma_x*sigma_xp, sigma_xp**2, sigma_y**2, sigma_y*sigma_yp, sigma_yp**2, sigma_z**2, sigma_z*sigma_zp, sigma_zp**2]
#envelope = np.array([1, 0, 0, 1, 0, 0, 1, 0, 0])
#
#
### Lattice construction
#spaceChargeOn = 1
#drift = Drift('drift', 1, spaceChargeOn, multipart, twiss)
#
## K = sqrt(e*g/p) ,from ref E.
#quad = Quad('quad', 0.1, 1, spaceChargeOn, multipart, twiss)
#
#lattice = Lattice('lattice')
#lattice.appendElement(drift)
#lattice.appendElement(quad)
#print "latticeprint: \n" + lattice.printLattice()
##print "multipart before latt eval: " + str(multipart)
#partres, envres = lattice.evaluate(multipart,envelope) ### changes multipart and envelope!!!!!!!!!!!!!!!!!!
##print "multipart2 after latt eval: " + str(multipart2)
#
#print "partres: " + str(partres)
#print "envres: " + str(envres)
#
## Using the Lie algebra "raw"
#LieAlgRes = (quaddefocusNumFuns[0](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[1](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[2](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[3](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[4](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[5](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]))
##print "DiffAlgRes: " + str(diffAlgRes) # MATCHES the linear calculation!!!!! :)
#
## Using the differential algebra through DiffAlgElement
#LieAlgElemQuad = LieAlgElement("LieAlgElemQuad", LA, quadhamdefocus, quaddefocuskval, quaddefocuslval, order, 0, multipart2, envelope)
#LieAlgElemQuadpartres, LieAlgElemQuadenvres = LieAlgElemQuad.evaluate(multipart2,envelope) # not the same as raw! Since the element is split! But if n = 1 in the class then there is a very good match!!!! :)
#print "LieAlgElemQuadRes: " + str(LieAlgElemQuadpartres)


### Tons of particles
print "Tons of particles..."
#nbrOfParticles = 10
#multiparttonsMat = straight(nbrOfParticles)
#multiparttonsDiffAlg = straight(nbrOfParticles)
#
##print "multiparttonsMat: \n" + str(multiparttonsMat)
#
#tonsMatPartRes, tonsMatEnvRes = lattice.evaluate(multiparttonsMat,envelope)
#tonsLieAlgPartRes, tonsLieAlgEnvRes = LieAlgElemQuad.evaluate(multiparttonsDiffAlg,envelope)
#
##print "tonsMatPartRes: \n" + str(tonsMatPartRes)
##print "tonsLieAlgPartRes: \n" + str(tonsLieAlgPartRes)
#
## Checking the others
#scannedparts = scanned(nbrOfParticles)
##print "scannedparts: \n" + str(scannedparts)
#
#randomedparts = randomed(nbrOfParticles)
##print "randomedparts: \n" + str(randomedparts)
#
#gaussianparts = gaussian(nbrOfParticles)
##print "gaussianparts: \n" + str(gaussianparts)
#
#gaussianTwiss3Dparts = gaussianTwiss3D(nbrOfParticles, twiss)
##print "gaussianTwiss3Dparts: \n" + str(gaussianTwiss3Dparts)

### Leapfrog
print "Leapfrog..."
#x_0 = 100.0
#v_0 = 0.0
#def F(x):
#    return -x
#
## F should perhaps be numfuns: xpNumFun, ypNumFun, zpNumFun
## F is Force!
#
#L = 4
#n = 4*10*10000
#h = L/n
#
#x_of_i, v_of_i = leapfrog(x_0, v_0, F, h, n)
##print "Results:"
#print "x_of_i[-1]: " + str(x_of_i[-1])
#print "v_of_i[-1]: " + str(v_of_i[-1])



### IOHandling
print "IOHandling..."
## multipart
#filenameMultipart = "savedParticles"
#saveMultipart(filenameMultipart, multipart2)
#print "multipart2 stored in " + filenameMultipart
#
#loadedmultipart = loadMultipart(filenameMultipart + ".npy")
#print "loaded particles: \n" + str(loadedmultipart)
#
## twiss
#filenameTwiss = "savedTwiss"
#saveTwiss(filenameTwiss, twiss)
#print "twiss stored in " + filenameTwiss
#
#loadedtwiss = loadMultipart(filenameTwiss + ".npy")
#print "loaded twiss: \n" + str(loadedtwiss)
#
## envelope
#filenameEnvelope = "savedEnvelope"
#saveEnvelope(filenameEnvelope, envelope)
#print "envelope stored in " + filenameEnvelope
#
#loadedenvelope = loadEnvelope(filenameEnvelope + ".npy")
#print "loaded envelope: \n" + str(loadedenvelope)
#
## lattice
#filenameLattice = "savedLattice.npy"
#saveLattice(filenameLattice, lattice)
#print "lattice stored in " + filenameLattice
#
#loadedlattice = loadLattice(filenameLattice)
#print "loaded lattice: \n" + loadedlattice.printLattice()

#filename = "data/" + "saved"
#saveAll(filename, multipart2, twiss, envelope, lattice)
#multipartload, twissload, envelopeload, latticeload = loadAll(filename)
#print "multipartload: \n" + str(multipartload)
#print "twissload: \n" + str(twissload)
#print "envelopeload: \n" + str(envelopeload)
#print "latticeload: \n" + str(latticeload.printLattice())
#

### Compare with old code
print "Compare with old code/ compare space charges..."
## Load the same particles (old code comp)
datafilepart = "../data/" + "inpart1000" + ".txt"
datafiletwiss = "../data/" + "intwiss" + ".txt"
#multipartfromold, twissfromold = loadSummer2015Format(datafilepart, datafiletwiss)
#multipartfromoldcopy = copy.deepcopy(multipartfromold)
#print "multipartfromoldcopy: \n" + str(multipartfromoldcopy)

# sets z as x from old format
multipartfromold, twissfromold = loadSummer2015Formatzasx(datafilepart, datafiletwiss)
multipartfromoldcopy = copy.deepcopy(multipartfromold)
twissfromoldcopy = copy.deepcopy(twissfromold)


#nbrOfParticles = 1000
#multipartfromold = gaussianTwiss3D(nbrOfParticles, twissfromold)
#multipartfromoldcopy = copy.deepcopy(multipartfromold)
#print "multipartfromoldcopy: \n" + str(multipartfromoldcopy)

spaceChargeOnInComp = 0

E = 2e9*constants.e # 2GeV to joule from ref F.
freq = 704.42e6 # (Hz) from ref. F

rf_lambda = constants.c/freq  # beam data needed
m = constants.m_p
beta = betaFromE(m, E)
q = constants.e
beamdata = [beta, rf_lambda, m, q]

#envelopeInComp = np.array([1, 0, 0, 1, 0, 0, 1, 0, 0])
envelopeInComp = envelopeFromMultipart(multipartfromold)
#print "envelopeInComp: " + str(envelopeInComp)

nbrOfSplits = 1

## the lattice will be a FODSO cell (Focusing Quad, Drift, Defocusing Quad, Sextupole, Drift)
compLattice = Lattice('compLattice')

cavityName = "cavity"
cavityLength = 1.0
cavity = Cavity(cavityName, cavityLength)

#fQName = "fQ"
#fQuadLength = 0.4
#fQuadStrength = -0.8 # this is k
#fQ = Quad(fQName, fQuadStrength, fQuadLength, spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
##compLattice.appendElement(fQ)
#
#driftName = "drift"
#driftLength = 1.0
#compDrift = Drift(driftName, driftLength, spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
##compLattice.appendElement(compDrift)
#
#dQName = "dQ"
#dQuadLength = 0.4
#dQuadStrength = 0.8
#dQ = Quad(dQName, dQuadStrength, dQuadLength, spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
##compLattice.appendElement(dQ)
#
#sextuName = "sextu"
#hamToUse = "sextupoleham"
#sextuLength = 0.3
#sextuStrength = 0.6
#compOrder = 6
#sextu = LieAlgElement(sextuName, hamToUse, sextuStrength, sextuLength, compOrder, spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
#compLattice.appendElement(sextu)
#
##compLattice.appendElement(compDrift)
#
#dipoleName = "dipole"
#dipoleRho = 5
#dipoleAlpha = math.pi/4 
#dipolen = 0.5
#compDipole = Dipole(dipoleName, dipoleRho, dipoleAlpha, dipolen, spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
#compLattice.appendElement(compDipole)
#
#compLattice.appendElement(fQ)
#compLattice.appendElement(compDrift)
#compLattice.appendElement(dQ)
#compLattice.appendElement(compDrift)
#
#compLattice.appendElement(fQ)
#compLattice.appendElement(compDrift)
#compLattice.appendElement(dQ)
#compLattice.appendElement(compDrift)
#
#compLattice.appendElement(fQ)
#compLattice.appendElement(compDrift)
#compLattice.appendElement(dQ)
#compLattice.appendElement(compDrift)

#print compLattice.printLattice()
#parsedLattice = parseLatticeString(compLattice.printLattice(), spaceChargeOnInComp, multipartfromold, twissfromold, beamdata, nbrOfSplits)
#print "parsedLattice: \n" + parsedLattice.printLattice()

#saveLatticeString("../data/" + "savedlatticestringonlysextu" + ".npy", parsedLattice)

#loadedLatticeString = loadLatticeString("../data/" + "savedlatticestringwithsextu" + ".npy")
#print "loadedLatticeString: \n" + loadedLatticeString

#saveLattice("../data/" + "savedlatticewithsextu" + ".npy", compLattice)

## Calculate
#partresInComp, envresInComp, twissresInComp = compLattice.evaluate(multipartfromold,envelopeInComp,twissfromold) # Does eval still change input? yes
#partresInComp, envresInComp, twissresInComp = parsedLattice.evaluate(multipartfromold,envelopeInComp,twissfromold) # Does eval still change input? yes

#saveSummer2015Format("../data/" + "outpartFODSOspaceChargetesttest" + ".txt","../data/" + "outtwiss" + ".txt",partresInComp, twissfromold)

#plotEverything(multipartfromoldcopy, twissfromoldcopy, partresInComp)

### Plotting
print "Plotting..."
#twissin = np.array([0.0, 10.3338028723, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06])
#twissout = twissin
#envelopein = np.array([1, 0, 0, 1, 0, 0, 1, 0, 0])
#envelopeout = envelopein
#
#nbrOfParticlestest = 100
#
#multipartin = gaussianTwiss3D(nbrOfParticlestest, twissout)
#multipartout = multipartin
#
#multipartout, envelopeout = lattice.evaluate(multipartout, envelopeout)
#
##plotEverything(multipartin, twissin, multipartout)

# references
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition
# A. 7.2. Space Charge Impulses in simulatingbeamswithellipsoidalsymmetry-secondedition
# B. A MODIFIED QUADSCAN TECHNIQUE FOR EMITTANCE.pdf
# C. Accelerator-Recipies.pdf by E. Laface
# D. The leapfrog method and other symplectic algorithms for integrating Newtons laws of motion Peter Young Dated April 21 2014
# E. ESS Linac simulator
# F. WEPEA 040