import numpy as np
from accelerator import Lattice, Element, LinearElement, Quad, Drift, DifferentialAlgebra, DiffAlgElement
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
from particleFactory import straight, scanned, randomed, gaussian


sigma = 0.001 # the standard deviation that the user will enter
#epsilon = sqrt(sigmax**2*sigmaxp**2-sigmaxxp**2)

## Symbols
qx = Symbol('qx')
qy = Symbol('qy')
qz = Symbol('qz')
px = Symbol('px')
py = Symbol('py')
pz = Symbol('pz')
l = Symbol('l') # arbitrary length
k = Symbol('k')

## Hamiltonians
driftham = -l/2*(px**2 + py**2 + pz**2)
quadham = -l/2*(k**2*(qx**2-qy**2)+px**2+py**2+pz**2) # replace k with -k for defocus. Without quad term in z dir
quadhamdefocus = -l/2*(-k**2*(qx**2-qy**2)+px**2+py**2+pz**2) # replace k with -k for defocus. Without quad term in z dir
sextupoleham = -l/2*(2/3*k**2*(qx**3-3*qx*qy**2)+(px**2+py**2)) # should the ps' perhaps be divided by 2 as in nonlinear2013_3.pdf? That division is assumed to be the l/2 in the beginning, . k is actually k**2
octupoleham = -l/2*(2/4*k*(qx**4-6*qx**2*qy**2+qy**4)+(px**2+py**2)) # same decision as above

## Non-linear
order = 5 # user sets this

DA = DifferentialAlgebra()

quaddefocuskval = 0.1
quaddefocuslval = 1
quaddefocusNumFuns = DA.hamToNumFuns(quadhamdefocus, quaddefocuskval, quaddefocuslval, order)

#print "quaddefocusNumFuns:" + str(quaddefocusNumFuns)

#diffAlgRes = (quaddefocusNumFuns[0](particle1), quaddefocusNumFuns[1](particle1), quaddefocusNumFuns[2](particle1), quaddefocusNumFuns[3](particle1), quaddefocusNumFuns[4](particle1), quaddefocusNumFuns[5](particle1))

#sextukval = 
#sextupoleNumFuns = DA.hamToNumFuns(sextupoleham, kval, lval, order)

## Multiparticle
x = 0.0
xp =1.0
y = 0.0
yp =1.0
z = 0.0
zp =0.0
s = 1.0 # from ref 1 in section 2.3
zvector = np.array([x, xp, y, yp, z, zp])
particle1 = np.array([zvector, s])
particle2 = np.array([-zvector, s])
multipart = np.array([particle1, particle2])
multipart2 = np.array([particle1, particle2])
#print "len(np.atleast_1d(multipart))" + str(len(np.atleast_1d(multipart)))
#print "multipart[0][0:6]" + str(multipart[0][0:6])
#print "multipart[1] (s)" + str(multipart[1])
print "multipart" + str(multipart)
#print "len(multipart) " + str(len(multipart))
# envelope comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
envelope = np.array([0.0, 10.3338028723, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06])

## Lattice construction
spaceChargeOn = 0
drift = Drift('drift', 1, spaceChargeOn, multipart, envelope)

quad = Quad('quad', 0.1, 1, spaceChargeOn, multipart, envelope)

lattice = Lattice('lattice')
#lattice.appendElement(drift)
lattice.appendElement(quad)
#print "multipart before latt eval: " + str(multipart)
partres, envres = lattice.evaluate(multipart,envelope) ### changes multipart and envelope!!!!!!!!!!!!!!!!!!
#print "multipart2 after latt eval: " + str(multipart2)

print "partres: " + str(partres)

# Using the differential algebra "raw"
diffAlgRes = (quaddefocusNumFuns[0](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[1](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[2](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[3](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[4](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]), quaddefocusNumFuns[5](particle1[0][0],particle1[0][1],particle1[0][2],particle1[0][3],particle1[0][4],particle1[0][5]))
print "DiffAlgRes: " + str(diffAlgRes) # MATCHES the linear calculation!!!!! :)

# Using the differential algebra through DiffAlgElement
diffAlgElemQuad = DiffAlgElement("diffAlgElemQuad", DA, quadhamdefocus, quaddefocuskval, quaddefocuslval, order, 0, multipart2, envelope)
diffAlgElemQuadpartres, diffAlgElemQuadenvres = diffAlgElemQuad.evaluate(multipart2,envelope) # not the same as raw! Since the element is split! But if n = 1 in the class then there is a very good match!!!! :)
print "diffAlgElemQuadRes: " + str(diffAlgElemQuadpartres)


### Tons of particles
print "Tons of particles..."
nbrOfParticles = 10
multiparttonsMat = straight(nbrOfParticles)
multiparttonsDiffAlg = straight(nbrOfParticles)

#print "multiparttonsMat: \n" + str(multiparttonsMat)

tonsMatPartRes, tonsMatEnvRes = lattice.evaluate(multiparttonsMat,envelope)
tonsDiffAlgPartRes, tonsDiffAlgEnvRes = diffAlgElemQuad.evaluate(multiparttonsDiffAlg,envelope)

#print "tonsMatPartRes: \n" + str(tonsMatPartRes)
#print "tonsDiffAlgPartRes: \n" + str(tonsDiffAlgPartRes)

# Checking the others
scannedparts = scanned(nbrOfParticles)
#print "scannedparts: \n" + str(scannedparts)

randomedparts = randomed(nbrOfParticles)
#print "randomedparts: \n" + str(randomedparts)

gaussianparts = gaussian(nbrOfParticles)
#print "gaussianparts: \n" + str(gaussianparts)



## References
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition