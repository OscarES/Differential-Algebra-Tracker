import numpy as np
from accelerator import Lattice, Element, LinearElement, Quad, Drift


sigma = 0.001 # the standard deviation that the user will enter
#epsilon = sqrt(sigmax**2*sigmaxp**2-sigmaxxp**2)




## Non-linear
order = 5 # user sets this

## Multiparticle
x = np.array([0.0])
xp = np.array([1.0])
y = np.array([0.0])
yp = np.array([1.0])
z = np.array([0.0])
zp = np.array([0.0])
s = 1.0 # from ref 1 in section 2.3
zvector = np.array([x, xp, y, yp, z, zp])
particle1 = np.array([zvector, s])
particle2 = np.array([-zvector, s])
multipart = np.array([particle1, particle2])
#print "len(np.atleast_1d(multipart))" + str(len(np.atleast_1d(multipart)))
#print "multipart[0][0:6]" + str(multipart[0][0:6])
#print "multipart[1] (s)" + str(multipart[1])
print "multipart" + str(multipart)
#print "len(multipart) " + str(len(multipart))
# envelope comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
envelope = np.array([0.0, 10.3338028723, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06, -3.331460652e-16, 8.85901414121, 1e-06])

## Lattice construction
spaceChargeOn = 1
drift = Drift('drift', 1, spaceChargeOn, multipart, envelope)

quad = Quad('quad', 0.1, 1, spaceChargeOn, multipart, envelope)

lattice = Lattice('lattice')
lattice.appendElement(drift)
lattice.appendElement(quad)
partres, envres = lattice.evaluate(multipart,envelope)

print "partres: " + str(partres)
#print "partres[0]: " + str(partres[0]) # don't have to worry about the dtype=object
#print "partres[0][0]: " + str(partres[0][0]) # don't have to worry about the dtype=object
## check if isLinear() works
#print drift.isLinear()
#if 1:
#    print "1 is yes"
#
#if not 0:
#    print "0 is no"

## References
# 1. simulatingbeamswithellipsoidalsymmetry-secondedition