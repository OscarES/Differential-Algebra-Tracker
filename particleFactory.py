import numpy as np
import random

def straight(nbrOfParticles):
    x = np.linspace(-1,1,nbrOfParticles) #x and xp for when xp is 0
    xp = np.linspace(0,0,nbrOfParticles)
    y = np.linspace(-1,1,nbrOfParticles) #y and yp for when yp is 0
    yp = np.linspace(0,0,nbrOfParticles)

    z = np.linspace(0,0,nbrOfParticles)
    zp = np.linspace(0,0,nbrOfParticles)

    s = np.linspace(0,0,nbrOfParticles)

    bigMatrix = np.array([x, xp, y, yp, z, zp])

    multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

    return multipart

    #return x,xp,y,yp,z,zp

def scanned(nbrOfParticles):
    x = np.linspace(-0.05,0.05,nbrOfParticles)     # x and xp for when xp is scanned
    xp = np.linspace(-0.0001,0.0001,nbrOfParticles)
    y = np.linspace(-0.05,0.05,nbrOfParticles)     # y and yp for when yp is scanned
    yp = np.linspace(-0.0001,0.0001,nbrOfParticles)

    z = np.linspace(0,0,nbrOfParticles)
    zp = np.linspace(0,0,nbrOfParticles)

    s = np.linspace(0,0,nbrOfParticles)

    bigMatrix = np.array([x, xp, y, yp, z, zp])

    multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

    return multipart

    #return x,xp,y,yp,z,zp
                           
def randomed(nbrOfParticles):
    x = [random.uniform(-0.05, 0.05) for _ in xrange(nbrOfParticles)]
    xp = [random.uniform(-0.00001, 0.00001) for _ in xrange(nbrOfParticles)]
    y = [random.uniform(-0.05, 0.05) for _ in xrange(nbrOfParticles)]
    yp = [random.uniform(-0.00001, 0.00001) for _ in xrange(nbrOfParticles)]

    z = np.linspace(0,0,nbrOfParticles)
    zp = np.linspace(0,0,nbrOfParticles)

    s = np.linspace(0,0,nbrOfParticles)

    bigMatrix = np.array([x, xp, y, yp, z, zp])

    multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

    return multipart

    #return x,xp,y,yp,z,zp

def gaussian(nbrOfParticles):
    x = np.random.normal(0,0.001,nbrOfParticles)
    xp = np.random.normal(0,0.000001,nbrOfParticles)
    y = np.random.normal(0,0.001,nbrOfParticles)
    yp = np.random.normal(0,0.000001,nbrOfParticles)

    z = np.linspace(0,0,nbrOfParticles)
    zp = np.linspace(0,0,nbrOfParticles)

    s = np.linspace(0,0,nbrOfParticles)

    bigMatrix = np.array([x, xp, y, yp, z, zp])

    multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

    return multipart

    #return x,xp,y,yp,z,zp

def gaussiantwiss(nbrOfParticles, alpha, beta, epsilon):
    xi = np.random.normal(0,1,nbrOfParticles)
    xip = np.random.normal(0,1,nbrOfParticles)

    M = np.array([
        [1/sqrt(beta*epsilon), 0],
        [alpha/sqrt(beta*epsilon), sqrt(beta/epsilon)]
        ])

    Minv = np.linalg.inv(M)

    x = np.zeros(nbrOfParticles)
    xp = np.zeros(nbrOfParticles)

    for i in range(nbrOfParticles):
        x[i] = Minv[0,0]*xi[i] + Minv[0,1]*xip[i]
        xp[i] = Minv[1,0]*xi[i] + Minv[1,1]*xip[i]

    return x,xp