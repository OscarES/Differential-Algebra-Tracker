from __future__ import division # needed for 1/2 = 0.5
import math
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

def gaussianTwiss3D(nbrOfParticles, twiss):
    alpha_x = twiss[0]
    beta_x = twiss[1]
    epsilon_rms_x = twiss[2]
    alpha_y = twiss[3]
    beta_y = twiss[4]
    epsilon_rms_y = twiss[5]
    alpha_z = twiss[6]
    beta_z = twiss[7]
    epsilon_rms_z = twiss[8]

    x, xp = gaussianTwiss1D(nbrOfParticles, alpha_x, beta_x, epsilon_rms_x)
    y, yp = gaussianTwiss1D(nbrOfParticles, alpha_y, beta_y, epsilon_rms_y)
    z, zp = gaussianTwiss1D(nbrOfParticles, alpha_z, beta_z, epsilon_rms_z)

    s = np.linspace(0,0,nbrOfParticles)

    bigMatrix = np.array([x, xp, y, yp, z, zp])

    multipart = [[bigMatrix[:,i], s[i]] for i in xrange(nbrOfParticles)]

    return multipart

def gaussianTwiss1D(nbrOfParticles, alpha, beta, epsilon):
    xi = np.random.normal(0,1,nbrOfParticles)
    xip = np.random.normal(0,1,nbrOfParticles)

    M = np.array([
        [1/np.sqrt(beta*epsilon), 0],
        [alpha/np.sqrt(beta*epsilon), np.sqrt(beta/epsilon)]
        ])

    Minv = np.linalg.inv(M)

    x = np.zeros(nbrOfParticles)
    xp = np.zeros(nbrOfParticles)

    for i in range(nbrOfParticles):
        x[i] = Minv[0,0]*xi[i] + Minv[0,1]*xip[i]
        xp[i] = Minv[1,0]*xi[i] + Minv[1,1]*xip[i]

    return x,xp

def envelopeFromMultipart(multipart):
    x = [multipart[i][0][0] for i in xrange(len(multipart))]
    xp = [multipart[i][0][1] for i in xrange(len(multipart))]
    y = [multipart[i][0][2] for i in xrange(len(multipart))]
    yp = [multipart[i][0][3] for i in xrange(len(multipart))]
    z = [multipart[i][0][4] for i in xrange(len(multipart))]
    zp = [multipart[i][0][5] for i in xrange(len(multipart))]

    #sigma_x = np.mean(x) # sigma isn't std! It is mean instead!
    #sigma_xp = np.mean(xp)
    #sigma_y = np.mean(y)
    #sigma_yp = np.mean(yp)
    #sigma_z = np.mean(z)
    #sigma_zp = np.mean(zp)

    #sigma_x = np.mean(map(abs, x)) # sigma isn't std! It is abs mean perhaps!
    #sigma_xp = np.mean(map(abs, xp))
    #sigma_y = np.mean(map(abs, y))
    #sigma_yp = np.mean(map(abs, yp))
    #sigma_z = np.mean(map(abs, z))
    #sigma_zp = np.mean(map(abs, zp))



    # Convert to float so that std won't complain
    x = np.array(x)
    xp = np.array(xp)
    y = np.array(y)
    yp = np.array(yp)
    z = np.array(z)
    zp = np.array(zp)

    x = x.astype('float')
    xp = xp.astype('float')
    y = y.astype('float')
    yp = yp.astype('float')
    z = z.astype('float')
    zp = zp.astype('float')

    sigma_x = np.std(x) # sigma is std
    sigma_xp = np.std(xp)
    sigma_y = np.std(y)
    sigma_yp = np.std(yp)
    sigma_z = np.std(z)
    sigma_zp = np.std(zp)

    envelope = [sigma_x**2, sigma_x*sigma_xp, sigma_xp**2, sigma_y**2, sigma_y*sigma_yp, sigma_yp**2, sigma_z**2, sigma_z*sigma_zp, sigma_zp**2]
    return envelope