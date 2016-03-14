import matplotlib.pyplot as plt

def plotEverything(multipartin,twiss,multipartout, envlist):#,envx,envy):
    xin = [multipartin[i][0][0] for i in xrange(len(multipartin))]
    xpin = [multipartin[i][0][1] for i in xrange(len(multipartin))]
    yin = [multipartin[i][0][2] for i in xrange(len(multipartin))]
    ypin = [multipartin[i][0][3] for i in xrange(len(multipartin))]
    zin = [multipartin[i][0][4] for i in xrange(len(multipartin))]
    zpin = [multipartin[i][0][5] for i in xrange(len(multipartin))]

    alpha_x = twiss[0]
    beta_x = twiss[1]
    epsilon_rms_x = twiss[2]
    alpha_y = twiss[3]
    beta_y = twiss[4]
    epsilon_rms_y = twiss[5]
    alpha_z = twiss[6]
    beta_z = twiss[7]
    epsilon_rms_z = twiss[8]

    envx = [envlist[i][0][0] for i in xrange(len(envlist))]
    envy = [envlist[i][0][3] for i in xrange(len(envlist))]
    s = [envlist[i][1] for i in xrange(len(envlist))]

    xo = [multipartout[i][0][0] for i in xrange(len(multipartout))]
    xpo = [multipartout[i][0][1] for i in xrange(len(multipartout))]
    yo = [multipartout[i][0][2] for i in xrange(len(multipartout))]
    ypo = [multipartout[i][0][3] for i in xrange(len(multipartout))]
    zo = [multipartout[i][0][4] for i in xrange(len(multipartout))]
    zpo = [multipartout[i][0][5] for i in xrange(len(multipartout))]

    plt.figure(0)
    ax1 = plt.subplot2grid((4,3), (0,0))
    plt.plot(xin,xpin,'ro')
    plt.title('Initial value in x')
    plt.xlabel('x')
    plt.ylabel('xp')
    ax2 = plt.subplot2grid((4,3), (0,1))
    plt.plot(yin,ypin,'bo')
    plt.title('Initial value in y')
    plt.xlabel('y')
    plt.ylabel('yp')
    ax3 = plt.subplot2grid((4,3), (0,2))
    plt.plot(xin,yin,'go')
    plt.title('Initial values x and y')
    plt.xlabel('x')
    plt.ylabel('y')
    ax4 = plt.subplot2grid((4,3), (1, 0), colspan=3)
    plt.plot(s,envx,'ro')
    plt.title('Envelope in sigma_x^2 by s')
    plt.xlabel('s')
    plt.ylabel('Envelope in sigma_x^2')
    ax5 = plt.subplot2grid((4,3), (2, 0), colspan=3)
    plt.plot(s,envy,'bo')
    plt.title('Envelope in sigma_y^2 by s')
    plt.xlabel('s')
    plt.ylabel('Envelope in sigma_y^2')
    ax6 = plt.subplot2grid((4,3), (3, 0))
    plt.plot(xo,xpo,'ro')
    plt.title('Values after FODO in x')
    plt.xlabel('x')
    plt.ylabel('xp')
    ax7 = plt.subplot2grid((4,3), (3, 1))
    plt.plot(yo,ypo,'bo')
    plt.title('Values after FODO in y')
    plt.xlabel('y')
    plt.ylabel('yp')
    ax8 = plt.subplot2grid((4,3), (3, 2))
    plt.plot(xo,yo,'go')
    plt.title('Values after FODO in y and x')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.suptitle("Plots")
    plt.show()

def plotEnvelope(envx,envy):
    plt.figure(0)
    ax4 = plt.subplot2grid((2,3), (0, 0), colspan=3)
    plt.plot(envx[:,0],envx[:,1],'ro')
    plt.title('Envelope in x by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in x')
    ax5 = plt.subplot2grid((2,3), (1, 0), colspan=3)
    plt.plot(envy[:,0],envy[:,1],'bo')
    plt.title('Envelope in y by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in y')

    plt.show()

## Plot def
def plotPhaseSpace(x,xp,y,yp):
    plt.subplot(121)
    plt.plot(x,xp,'ro')
    plt.xlabel('x')
    plt.ylabel('xp')
    
    plt.subplot(122)
    plt.plot(y,yp,'ro')
    plt.xlabel('y')
    plt.ylabel('yp')

    plt.show()