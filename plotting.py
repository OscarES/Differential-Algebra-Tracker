from __future__ import division # needed for 1/2 = 0.5
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from scipy import constants
from particleFactory import envelopeFromMultipart

def plotEverything(multipartin,twiss,multipartout, envlist, multipartafterall):#,envx,envy):
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
    envz = [envlist[i][0][6] for i in xrange(len(envlist))]
    s = [envlist[i][1] for i in xrange(len(envlist))]

    ## begin ellipse, see appendix C in Wille
    # initial
    gamma_x = (1+alpha_x**2)/beta_x
    
    #sigma_xp = np.sqrt(epsilon_rms_x*gamma_x)# b, from fig 3.23
    #sigma_x = np.sqrt(epsilon_rms_x*beta_x) # a
    sigma_xp = np.sqrt(envlist[0][0][2]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_x = np.sqrt(envx[0])
    three_sigma_x = 3*sigma_x
    three_sigma_xp = 3*sigma_xp
    angle_x = angleFrom_i_and_ip(xin, xpin) # in radians
    ellipse_x = Ellipse((0,0), 2*three_sigma_x, 2*three_sigma_xp, angle_x*180/constants.pi) # 2* is because matplotlibs height is from bottom to top and not center to top...
    ellipse_x.set_facecolor('none')
    ellipse_x.set_edgecolor((0,0,0))


    gamma_y = (1+alpha_y**2)/beta_y
    #sigma_yp = np.sqrt(epsilon_rms_y*gamma_y)# b, from fig 3.23
    #sigma_y = np.sqrt(epsilon_rms_y*beta_y) # a
    sigma_yp = np.sqrt(envlist[0][0][5]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_y = np.sqrt(envy[0])
    three_sigma_y = 3*sigma_y
    three_sigma_yp = 3*sigma_yp
    angle_y = angleFrom_i_and_ip(yin, ypin)
    ellipse_y = Ellipse((0,0), 2*three_sigma_y, 2*three_sigma_yp, angle_y*180/constants.pi) # 2* is because matplotlibs height is from bottom to top and not center to top...
    ellipse_y.set_facecolor('none')
    ellipse_y.set_edgecolor((0,0,0))

    ellipse_xy = Ellipse((0,0), 2*three_sigma_x, 2*three_sigma_y, 0)
    ellipse_xy.set_facecolor('none')
    ellipse_xy.set_edgecolor((0,0,0))

    gamma_z = (1+alpha_z**2)/beta_z
    sigma_zp = np.sqrt(envlist[0][0][8]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_z = np.sqrt(envz[0])
    ## end initial ellipses

    xo = [multipartout[i][0][0] for i in xrange(len(multipartout))]
    xpo = [multipartout[i][0][1] for i in xrange(len(multipartout))]
    yo = [multipartout[i][0][2] for i in xrange(len(multipartout))]
    ypo = [multipartout[i][0][3] for i in xrange(len(multipartout))]
    zo = [multipartout[i][0][4] for i in xrange(len(multipartout))]
    zpo = [multipartout[i][0][5] for i in xrange(len(multipartout))]

    xallo = [multipartafterall[i][0][0] for i in xrange(len(multipartafterall))]
    xpallo = [multipartafterall[i][0][1] for i in xrange(len(multipartafterall))]
    #print "xallo: \n" + str(xallo)
    #print "xpallo: \n" + str(xpallo)

    #print "xo: \n" + str(xo)
    #print "xpo: \n" + str(xpo)
    # Cull inf values
    #if len(xo) > 1:
    #    xo[xo > 1e6] = 0
    #    xpo[xpo > 1e6] = 0
    #    xallo[xo > 1e6] = 0
    #    xpallo[xpo > 1e6] = 0

    # after ellipses
    envelopeAfter = envelopeFromMultipart(multipartout)

    #angle_x_after = angleFrom_i_and_ip(xo, xpo)

    sigma_xTimessigma_xp_after = envelopeAfter[1]

    sigma_xp_after = np.sqrt(envelopeAfter[2]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_x_after = np.sqrt(envelopeAfter[0])
    three_sigma_x_after = 3*sigma_x_after
    three_sigma_xp_after = 3*sigma_xp_after

    #ellipse_x_after = Ellipse((0,0), 2*three_sigma_x_after, 2*three_sigma_xp_after, angle_x_after*180/constants.pi) # 2* is because matplotlibs height is from bottom to top and not center to top....
    #ellipse_x_after.set_facecolor('none')
    #ellipse_x_after.set_edgecolor((0,0,0))

    # y
    #angle_y_after = angleFrom_i_and_ip(yo, ypo)

    sigma_yTimessigma_yp_after = envelopeAfter[4]

    sigma_yp_after = np.sqrt(envelopeAfter[5]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_y_after = np.sqrt(envelopeAfter[3])
    three_sigma_y_after = 3*sigma_y_after
    three_sigma_yp_after = 3*sigma_yp_after

    #ellipse_y_after = Ellipse((0,0), 2*three_sigma_y_after, 2*three_sigma_yp_after, angle_y_after*180/constants.pi) # 2* is because matplotlibs height is from bottom to top and not center to top....
    #ellipse_y_after.set_facecolor('none')
    #ellipse_y_after.set_edgecolor((0,0,0))

    # z
    #angle_z_after = angleFrom_i_and_ip(zo, zpo)

    sigma_zTimessigma_zp_after = envelopeAfter[7]

    sigma_zp_after = np.sqrt(envelopeAfter[8]) # sqrt since env has sigma**2. 3*sigma will cover 99.7%. 
    sigma_z_after = np.sqrt(envelopeAfter[6])
    three_sigma_z_after = 3*sigma_z_after
    three_sigma_zp_after = 3*sigma_zp_after

    #ellipse_z_after = Ellipse((0,0), 2*three_sigma_z_after, 2*three_sigma_zp_after, angle_z_after*180/constants.pi) # 2* is because matplotlibs height is from bottom to top and not center to top....
    #ellipse_z_after.set_facecolor('none')
    #ellipse_z_after.set_edgecolor((0,0,0))

    # xy
    ellipse_xy_after = Ellipse((0,0), 2*three_sigma_x_after, 2*three_sigma_y_after, 0)
    ellipse_xy_after.set_facecolor('none')
    ellipse_xy_after.set_edgecolor((0,0,0))
    # end after ellipses

    plt.figure(0)

    ax1 = plt.subplot2grid((4,4), (0,0))
    #ax1.add_artist(ellipse_x)
    ax1.plot(xin,xpin,'ro', zorder=1)
    if not sigma_x == 0:
        ax1.set_xlim(-5*sigma_x, 5*sigma_x)
    if not sigma_xp == 0:
        ax1.set_ylim(-5*sigma_xp, 5*sigma_xp)
    plt.title('Initial values in x')
    plt.xlabel('x [m]')
    plt.ylabel('x\' []')

    ax2 = plt.subplot2grid((4,4), (0,1))
    #ax2.add_artist(ellipse_y)
    ax2.plot(yin,ypin,'bo', zorder=1)
    if not sigma_y == 0:
        ax2.set_xlim(-5*sigma_y, 5*sigma_y)
    if not sigma_yp == 0:
        ax2.set_ylim(-5*sigma_yp, 5*sigma_yp)
    plt.title('Initial values in y')
    plt.xlabel('y [m]')
    plt.ylabel('y\' []')

    ax3 = plt.subplot2grid((4,4), (0, 2))
    #ax8.add_artist(ellipse_y_after)
    ax3.plot(zin,zpin,'yo', zorder=1)
    if not sigma_z == 0:
        ax3.set_xlim(-5*sigma_z, 5*sigma_z)
    if not sigma_zp == 0:
        ax3.set_ylim(-5*sigma_zp, 5*sigma_zp)
    #ax3.set_xlim(-4e-3, 4e-3)
    #ax3.set_ylim(-4e-3, 4e-3)
    plt.title('Values before lattice in z')
    plt.xlabel('z [m]')
    plt.ylabel('$\delta$ []')

    ax4 = plt.subplot2grid((4,4), (0,3))
    ax4.add_artist(ellipse_xy)
    ax4.plot(xin,yin,'go', zorder=1)
    if not sigma_x == 0:
        ax4.set_xlim(-5*sigma_x, 5*sigma_x)
    if not sigma_y == 0:
        ax4.set_ylim(-5*sigma_y, 5*sigma_y)
    plt.title('Initial values x and y')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    ax5 = plt.subplot2grid((4,4), (1, 0), colspan=3)
    ax5.plot(s,envx,'ro')
    plt.title('Envelope in sigma_x^2 by s')
    plt.xlabel('s [m]')
    plt.ylabel('Envelope in sigma_x^2')

    ax6 = plt.subplot2grid((4,4), (2, 0), colspan=3)
    ax6.plot(s,envy,'bo')
    plt.title('Envelope in sigma_y^2 by s')
    plt.xlabel('s [m]')
    plt.ylabel('Envelope in sigma_y^2')

    ax7 = plt.subplot2grid((4,4), (2, 3))
    #ax7.add_artist(ellipse_x_after)
    ax7.plot(xallo,xpallo,'r.', zorder=1)
    if not sigma_x_after == 0:
        ax7.set_xlim(-5*sigma_x_after, 5*sigma_x_after)
    if not sigma_xp_after == 0:
        ax7.set_ylim(-5*sigma_xp_after, 5*sigma_xp_after)
    #ax7.set_xlim(-0.004, 0.004)
    #ax7.set_ylim(-0.004, 0.004)
    plt.title('Values after all elems in lattice in x')
    plt.xlabel('x [m]')
    plt.ylabel('x\' []')

    ax8 = plt.subplot2grid((4,4), (3, 0))
    #ax8.add_artist(ellipse_x_after)
    ax8.plot(xo,xpo,'ro', zorder=1)
    if not sigma_x_after == 0:
        ax8.set_xlim(-5*sigma_x_after, 5*sigma_x_after)
    if not sigma_xp_after == 0:
        ax8.set_ylim(-5*sigma_xp_after, 5*sigma_xp_after)
    #ax8.set_xlim(-0.004, 0.004)
    #ax8.set_ylim(-0.004, 0.004)
    plt.title('Values after lattice in x')
    plt.xlabel('x [m]')
    plt.ylabel('x\' []')

    ax9 = plt.subplot2grid((4,4), (3, 1))
    #ax9.add_artist(ellipse_y_after)
    ax9.plot(yo,ypo,'bo', zorder=1)
    if not sigma_y_after == 0:
        ax9.set_xlim(-5*sigma_y_after, 5*sigma_y_after)
    if not sigma_yp_after == 0:
        ax9.set_ylim(-5*sigma_yp_after, 5*sigma_yp_after)
    plt.title('Values after lattice in y')
    plt.xlabel('y [m]')
    plt.ylabel('y\' []')

    ax10 = plt.subplot2grid((4,4), (3, 2))
    #ax10.add_artist(ellipse_y_after)
    ax10.plot(zo,zpo,'yo', zorder=1)
    if not sigma_z_after == 0:
        ax10.set_xlim(-5*sigma_z_after, 5*sigma_z_after)
    if not sigma_zp_after == 0:
        ax10.set_ylim(-5*sigma_zp_after, 5*sigma_zp_after)
    plt.title('Values after lattice in z')
    plt.xlabel('z [m]')
    plt.ylabel('$\delta$ []')

    ax11 = plt.subplot2grid((4,4), (3, 3))
    ax11.add_artist(ellipse_xy_after)
    ax11.plot(xo,yo,'go', zorder=1)
    if not sigma_x_after == 0:
        ax11.set_xlim(-5*sigma_x_after, 5*sigma_x_after)
    if not sigma_y_after == 0:
        ax11.set_ylim(-5*sigma_y_after, 5*sigma_y_after)
    plt.title('Values after lattice in y and x')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    plt.suptitle("Plots")

    left  = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.2   # the amount of width reserved for blank space between subplots
    hspace = 0.4   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left, bottom, right, top, wspace, hspace) # Increases the space between subplots


    plt.show()

def plot_x_before_and_after(multipartin, multipartout):

    xin = [multipartin[i][0][0] for i in xrange(len(multipartin))]
    xpin = [multipartin[i][0][1] for i in xrange(len(multipartin))]

    xo = [multipartout[i][0][0] for i in xrange(len(multipartout))]
    xpo = [multipartout[i][0][1] for i in xrange(len(multipartout))]

    matplotlib.rcParams.update({'font.size': 31})
    plt.figure(0)

    ax1 = plt.subplot2grid((1,2), (0,0))
    ax1.plot(xin,xpin,'ro', zorder=1)
    plt.title('Initial values in x')
    plt.xlabel('x [m]')
    plt.ylabel('x\' []')

    ax2 = plt.subplot2grid((1,2), (0,1))
    ax2.plot(xo,xpo,'ro', zorder=1)
    plt.title('Values after lattice in x')
    plt.xlabel('x [m]')
    plt.ylabel('x\' []')

    plt.show()

def plot_z_before_and_after(multipartin, multipartout):

    zin = [multipartin[i][0][4] for i in xrange(len(multipartin))]
    zpin = [multipartin[i][0][5] for i in xrange(len(multipartin))]

    zo = [multipartout[i][0][4] for i in xrange(len(multipartout))]
    zpo = [multipartout[i][0][5] for i in xrange(len(multipartout))]

    matplotlib.rcParams.update({'font.size': 31})
    plt.figure(0)

    ax1 = plt.subplot2grid((1,2), (0,0))
    ax1.plot(zin,zpin,'ro', zorder=1)
    plt.title('Initial values in z')
    plt.xlabel('z [m]')
    plt.ylabel('$\delta$ []')

    ax2 = plt.subplot2grid((1,2), (0,1))
    ax2.plot(zo,zpo,'ro', zorder=1)
    plt.title('Values after lattice in z')
    plt.xlabel('z [m]')
    plt.ylabel('$\delta$ []')

    plt.show()

def angleFrom_i_and_ip(i,ip):
    if len(i) == 1:
        return 0
    data = np.array([i, ip])
    covM = np.cov(data)

    eigenvalues, evectors = np.linalg.eig(covM)
    largest_eigenvalue_index = np.argmax(eigenvalues)

    eigenvector_with_largest_eigenvalue = evectors[:,largest_eigenvalue_index]

    angleFromCov = np.arctan(eigenvector_with_largest_eigenvalue[1]/eigenvector_with_largest_eigenvalue[0])
    return angleFromCov

def plotEnvelope(envx,envy):
    plt.figure(0)
    ax4 = plt.subplot2grid((2,3), (0, 0), colspan=3)
    plt.plot(envx[:,0],envx[:,1],'ro')
    plt.title('Envelope in x by z')
    plt.xlabel('z [m]')
    plt.ylabel('Envelope in x [m]')
    ax5 = plt.subplot2grid((2,3), (1, 0), colspan=3)
    plt.plot(envy[:,0],envy[:,1],'bo')
    plt.title('Envelope in y by z')
    plt.xlabel('z [m]')
    plt.ylabel('Envelope in y [m]')

    plt.show()

## Plot def
def plotPhaseSpace(x,xp,y,yp):
    plt.subplot(121)
    plt.plot(x,xp,'ro')
    plt.xlabel('x [m]')
    plt.ylabel('xp [rad]')
    
    plt.subplot(122)
    plt.plot(y,yp,'ro')
    plt.xlabel('y [m]')
    plt.ylabel('yp [rad]')

    plt.show()