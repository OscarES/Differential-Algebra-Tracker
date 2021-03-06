import math

import sys
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import GLU
from PyQt5 import QtGui
from PyQt5.QtOpenGL import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5 import QtCore
from PyQt5 import QtOpenGL
from PyQt5.QtWidgets import QMainWindow, QAction, qApp, QApplication, QInputDialog
from IOHandler import saveLattice, loadLattice, loadSummer2015Formatzasx, saveBeamdata, loadBeamdata, saveTwiss, loadTwiss, saveMultipart, loadMultipart, saveEnvelope # loadLattice, saveLattice
import numpy as np
from scipy import *
from facility import *
from relativity import *
import re

class LatticeOverviewWidget(QGLWidget):
    '''
    Widget for giving an overview of the lattice
    '''
    def __init__(self, facility, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.setMinimumSize(600, 500)

        self.z = 2.0

        self.w_pressed = 0
        self.a_pressed = 0
        self.s_pressed = 0
        self.d_pressed = 0

        self.cameraPos = np.array([0.0, 0.0, 1.0])
        self.cameraTarget = np.array([0.0, 0.0, 0.0])
        diff = self.cameraPos - self.cameraTarget
        self.cameraDirection = diff/np.linalg.norm(diff)
        self.cameraSpeed = 0.5
        self.cameraSpeed_s = 0.5

        self.up = np.array([0.0, 1.0, 0.0])
        tempCross = np.cross(self.up, self.cameraDirection)
        self.cameraRight = tempCross/np.linalg.norm(tempCross)
        self.cameraUp = np.cross(self.cameraDirection, self.cameraRight)

        # extract from lattice
        self.facility = facility

    def loadLattice(self):
        lattticeString = self.facility.printLattice()
        self.elements = []
        nextWillBeL = 0
        print lattticeString
        for line in lattticeString.split():
            if nextWillBeL:
                self.elements.append([tempword, float(line)])
                self.cameraSpeed_s = float(line)
                nextWillBeL = 0
            if line == "drift" or line == "dipole" or line == "quad" or line == "liealgelem" or line == "rotation" or line == "cavity" or line == "cavmat":  
                tempword = line
            if line == "L:":
                nextWillBeL = 1

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0,  150))

        self.loadLattice()

        self.array_of_blocks = []
        for elem in self.elements:
            if elem[0] == "drift":
                color = [0, 0, 0]
            elif elem[0] == "dipole":
                color = [0, 1, 0]
            elif elem[0] == "quad":
                color = [1, 0, 0]
            elif elem[0] == "liealgelem":
                color = [0, 0, 1]
            elif elem[0] == "rotation":
                color = [0, 0, 0]
            elif elem[0] == "cavity" or elem[0] == "cavmat":
                color = [1, 1, 0]

            length = elem[1]

            if elem[0] == "liealgelem":
                VtxArray, IdxArray, ClrArray = self.createHexaBlock(length, color)
            elif elem[0] == "cavity" or elem[0] == "cavmat":
                VtxArray, IdxArray, ClrArray = self.createCylinder(length, color)
            else:
                VtxArray, IdxArray, ClrArray = self.createGeomBlock(length, color)
            block = [VtxArray, IdxArray, ClrArray, length, elem[0]]
            self.array_of_blocks.append(block)

        glEnable(GL_DEPTH_TEST)

    def resizeGL(self, width, height):
        if height == 0: height = 1

        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def lookAt(self, right, up, direction, position):
        a = np.array([
                    [right[0], right[1], right[2], 0.0],
                    [up[0], up[1], up[2], 0.0],
                    [direction[0], direction[1], direction[2], 0.0],
                    [0.0, 0.0, 0.0, 1.0]
                    ])
        b = np.array([
                    [1.0, 0.0, 0.0, -position[0]],
                    [0.0, 1.0, 0.0, -position[1]],
                    [0.0, 0.0, 1.0, -position[2]],
                    [0.0, 0.0, 0.0, 1.0]
                    ])
        return a*b

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        if self.w_pressed:
            self.cameraPos += self.cameraDirection*self.cameraSpeed
            self.w_pressed = 0
        if self.s_pressed:
            self.cameraPos -= self.cameraDirection*self.cameraSpeed
            self.s_pressed = 0
        if self.a_pressed:
            tempCross = np.cross(self.cameraDirection, self.up)
            normTempCross = tempCross/np.linalg.norm(tempCross)
            self.cameraPos -= normTempCross*self.cameraSpeed_s
            self.a_pressed = 0
        if self.d_pressed:
            tempCross = np.cross(self.cameraDirection, self.up)
            normTempCross = tempCross/np.linalg.norm(tempCross)
            self.cameraPos += normTempCross*self.cameraSpeed_s
            self.d_pressed = 0

        self.zsofar = 0
        for block in self.array_of_blocks:
            self.elemPaint(block)

        # draw when there are no elements
        if self.array_of_blocks == []:
            glLoadIdentity()
            glEnableClientState(GL_VERTEX_ARRAY)
            glEnableClientState(GL_COLOR_ARRAY)
            glVertexPointerf([0,0,0])
            glColorPointerf([0,0,0])
            glDrawElementsui(GL_TRIANGLES, [0,0,0])

        
    def elemPaint(self, elem):
        glLoadIdentity()
        glTranslate(self.zsofar, 0.0, -5.0)
        self.zsofar = self.zsofar + elem[3]
        glTranslate(self.cameraPos[0], self.cameraPos[1], self.cameraPos[2])
        glRotate(90, 0.0, 1.0, 0.0)
        elemtype = elem[4]
        if elemtype == "quad":
            glRotate(45, 0.0, 0.0, 1.0)
        elif elemtype == "liealgelem" or elemtype == "cavity" or elemtype == "cavmat":
            glTranslate(0.5, 0.5, 0.0)
        glTranslate(-0.5, -0.5, -0.5)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointerf(elem[0])
        glColorPointerf(elem[2])
        if elemtype == "liealgelem" or elemtype == "cavity" or elemtype == "cavmat":
            glDrawElementsui(GL_TRIANGLES, elem[1])
        else:
            glDrawElementsui(GL_QUADS, elem[1])

    # this will create rectangular blocks
    def createGeomBlock(self, z, color):
        blockVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, z],
                 [1.0, 0.0, z],
                 [1.0, 1.0, z],
                 [0.0, 1.0, z]])
        blockIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        blockClrArray = np.zeros((8,3))
        blockClrArray[:,0] = color[0]
        blockClrArray[:,1] = color[1]
        blockClrArray[:,2] = color[2]
        return blockVtxArray, blockIdxArray, blockClrArray

    # this will create the sextupole blocks
    def createHexaBlock(self, z, color):
        blockVtxArray = array(
                [[0.866025, 0.5, 0.0],      # 0
                 [0.0, 1.0, 0.0],           # 1
                 [-0.866025, 0.5, 0.0],     # 2
                 [-0.866025, -0.5, 0.0],    # 3
                 [0.0,-1.0,0.0],            # 4
                 [0.866025, -0.5,0.0],      # 5
                 [0.0,0.0,0.0],             # 6
                 [0.866025, 0.5, z],        # 7
                 [0.0, 1.0, z],             # 8
                 [-0.866025, 0.5, z],       # 9
                 [-0.866025, -0.5, z],      # 10
                 [0.0,-1.0,z],              # 11
                 [0.866025, -0.5,z],        # 12
                 [0.0,0.0,z]])              # 13
        blockIdxArray = [
                1, 0, 6,
                2, 1, 6,
                3, 2, 6,
                4, 3, 6,
                5, 4, 6,
                0, 5, 6,
                0, 1, 7,
                1, 2, 8,
                2, 3, 9,
                3, 4, 10,
                4, 5, 11,
                5, 0, 12,
                1, 8, 7,
                2, 9, 8,
                3, 10, 9,
                4, 11, 10,
                5, 12, 11,
                0, 7, 12,
                7, 8, 13,
                8, 9, 13,
                9, 10, 13,
                10, 11, 13,
                11, 12, 13,
                12, 7, 13]
        blockClrArray = np.zeros((14,3))
        blockClrArray[:,0] = color[0]
        blockClrArray[:,1] = color[1]
        blockClrArray[:,2] = color[2]
        return blockVtxArray, blockIdxArray, blockClrArray

    # this will create the cavity cylinder
    def createCylinder(self, z, color):
        # n is range(0,16)
        # x from [cos((i+1)/16*2*constants.pi).evalf() for i in n]
        # y from [sin((i+1)/16*2*constants.pi).evalf() for i in n]
        blockVtxArray = array(
                [[0.923880, 0.382683, 0.0],     # 0 ok
                 [0.707107, 0.707107, 0.0],     # 1 ok
                 [0.382683, 0.923880, 0.0],     # 2 ok
                 [0.0, 1.0, 0.0],               # 3 ok
                 [-0.382683, 0.923880, 0.0],    # 4 ok
                 [-0.707107, 0.707107, 0.0],    # 5 ok
                 [-0.923880, 0.382683, 0.0],    # 6 ok
                 [-1.0, 0.0, 0.0],              # 7 ok
                 [-0.923880, -0.382683, 0.0],   # 8 ok
                 [-0.707107, -0.707107, 0.0],   # 9 ok
                 [-0.382683, -0.923880, 0.0],   # 10 ok
                 [0.0,-1.0, 0.0],               # 11 ok
                 [0.382683, -0.923880, 0.0],    # 12 ok
                 [0.707107, -0.707107, 0.0],    # 13 ok
                 [0.923880, -0.382683, 0.0],    # 14 ok
                 [1.0, 0.0, 0.0],               # 15 ok, above x and y ok
                 [0.0, 0.0, 0.0],               # 16 middle back ok
                 [0.923880, 0.382683, z],       # 17 ok
                 [0.707107, 0.707107, z],       # 18 ok
                 [0.382683, 0.923880, z],       # 19 ok
                 [0.0, 1.0, z],                 # 20 ok
                 [-0.382683, 0.923880, z],      # 21 ok
                 [-0.707107, 0.707107, z],      # 22 ok
                 [-0.923880, 0.382683, z],      # 23 ok
                 [-1.0, 0.0, z],                # 24 ok
                 [-0.923880, -0.382683, z],     # 25 ok
                 [-0.707107, -0.707107, z],     # 26 ok
                 [-0.382683, -0.923880, z],     # 27 ok
                 [0.0,-1.0, z],                 # 28 ok
                 [0.382683, -0.923880, z],      # 29 ok
                 [0.707107, -0.707107, z],      # 30 ok
                 [0.923880, -0.382683, z],      # 31 ok
                 [1.0, 0.0, z],                 # 32 ok
                 [0.0, 0.0, z]])                # 33 middle front ok
        blockIdxArray = [
                1, 0, 16, # A begin ok
                2, 1, 16,
                3, 2, 16,
                4, 3, 16,
                5, 4, 16,
                6, 5, 16,
                7, 6, 16,
                8, 7, 16,
                9, 8, 16,
                10, 9, 16,
                11, 10, 16,
                12, 11, 16,
                13, 12, 16,
                14, 13, 16,
                15, 14, 16,
                0, 15, 16, # A end ok
                0, 1, 17, # B begin ok
                1, 2, 18,
                2, 3, 19,
                3, 4, 20,
                4, 5, 21,
                5, 6, 22,
                6, 7, 23,
                7, 8, 24,
                8, 9, 25,
                9, 10, 26,
                10, 11, 27,
                11, 12, 28,
                12, 13, 29,
                13, 14, 30,
                14, 15, 31,
                15, 0, 32, # B end ok
                1, 18, 17, # C begin ok
                2, 19, 18,
                3, 20, 19,
                4, 21, 20,
                5, 22, 21,
                6, 23, 22,
                7, 24, 23,
                8, 25, 24,
                9, 26, 25,
                10, 27, 26,
                11, 28, 27,
                12, 29, 28,
                13, 30, 29,
                14, 31, 30,
                15, 32, 31,
                0, 17, 32, # C end ok
                17, 18, 33, # D begin
                18, 19, 33,
                19, 20, 33,
                20, 21, 33,
                21, 22, 33,
                22, 23, 33,
                23, 24, 33,
                24, 25, 33,
                25, 26, 33,
                26, 27, 33,
                27, 28, 33,
                28, 29, 33,
                29, 30, 33,
                30, 31, 33,
                31, 32, 33,
                32, 17, 33] # D end
        blockClrArray = np.zeros((34,3))
        blockClrArray[:,0] = color[0]
        blockClrArray[:,1] = color[1]
        blockClrArray[:,2] = color[2]
        return blockVtxArray, blockIdxArray, blockClrArray

    # If one clicks on the openGL panel it gets the focus, also this will be called if focus on opengl is required
    def mousePressEvent(self, bla):
        self.setFocus()
        self.updateGL()


class BeamEditor(QGroupBox):
    '''
    Widget for editing the input beam
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)

        self.parent = parent
        self.setTitle("Beam Editor")

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        ## Beamdata
        # beta
        self.textBeta = QLabel("Relativistic beta []:")
        grid.addWidget(self.textBeta, 0, 0)
        
        self.enterBeta = QLineEdit()
        grid.addWidget(self.enterBeta, 0, 1)

        # rf_lambda
        self.textLambda = QLabel("RF-Lambda [m]:")
        grid.addWidget(self.textLambda, 1, 0)
        
        self.enterLambda = QLineEdit()
        grid.addWidget(self.enterLambda, 1, 1)

        # m
        self.textMass = QLabel("Mass of a particle [kg]:")
        grid.addWidget(self.textMass, 2, 0)
        
        self.enterMass = QLineEdit() # should be a ComboBox with m_p, m_e
        grid.addWidget(self.enterMass, 2, 1)

        # q
        self.textCharge = QLabel("Charge of a particle [C]:")
        grid.addWidget(self.textCharge, 3, 0)
        
        self.enterCharge = QLineEdit() # should be a comboBox with e and -e
        grid.addWidget(self.enterCharge, 3, 1)

        # E
        self.textEnergy = QLabel("Energy [J]:")
        grid.addWidget(self.textEnergy, 4, 0)
        
        self.enterEnergy = QLineEdit() # should be a QLabel and show the calculated value from beta
        grid.addWidget(self.enterEnergy, 4, 1)

        # nbrOfParticles
        self.textNbrOfParticles = QLabel("# of particles []:")
        grid.addWidget(self.textNbrOfParticles, 5, 0)
        
        self.enterNbrOfParticles = QLineEdit()
        grid.addWidget(self.enterNbrOfParticles, 5, 1)

        # Current
        self.textCurrent = QLabel("Current [A]:")
        grid.addWidget(self.textCurrent, 6, 0)
        
        self.enterCurrent = QLineEdit()
        grid.addWidget(self.enterCurrent, 6, 1)

        # Set the fields to the default parameters
        self.writeInBeamdataFields()

        # Connect beta and energy fields (must come after beamdata's default parameters has been set)
        self.enterBeta.textChanged.connect(self.setEnergyFromBeta)
        self.enterEnergy.textChanged.connect(self.setBetaFromEnergy)

        ## Twiss
        # twiss comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
        # alpha_x
        self.textAlpha_x = QLabel("Alpha_x []:")
        grid.addWidget(self.textAlpha_x, 0, 3)
        
        self.enterAlpha_x = QLineEdit()
        grid.addWidget(self.enterAlpha_x, 0, 4)

        # beta_x
        self.textBeta_x = QLabel("Beta_x [m]:")
        grid.addWidget(self.textBeta_x, 1, 3)
        
        self.enterBeta_x = QLineEdit()
        grid.addWidget(self.enterBeta_x, 1, 4)
        
        # epsilon_rms_x
        self.textEpsilon_rms_x = QLabel("Epsilon_rms_x [m]:")
        grid.addWidget(self.textEpsilon_rms_x, 2, 3)
        
        self.enterEpsilon_rms_x = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_x, 2, 4)
        
        # alpha_y
        self.textAlpha_y = QLabel("Alpha_y []:")
        grid.addWidget(self.textAlpha_y, 3, 3)
        
        self.enterAlpha_y = QLineEdit()
        grid.addWidget(self.enterAlpha_y, 3, 4)
        
        # beta_y
        self.textBeta_y = QLabel("Beta_y [m]:")
        grid.addWidget(self.textBeta_y, 4, 3)
        
        self.enterBeta_y = QLineEdit()
        grid.addWidget(self.enterBeta_y, 4, 4)
        
        # epsilon_rms_y
        self.textEpsilon_rms_y = QLabel("Epsilon_rms_y [m]:")
        grid.addWidget(self.textEpsilon_rms_y, 5, 3)
        
        self.enterEpsilon_rms_y = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_y, 5, 4)
        
        # alpha_z
        self.textAlpha_z = QLabel("Alpha_z []:")
        grid.addWidget(self.textAlpha_z, 6, 3)
        
        self.enterAlpha_z = QLineEdit()
        grid.addWidget(self.enterAlpha_z, 6, 4)
        
        # beta_z
        self.textBeta_z = QLabel("Beta_z [m]:")
        grid.addWidget(self.textBeta_z, 7, 3)
        
        self.enterBeta_z = QLineEdit()
        grid.addWidget(self.enterBeta_z, 7, 4)
        
        # epsilon_rms_z
        self.textEpsilon_rms_z = QLabel("Epsilon_rms_z [m]:")
        grid.addWidget(self.textEpsilon_rms_z, 8, 3)
        
        self.enterEpsilon_rms_z = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_z, 8, 4)

        # Set the fields to the default parameters
        self.writeInTwissFields()

        # buttons
        useBeamdataButton = QPushButton("Use Beamdata")
        useBeamdataButton.clicked.connect(self.useBeamdataInput)
        grid.addWidget(useBeamdataButton,9,1)

        useTwissButton = QPushButton("Use Twiss")
        useTwissButton.clicked.connect(self.useTwissInput)
        grid.addWidget(useTwissButton,9,4)

        generateGridpartButton = QPushButton("Generate particle grid\n(# of particles squared)")
        generateGridpartButton.clicked.connect(self.generateGridpart)
        grid.addWidget(generateGridpartButton,8,5)

        generateMultipartButton = QPushButton("Generate Multiparticles")
        generateMultipartButton.clicked.connect(self.generateMultipart)
        grid.addWidget(generateMultipartButton,9,5)

        saveBeamdataButton = QPushButton("Save Beamdata")
        saveBeamdataButton.clicked.connect(self.saveBeamdata)
        grid.addWidget(saveBeamdataButton,10,1)

        saveTwissButton = QPushButton("Save Twiss")
        saveTwissButton.clicked.connect(self.saveTwiss)
        grid.addWidget(saveTwissButton,10,4)

        saveMultipartButton = QPushButton("Save Multiparticles")
        saveMultipartButton.clicked.connect(self.saveMultipart)
        grid.addWidget(saveMultipartButton,10,5)

        loadBeamdataButton = QPushButton("Load Beamdata")
        loadBeamdataButton.clicked.connect(self.loadBeamdata)
        grid.addWidget(loadBeamdataButton,11,1)

        loadTwissButton = QPushButton("Load Twiss")
        loadTwissButton.clicked.connect(self.loadTwiss)
        grid.addWidget(loadTwissButton,11,4)

        loadMultipartButton = QPushButton("Load Multiparticles")
        loadMultipartButton.clicked.connect(self.loadMultipart)
        grid.addWidget(loadMultipartButton,11,5)

    def setEnergyFromBeta(self):
        if self.enterEnergy.hasFocus():
            return # Don't change value if enterEnergy field has focus
        self.valueOfBeta = self.enterBeta.text()
        try:
            beta = float(self.valueOfBeta)
        except:
            return # Only run the setBeta when the value is a valid float

        self.valueOfMass = self.enterMass.text()
        try:
            m_0 = float(self.valueOfMass)
        except:
            return # Only run the setBeta when the value is a valid float
        gamma = gammaFromBeta(beta)
        E = EFromBeta(m_0, beta)
        self.enterEnergy.setText(str(E))

    def setBetaFromEnergy(self):
        if self.enterBeta.hasFocus():
            return # Don't change value if enterBeta field has focus
        self.valueOfEnergy = self.enterEnergy.text()
        try:
            E = float(self.valueOfEnergy)
        except:
            return # Only run the setBeta when the value is a valid float

        self.valueOfMass = self.enterMass.text()
        try:
            m_0 = float(self.valueOfMass)
        except:
            return # Only run the setBeta when the value is a valid float
        beta = betaFromE(m_0, E)
        self.enterBeta.setText(str(beta))

    def saveBeamdata(self):
        fname = QFileDialog.getSaveFileName(self, 'Save Beamdata file', '')
        try:
            saveBeamdata(fname[0],self.facility.getBeamdata()) # fname[0] is the path string
        except AttributeError:
            fnameasstring = ''

    def loadBeamdata(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Beamdata file', '')
        try:
            beamdata = loadBeamdata(fname[0])
            self.facility.setBeamdata(beamdata)
        except:
            print "Bad beamdata file!"

    def saveTwiss(self):
        fname = QFileDialog.getSaveFileName(self, 'Save Twiss file', '')
        try:
            saveTwiss(fname[0],self.facility.getTwiss()) # fname[0] is the path string
        except AttributeError:
            fnameasstring = ''

    def loadTwiss(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Twiss file', '')
        try:
            twiss = loadTwiss(fname[0])
            self.facility.setTwiss(twiss)
        except:
            print "Bad twiss file!"

    def generateMultipart(self):
        nbrOfParticles = self.getNbrOfParticles()
        twiss = self.getTwissFromInput()
        self.facility.generateMultipart(nbrOfParticles, twiss)
        return

    def generateGridpart(self):
        nbrOfParticles = self.getNbrOfParticles()
        self.facility.generateGridpart(nbrOfParticles)
        return

    def getNbrOfParticles(self):
        self.valueOfNbrOfParticles = self.enterNbrOfParticles.text()
        try:
            nbrOfParticles = int(self.valueOfNbrOfParticles)
        except:
            print "nbrOfParticles not integer, set to 10 instead"
            nbrOfParticles = 10
        return nbrOfParticles

    # [beta, rf_lambda, m, q, E, nbrOfParticles]
    def getBeamdataFromInput(self):
        self.valueOfBeta = self.enterBeta.text()
        self.valueOfLambda = self.enterLambda.text()
        self.valueOfMass = self.enterMass.text()
        self.valueOfCharge = self.enterCharge.text()
        self.valueOfEnergy = self.enterEnergy.text()
        self.valueOfNbrOfParticles = self.enterNbrOfParticles.text()
        self.valueOfCurrent = self.enterCurrent.text()

        try:
            beta = float(self.valueOfBeta)
        except:
            print "Beta not a number, set to 0.0 instead"
            beta = 0.0
        try:
            rf_lambda = float(self.valueOfLambda)
        except:
            print "Rf_lamda not a number, set to 0.0 instead"
            rf_lambda = 0.0
        try:
            m = float(self.valueOfMass)
        except:
            print "Mass not a number, set to 1.0 instead"
            m = 1.0
        try:
            q = float(self.valueOfCharge)
        except:
            print "Charge not a number, set to 0.0 instead"
            q = 0.0
        try:
            E = float(self.valueOfEnergy)
        except:
            print "Energy not a number, set to 0.0 instead"
            E = 0.0
        try:
            nbrOfParticles = int(self.valueOfNbrOfParticles)
        except:
            print "NbrOfParticles not a number, set to 1 instead"
            nbrOfParticles = 1
        try:
            I = float(self.valueOfCurrent)
        except:
            print "Current not a number, set to 0 instead"
            I = 0
        beamdata = [beta, rf_lambda, m, q, E, nbrOfParticles, I]
        return beamdata

    def useBeamdataInput(self):
        beamdata = self.getBeamdataFromInput()
        self.facility.setBeamdata(beamdata)

    def writeInBeamdataFields(self):
        beamdata = self.facility.getBeamdata()
        self.enterBeta.setText(str(beamdata[0]))
        self.enterLambda.setText(str(beamdata[1]))
        self.enterMass.setText(str(beamdata[2]))
        self.enterCharge.setText(str(beamdata[3]))
        self.enterEnergy.setText(str(beamdata[4]))
        self.enterNbrOfParticles.setText(str(beamdata[5]))
        self.enterCurrent.setText(str(beamdata[6]))

    def getTwissFromInput(self):
        self.valueOfAlpha_x = self.enterAlpha_x.text()
        self.valueOfBeta_x = self.enterBeta_x.text()
        self.valueOfEpsilon_rms_x = self.enterEpsilon_rms_x.text()
        self.valueOfAlpha_y = self.enterAlpha_y.text()
        self.valueOfBeta_y = self.enterBeta_y.text()
        self.valueOfEpsilon_rms_y = self.enterEpsilon_rms_y.text()
        self.valueOfAlpha_z = self.enterAlpha_z.text()
        self.valueOfBeta_z = self.enterBeta_z.text()
        self.valueOfEpsilon_rms_z = self.enterEpsilon_rms_z.text()
        try:
            Alpha_x = float(self.valueOfAlpha_x)
        except:
            print "Alpha_x not a number, set to 0.0 instead"
            Alpha_x = 0.0
        try:
            Beta_x = float(self.valueOfBeta_x)
        except:
            print "Beta_x not a number, set to 0.0 instead"
            Beta_x = 0.0
        try:
            Epsilon_rms_x = float(self.valueOfEpsilon_rms_x)
        except:
            print "Epsilon_rms_x not a number, set to 0.0 instead"
            Epsilon_rms_x = 0.0
        try:
            Alpha_y = float(self.valueOfAlpha_y)
        except:
            print "Alpha_y not a number, set to 0.0 instead"
            Alpha_y = 0.0
        try:
            Beta_y = float(self.valueOfBeta_y)
        except:
            print "Beta_y not a number, set to 0.0 instead"
            Beta_y = 0.0
        try:
            Epsilon_rms_y = float(self.valueOfEpsilon_rms_y)
        except:
            print "Epsilon_rms_y not a number, set to 0.0 instead"
            Epsilon_rms_y = 0.0
        try:
            Alpha_z = float(self.valueOfAlpha_z)
        except:
            print "Alpha_z not a number, set to 0.0 instead"
            Alpha_z = 0.0
        try:
            Beta_z = float(self.valueOfBeta_z)
        except:
            print "Beta_z not a number, set to 0.0 instead"
            Beta_z = 0.0
        try:
            Epsilon_rms_z = float(self.valueOfEpsilon_rms_z)
        except:
            print "Epsilon_rms_z not a number, set to 0.0 instead"
            Epsilon_rms_z = 0.0
        twiss = [Alpha_x, Beta_x, Epsilon_rms_x, Alpha_y, Beta_y, Epsilon_rms_y, Alpha_z, Beta_z, Epsilon_rms_z]
        return twiss

    def useTwissInput(self):
        twiss = self.getTwissFromInput()
        self.facility.setTwiss(twiss)

    def writeInTwissFields(self):
        twiss = self.facility.getTwiss()
        self.enterAlpha_x.setText(str(twiss[0]))
        self.enterBeta_x.setText(str(twiss[1]))
        self.enterEpsilon_rms_x.setText(str(twiss[2]))
        self.enterAlpha_y.setText(str(twiss[3]))
        self.enterBeta_y.setText(str(twiss[4]))
        self.enterEpsilon_rms_y.setText(str(twiss[5]))
        self.enterAlpha_z.setText(str(twiss[6]))
        self.enterBeta_z.setText(str(twiss[7]))
        self.enterEpsilon_rms_z.setText(str(twiss[8]))

    def saveMultipart(self):
        fname = QFileDialog.getSaveFileName(self, 'Save Multipart file', '')
        try:
            saveMultipart(fname[0],self.facility.getMultipart()) # fname[0] is the path string
        except AttributeError:
            fnameasstring = ''

    def loadMultipart(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Multipart file', '')
        try:
            multipart = loadMultipart(fname[0])
            self.facility.setMultipart(multipart)
        except Exception,e: print str(e) # prints the entire error
        #except:
        #    print "Bad multipart file!" + str(sys.exc_info()[-1].tb_lineno)
    

class LatticeEditor(QGroupBox):
    '''
    Widget for editing the lattice
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)

        self.parent = parent
        self.setTitle("Lattice Editor")

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        # load and save
        loadButton = QPushButton("Load Lattice")
        loadButton.clicked.connect(self.loadLattice)
        grid.addWidget(loadButton,0,0)

        saveButton = QPushButton("Save Lattice")
        saveButton.clicked.connect(self.saveLattice)
        grid.addWidget(saveButton,0,1)

        deleteButton = QPushButton("Delete Lattice")
        deleteButton.clicked.connect(self.deleteLattice)
        grid.addWidget(deleteButton,0,2)

        # SC and splits
        self.textNbrOfSplits = QLabel("# of Splits:")
        grid.addWidget(self.textNbrOfSplits, 1, 4)
        self.textNbrOfSplits.hide()
        
        self.enterNbrOfSplits = QLineEdit()
        grid.addWidget(self.enterNbrOfSplits, 1, 5)
        self.enterNbrOfSplits.hide()

        # SC radio buttons (they come after since they change the split buttons)
        self.textSC = QLabel("Space Charge:")
        grid.addWidget(self.textSC, 1, 0)
        
        self.radioGroup = QButtonGroup(self)

        self.radioNoSC = QRadioButton("None")
        self.radioGroup.addButton(self.radioNoSC)
        self.radioNoSC.toggled.connect(self.radioNoSC_clicked)
        grid.addWidget(self.radioNoSC, 1, 1)
        self.radioNoSC.toggle()

        self.radioSC1 = QRadioButton("Allen's SC")
        self.radioGroup.addButton(self.radioSC1)
        self.radioSC1.toggled.connect(self.radioSC1_clicked)
        grid.addWidget(self.radioSC1, 1, 2)

        self.radioSC2 = QRadioButton("Elliptic integral")
        self.radioGroup.addButton(self.radioSC2)
        self.radioSC2.toggled.connect(self.radioSC2_clicked)
        grid.addWidget(self.radioSC2, 1, 3)

        redoLatticeButton = QPushButton("Redo Lattice")
        redoLatticeButton.clicked.connect(self.redoLattice)
        grid.addWidget(redoLatticeButton,1,6)

        

        self.selectedElement = "Drift"

        self.elementSelector = QComboBox(self)
        self.elementSelector.addItem("Drift")
        self.elementSelector.addItem("Dipole")
        self.elementSelector.addItem("Quadrupole")
        self.elementSelector.addItem("LieDrift")
        self.elementSelector.addItem("Sextupole")
        self.elementSelector.addItem("Octupole")
        self.elementSelector.addItem("Sextupolerel")
        self.elementSelector.addItem("Sextupolekin")
        self.elementSelector.addItem("Octupolekin")
        self.elementSelector.addItem("Sextupolemat")
        self.elementSelector.addItem("Sextupolematema")
        self.elementSelector.addItem("Rotation")
        self.elementSelector.addItem("RF-Cavity") #  (Not implemented/working)
        self.elementSelector.addItem("Cavity Matrix") #  (need more work), bug now is that z' always goes to around -1 when it should be cetered around 0)
        self.elementSelector.addItem("Higher order element") #  (Not implemented/working)
        grid.addWidget(self.elementSelector, 2, 0)

        self.elementSelector.activated[str].connect(self.activatedElementSelector)
        self.selectedElement = "Drift"
        
        self.textName = QLabel("Name:")
        grid.addWidget(self.textName, 3, 0)

        self.enterName = QLineEdit()
        grid.addWidget(self.enterName, 3, 1)
        name = self.enterName.text()

        self.textL = QLabel("L [m]:")
        grid.addWidget(self.textL, 4, 0)
        
        self.enterL = QLineEdit()
        grid.addWidget(self.enterL, 4, 1)
        valueOfL = self.enterL.text()

        ## Quad and sextupole
        self.textK = QLabel("K [sqrt(T*C*s/kg)]:") # see calculations 2016-03-18 and TraceWin man s 102
        grid.addWidget(self.textK, 5, 0)
        self.textK.hide()

        self.enterK = QLineEdit()
        grid.addWidget(self.enterK, 5, 1)
        self.enterK.hide()

        ## Dipole
        self.textRho = QLabel("Rho [m]:")
        grid.addWidget(self.textRho, 4, 0)
        self.textRho.hide()

        self.enterRho = QLineEdit()
        grid.addWidget(self.enterRho, 4, 1)
        self.enterRho.hide()

        self.textAlpha = QLabel("Alpha [unitless angle]:")
        grid.addWidget(self.textAlpha, 5, 0)
        self.textAlpha.hide()

        self.enterAlpha = QLineEdit()
        grid.addWidget(self.enterAlpha, 5, 1)
        self.enterAlpha.hide()

        self.textn = QLabel("n []:")
        grid.addWidget(self.textn, 6, 0)
        self.textn.hide()

        self.entern = QLineEdit()
        grid.addWidget(self.entern, 6, 1)
        self.entern.hide()

        ## Sextupole
        self.textsK = QLabel("K [m^-3]:") # see calculations 2016-03-18 and TraceWin man s 102
        grid.addWidget(self.textsK, 5, 0)
        self.textsK.hide()

        self.entersK = QLineEdit()
        grid.addWidget(self.entersK, 5, 1)
        self.entersK.hide()

        ## Octupole
        self.textoK = QLabel("K [m^-4]:") # see calculations 2016-03-18 and TraceWin man s 102
        grid.addWidget(self.textoK, 5, 0)
        self.textoK.hide()

        self.enteroK = QLineEdit()
        grid.addWidget(self.enteroK, 5, 1)
        self.enteroK.hide()

        ## Order for Liealgelements
        self.textOrder = QLabel("Order []:")
        grid.addWidget(self.textOrder, 6, 0)
        self.textOrder.hide()

        self.enterOrder = QLineEdit()
        grid.addWidget(self.enterOrder, 6, 1)
        self.enterOrder.hide()

        ## Rotation
        self.textnu_x = QLabel("nu_x:")
        grid.addWidget(self.textnu_x, 4, 0)
        self.textnu_x.hide()

        self.enternu_x = QLineEdit()
        grid.addWidget(self.enternu_x, 4, 1)
        self.enternu_x.hide()

        self.textnu_y = QLabel("nu_y:")
        grid.addWidget(self.textnu_y, 5, 0)
        self.textnu_y.hide()

        self.enternu_y = QLineEdit()
        grid.addWidget(self.enternu_y, 5, 1)
        self.enternu_y.hide()

        ## Cavity Matrix
        self.texta = QLabel("a [m]:")
        grid.addWidget(self.texta, 5, 0)
        self.texta.hide()

        self.entera = QLineEdit()
        grid.addWidget(self.entera, 5, 1)
        self.entera.hide()

        self.textEfield_0 = QLabel("Efield_0 [V/m]:")
        grid.addWidget(self.textEfield_0, 6, 0)
        self.textEfield_0.hide()

        self.enterEfield_0 = QLineEdit()
        grid.addWidget(self.enterEfield_0, 6, 1)
        self.enterEfield_0.hide()

        self.textphi_0 = QLabel("phi_0 [rad]:")
        grid.addWidget(self.textphi_0, 7, 0)
        self.textphi_0.hide()

        self.enterphi_0 = QLineEdit()
        grid.addWidget(self.enterphi_0, 7, 1)
        self.enterphi_0.hide()

        ## Create Element
        createElementButton = QPushButton("Create Element")
        createElementButton.clicked.connect(self.createElement) # here arguments should be passed
        grid.addWidget(createElementButton, 8,1)

        ## Print Matrices
        printMatricesButton = QPushButton("Print Matrices")
        printMatricesButton.clicked.connect(self.printMatrices)
        grid.addWidget(printMatricesButton, 8,2)

    def radioNoSC_clicked(self, enabled):
        if enabled:
            self.spaceCharge = 0
            self.textNbrOfSplits.hide()
            self.enterNbrOfSplits.hide()

    def radioSC1_clicked(self, enabled):
        if enabled:
            self.spaceCharge = 1
            self.textNbrOfSplits.show()
            self.enterNbrOfSplits.show()

    def radioSC2_clicked(self, enabled):
        if enabled:
            self.spaceCharge = 2
            self.textNbrOfSplits.show()
            self.enterNbrOfSplits.show()

    def redoLattice(self):
        # send self.spaceCharge to facility which will send it to accelerator
        if self.spaceCharge == 0:
            self.nbrOfSplits = 1
        else:
            self.valueOfNbrOfSplits = self.enterNbrOfSplits.text()
            try:
                self.nbrOfSplits = int(self.valueOfNbrOfSplits)
            except:
                print "# of Splits not a number, set to 1 instead!"
                self.nbrOfSplits = 1
        self.facility.setSpaceChargeOnAndSplits(self.spaceCharge, self.nbrOfSplits)

    def activatedElementSelector(self, text):
        print text
        self.selectedElement = text

        # Clear away all boxes
        self.textK.hide()
        self.enterK.hide()
        self.textL.hide()
        self.enterL.hide()
        self.textRho.hide()
        self.enterRho.hide()
        self.textAlpha.hide()
        self.enterAlpha.hide()
        self.textn.hide()
        self.entern.hide()
        self.textsK.hide()
        self.entersK.hide()
        self.textoK.hide()
        self.enteroK.hide()
        self.textOrder.hide()
        self.enterOrder.hide()
        self.textnu_x.hide()
        self.enternu_x.hide()
        self.textnu_y.hide()
        self.enternu_y.hide()
        self.texta.hide()
        self.entera.hide()
        self.textEfield_0.hide()
        self.enterEfield_0.hide()
        self.textphi_0.hide()
        self.enterphi_0.hide()

        if text == "Drift":
            self.textL.show()
            self.enterL.show()
        elif text == "Quadrupole":
            self.textK.show()
            self.enterK.show()
            self.textL.show()
            self.enterL.show()
        elif text == "Dipole":
            self.textRho.show()
            self.enterRho.show()
            self.textAlpha.show()
            self.enterAlpha.show()
            self.textn.show()
            self.entern.show()
        if text == "LieDrift":
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Sextupole":
            self.textsK.show()
            self.entersK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Octupole":
            self.textoK.show()
            self.enteroK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Sextupolerel":
            self.textsK.show()
            self.entersK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Sextupolekin":
            self.textsK.show()
            self.entersK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Octupolekin":
            self.textoK.show()
            self.enteroK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()
        elif text == "Sextupolemat":
            self.textsK.show()
            self.entersK.show()
            self.textL.show()
            self.enterL.show()
        elif text == "Sextupolematema":
            self.textsK.show()
            self.entersK.show()
            self.textL.show()
            self.enterL.show()
        elif text == "Rotation":
            self.textnu_x.show()
            self.enternu_x.show()
            self.textnu_y.show()
            self.enternu_y.show()
        elif text == "RF-Cavity":
            self.textL.show()
            self.enterL.show()
        elif text == "Cavity Matrix":
            self.textL.show()
            self.enterL.show()
            self.texta.show()
            self.entera.show()
            self.textEfield_0.show()
            self.enterEfield_0.show()
            self.textphi_0.show()
            self.enterphi_0.show()

    def createElement(self):
        ## Take the inputs in fields
        name = self.enterName.text()

        if not self.enterL.isHidden():
            valueOfL = self.enterL.text()
            try:
                L = float(valueOfL)
            except:
                print "Not a number! L set to 0.0"
                L = 0.0

        if not self.enterK.isHidden():
            valueOfK = self.enterK.text()
            try:
                K = float(valueOfK)
            except:
                print "Not a number! K set to 0.0"
                K = 0.0

        ## Dipole
        if not self.enterRho.isHidden():
            valueOfRho = self.enterRho.text()
            try:
                Rho = float(valueOfRho)
            except:
                print "Not a number! Rho set to 0.0"
                Rho = 0.0

        if not self.enterAlpha.isHidden():
            valueOfAlpha = self.enterAlpha.text()
            try:
                Alpha = float(valueOfAlpha)
            except:
                print "Not a number! Alpha set to 0.0"
                Alpha = 0.0

        if not self.entern.isHidden():
            valueOfn = self.entern.text()
            try:
                n = float(valueOfn)
            except:
                print "Not a number! n set to 0.0"
                n = 0.0

        ## Sextupole
        if not self.entersK.isHidden():
            valueOfsK = self.entersK.text()
            try:
                K = float(valueOfsK)
            except:
                print "Not a number! K set to 0.0"
                K = 0.0

        ## Octupole
        if not self.enteroK.isHidden():
            valueOfoK = self.enteroK.text()
            try:
                K = float(valueOfoK)
            except:
                print "Not a number! K set to 0.0"
                K = 0.0

        ## Order for LieAlg elements
        if not self.enterOrder.isHidden():
            valueOfOrder = self.enterOrder.text()
            try:
                Order = int(valueOfOrder)
            except:
                print "Not a number! Order set to 5"
                Order = 5

        ## Rotation
        if not self.enternu_x.isHidden():
            valueOfnu_x = self.enternu_x.text()
            try:
                nu_x = float(valueOfnu_x)
            except:
                print "Not a number! nu_x set to 0.246"
                nu_x = 0.246

        if not self.enternu_y.isHidden():
            valueOfnu_y = self.enternu_y.text()
            try:
                nu_y = float(valueOfnu_y)
            except:
                print "Not a number! nu_y set to 0.246"
                nu_y = 0.246

        ## Cavity Matrix
        if not self.entera.isHidden():
            valueOfa = self.entera.text()
            try:
                a = float(valueOfa)
            except:
                print "Not a number! a set to 0.76553527627 (gives k = pi, which is nice to combine with L = 1)"
                a = 0.76553527627

        if not self.enterEfield_0.isHidden():
            valueOfEfield_0 = self.enterEfield_0.text()
            try:
                Efield_0 = float(valueOfEfield_0)
            except:
                print "Not a number! Efield_0 set to 100"
                Efield_0 = 100

        if not self.enterphi_0.isHidden():
            valueOfphi_0 = self.enterphi_0.text()
            try:
                phi_0 = float(valueOfphi_0)
            except:
                print "Not a number! phi_0 set to 0"
                phi_0 = 0

        ## Make new elements
        if self.selectedElement == "Drift":
            self.facility.createDrift(name, L)
        elif self.selectedElement == "Dipole":
            self.facility.createDipole(name, Rho, Alpha, n)
        elif self.selectedElement == "Quadrupole":
            self.facility.createQuadrupole(name, K, L)
        elif self.selectedElement == "LieDrift":
            self.facility.createLieDrift(name, L, Order)
        elif self.selectedElement == "Sextupole":
            self.facility.createSextupole(name, K, L, Order)
        elif self.selectedElement == "Octupole":
            self.facility.createOctupole(name, K, L, Order)
        elif self.selectedElement == "Sextupolerel":
            self.facility.createSextupolerel(name, K, L, Order)
        elif self.selectedElement == "Sextupolekin":
            self.facility.createSextupolekin(name, K, L, Order)
        elif self.selectedElement == "Octupolekin":
            self.facility.createOctupolekin(name, K, L, Order)
        elif self.selectedElement == "Sextupolemat":
            self.facility.createSextupolemat(name, K, L)
        elif self.selectedElement == "Sextupolematema":
            self.facility.createSextupolematema(name, K, L)
        elif self.selectedElement == "Rotation":
            self.facility.createRotation(name, nu_x, nu_y)
        elif self.selectedElement == "RF-Cavity":
            Oscillations = 2
            AmplitudeA = 0
            AmplitudeB = 30 # 30 MeV / m
            E_0 = AmplitudeB
            Sigma = 1
            P = 3
            Ezofs = [Oscillations, AmplitudeA, AmplitudeB, E_0, Sigma, P]
            self.facility.createCavity(name, L, Ezofs)
        elif self.selectedElement == "Cavity Matrix":
            self.facility.createCavityMatrix(name, L, a, Efield_0, phi_0)
        else:
            return
        # Higher order (specify order)
        # Cavity
        

        self.parent.latticeoverview.initializeGL() # update the paint lattice in overview
        self.parent.parent.widget.latticeoverview.s_pressed = 1 # zoom out
        self.parent.parent.widget.latticeoverview.d_pressed = 1 # move with lattice
        self.parent.latticeoverview.mousePressEvent("bla") # repaint by setting focus

    def saveLattice(self):
        fname = QFileDialog.getSaveFileName(self, 'Save Lattice file', '')
        try:
            saveLattice(fname[0],self.facility.getLattice()) # fname[0] is the path string
        except AttributeError:
            fnameasstring = ''

    def loadLattice(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Lattice file', '')
        try:
            lattice = loadLattice(fname[0], self.facility)
            if lattice == 0:
                raise ValueError("loadLattice returned only 0")
            self.facility.setLattice(lattice)
        except:
            print "Lattice load failed!" + str(sys.exc_info()[-1].tb_lineno)
        self.parent.latticeoverview.initializeGL()
        self.parent.latticeoverview.paintGL()

    def deleteLattice(self):
        item, ok = QInputDialog.getItem(self, "Are you sure you really want to delete the lattice?", "options", ("Yes", "No"), 1, False)
        if item == "Yes":
            self.facility.deleteLattice()
            self.parent.latticeoverview.initializeGL()
            self.parent.parent.widget.latticeoverview.s_pressed = 1 # zoom out
            self.parent.parent.widget.latticeoverview.a_pressed = 1 # move with lattice
            self.parent.latticeoverview.mousePressEvent("bla") # repaint by setting focus
            self.parent.evalwidget.updateLaps()
        
    def printMatrices(self):
        print self.facility.printMatrices()
    
class EvalWidget(QWidget):
    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)
        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        EvalButton = QPushButton("Evaluate!")
        EvalButton.setMinimumSize(300,100)
        EvalButton.setStyleSheet("background-color:green")
        EvalButton.clicked.connect(self.evaluate)
        grid.addWidget(EvalButton,0,0)

        self.SaveResultButton = QPushButton("Save Results")
        self.SaveResultButton.clicked.connect(self.saveResults)
        grid.addWidget(self.SaveResultButton,0,1)
        self.SaveResultButton.hide()

        EvalButton = QPushButton("Evaluate with profiling!")
        EvalButton.setMinimumSize(200,50)
        EvalButton.setStyleSheet("background-color:yellow")
        EvalButton.clicked.connect(self.evaluateWithProfiling)
        grid.addWidget(EvalButton,0,2)

        self.textLaps = QLabel("# of laps:") # see calculations 2016-03-18 and TraceWin man s 102
        grid.addWidget(self.textLaps, 0, 3)

        self.enterLaps = QLineEdit()
        grid.addWidget(self.enterLaps, 0, 4)
        self.enterLaps.setText(str(self.facility.getLaps()))

    def evaluate(self):
        self.setLaps()
        self.facility.evaluate()
        self.facility.plotAfterEval()
        self.SaveResultButton.show()
        return

    def saveResults(self):
        resultmultipart, resultenvelope, resulttwiss, resultenvlist = self.facility.getResults()

        fnamemul = QFileDialog.getSaveFileName(self, 'Save Multipart file', '')
        try:
            saveMultipart(fnamemul[0],resultmultipart) # fname[0] is the path string
        except AttributeError:
            print "Could not save multipart results!"

        fnameenvelope = QFileDialog.getSaveFileName(self, 'Save Envelope file', '')
        try:
            saveEnvelope(fnameenvelope[0],resultenvelope) # fname[0] is the path string
        except AttributeError:
            print "Could not save envelope results!"
        
        fnametwiss = QFileDialog.getSaveFileName(self, 'Save Twiss file', '')
        try:
            saveTwiss(fnametwiss[0],resulttwiss) # fname[0] is the path string
        except AttributeError:
            print "Could not save twiss results!"

        fnameenvelopelist = QFileDialog.getSaveFileName(self, 'Save Envelopelist file', '')
        try:
            saveEnvelope(fnameenvelopelist[0],resultenvlist) # fname[0] is the path string
        except AttributeError:
            print "Could not save envelopelist results!"

        return

    def evaluateWithProfiling(self):
        self.setLaps()
        self.facility.evaluateWithProfiling()
        self.facility.plotAfterEval()
        self.SaveResultButton.show()
        return

    def setLaps(self):
        valueOfLaps = self.enterLaps.text()
        try:
            laps = int(valueOfLaps)
        except:
            print "Laps not an integer, set to 1 instead."
            laps = 1
        self.facility.setLaps(laps)

    def updateLaps(self):
        self.enterLaps.setText(str(self.facility.getLaps()))

# layout manager (aranges the different widgets)
class FormWidget(QWidget):
    def __init__(self, parent, facility):        
        #  Initialize GLUT and process user parameters
        glutInit(sys.argv);

        #  Request double buffered true color window with Z-buffer
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

        #  Enable Z-buffer depth test
        glEnable(GL_DEPTH_TEST)

        super(FormWidget, self).__init__(parent)

        self.parent = parent

        ## Layout
        self.layout = QGridLayout()

        # Lattice overview
        self.latticeoverview = LatticeOverviewWidget(facility, self)
        self.latticeoverview.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout.addWidget(self.latticeoverview, 0, 0, 3, 1)

        self.latticeoverview.setFocus() # starts with focus
        self.latticeoverview.mousePressEvent("bla") # for some reason I need a parameter in this overloaded function (set it to bla sends 0 as bla)

        # Beam editor
        self.beameditor = BeamEditor(self, facility)
        self.beameditor.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.layout.addWidget(self.beameditor, 0, 1)

        # Lattice editor
        self.latticeeditor = LatticeEditor(self, facility)
        self.latticeeditor.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.layout.addWidget(self.latticeeditor, 1, 1)

        # Evaluate
        self.evalwidget = EvalWidget(self, facility)
        self.evalwidget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.layout.addWidget(self.evalwidget, 2, 1)

        self.setLayout(self.layout) 

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)
        self.facility = Facility()

        self.widget = FormWidget(self, self.facility)
        self.setCentralWidget(self.widget)

        # Move to top left
        qr = self.frameGeometry()
        self.move(qr.topLeft())

    def keyPressEvent(self, e):
        
        if e.key() == Qt.Key_Escape:
            self.widget.latticeoverview.spin()

        # keys: w = 87, a = 65, s = 83, d = 68
        if e.key() not in range(256):
            return
        if chr(e.key()) == 'W':
            self.widget.latticeoverview.w_pressed = 1
        if chr(e.key()) == 'A':
            self.widget.latticeoverview.a_pressed = 1
        if chr(e.key()) == 'S':
            self.widget.latticeoverview.s_pressed = 1
        if chr(e.key()) == 'D':
            self.widget.latticeoverview.d_pressed = 1
        self.widget.latticeoverview.paintGL() # repaint
        self.widget.latticeoverview.updateGL() # show new painting

    def keyReleaseEvent(self, e):
        if e.key() not in range(256):
            return
        if chr(e.key()) == 'W':
            self.widget.latticeoverview.w_pressed = 0
        if chr(e.key()) == 'A':
            self.widget.latticeoverview.a_pressed = 0
        if chr(e.key()) == 'S':
            self.widget.latticeoverview.s_pressed = 0
        if chr(e.key()) == 'D':
            self.widget.latticeoverview.d_pressed = 0
        self.widget.latticeoverview.paintGL()
        self.widget.latticeoverview.updateGL()

    def resizeEvent(self, event):
        if self.widget is not None:
            self.widget.resizeEvent(event)

    


if __name__ == '__main__':
    app = QApplication(['DAT Widget Interface'])
    window = DATWidgetInterface()
    window.show()
    app.exec_()