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
from PyQt5.QtWidgets import QMainWindow, QAction, qApp, QApplication
from IOHandler import saveLattice, loadLattice, loadSummer2015Formatzasx, saveBeamdata, loadBeamdata, saveTwiss, loadTwiss, saveMultipart, loadMultipart # loadLattice, saveLattice
from accelerator import Lattice
from numpy import array
import numpy as np
from scipy import *
from facility import *
import re
from particleFactory import envelopeFromMultipart
from plotting import plotEverything

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

        self.up = np.array([0.0, 1.0, 0.0])
        tempCross = np.cross(self.up, self.cameraDirection)
        self.cameraRight = tempCross/np.linalg.norm(tempCross)
        self.cameraUp = np.cross(self.cameraDirection, self.cameraRight)

        # extract from lattice
        self.facility = facility
        self.lattice = self.facility.getLattice()

    def loadLattice(self):
        #self.lattice = self.facility.getLattice()
        #lattticeString = self.lattice.printLattice()
        lattticeString = self.facility.printLattice()
        self.elements = []
        nextWillBeL = 0
        for line in lattticeString.split():
            if nextWillBeL:
                self.elements.append([tempword, float(line)])
                nextWillBeL = 0
            if line == "drift" or line == "dipole" or line == "quad" or line == "liealgelem" or line == "cavity":  
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
            elif elem[0] == "cavity":
                color = [1, 1, 0]

            length = elem[1]

            if elem[0] == "liealgelem":
                VtxArray, IdxArray, ClrArray = self.createHexaBlock(length, color)
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

        cameraSpeed = 0.5
        if self.w_pressed:
            self.cameraPos += self.cameraDirection*cameraSpeed
            self.w_pressed = 0
        if self.s_pressed:
            self.cameraPos -= self.cameraDirection*cameraSpeed
            self.s_pressed = 0
        if self.a_pressed:
            tempCross = np.cross(self.cameraDirection, self.up)
            normTempCross = tempCross/np.linalg.norm(tempCross)
            self.cameraPos -= normTempCross*cameraSpeed
            self.a_pressed = 0
        if self.d_pressed:
            tempCross = np.cross(self.cameraDirection, self.up)
            normTempCross = tempCross/np.linalg.norm(tempCross)
            self.cameraPos += normTempCross*cameraSpeed
            self.d_pressed = 0

        self.zsofar = 0
        for block in self.array_of_blocks:
            self.elemPaint(block)

        
    def elemPaint(self, elem):
        glLoadIdentity()
        glTranslate(self.zsofar, 0.0, -5.0)
        self.zsofar = self.zsofar + elem[3]
        glTranslate(self.cameraPos[0], self.cameraPos[1], self.cameraPos[2])
        glRotate(90, 0.0, 1.0, 0.0)
        elemtype = elem[4]
        if elemtype == "quad":
            glRotate(45, 0.0, 0.0, 1.0)
        elif elemtype == "liealgelem":
            glTranslate(0.5, 0.5, 0.0)
        glTranslate(-0.5, -0.5, -0.5)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointerf(elem[0])
        glColorPointerf(elem[2])
        if elemtype == "liealgelem":
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

    # If one clicks on the openGL panel it gets the focus, also this will be called if focus on opengl is required
    def mousePressEvent(self, bla):
        self.setFocus()
        self.updateGL()


class BeamEditor(QWidget):
    '''
    Widget for editing the input beam
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)

        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        self.textBeamEditor = QLabel("Beam Editor")
        grid.addWidget(self.textBeamEditor, 0, 0)

        ## Beamdata
        # beta
        self.textBeta = QLabel("Relativistic beta:")
        grid.addWidget(self.textBeta, 1, 0)
        
        self.enterBeta = QLineEdit()
        grid.addWidget(self.enterBeta, 1, 1)
        valueOfBeta = self.enterBeta.text()

        # rf_lambda
        self.textLambda = QLabel("RF-Lambda [m]:")
        grid.addWidget(self.textLambda, 2, 0)
        
        self.enterLambda = QLineEdit()
        grid.addWidget(self.enterLambda, 2, 1)
        valueOfLambda = self.enterLambda.text()

        # m
        self.textMass = QLabel("Mass of a particle:")
        grid.addWidget(self.textMass, 3, 0)
        
        self.enterMass = QLineEdit() # should be a ComboBox with m_p, m_e
        grid.addWidget(self.enterMass, 3, 1)
        valueOfMass = self.enterMass.text()

        # q
        self.textCharge = QLabel("Charge of a particle:")
        grid.addWidget(self.textCharge, 4, 0)
        
        self.enterCharge = QLineEdit() # should be a comboBox with e and -e
        grid.addWidget(self.enterCharge, 4, 1)
        valueOfCharge = self.enterCharge.text()

        # E
        self.textEnergy = QLabel("Energy:")
        grid.addWidget(self.textEnergy, 5, 0)
        
        self.enterEnergy = QLineEdit() # should be a QLabel and show the calculated value from beta
        grid.addWidget(self.enterEnergy, 5, 1)
        valueOfEnergy = self.enterEnergy.text()

        # nbrOfParticles
        self.textNbrOfParticles = QLabel("# of particles:")
        grid.addWidget(self.textNbrOfParticles, 6, 0)
        
        self.enterNbrOfParticles = QLineEdit()
        grid.addWidget(self.enterNbrOfParticles, 6, 1)
        valueOfNbrOfParticles = self.enterNbrOfParticles.text()

        saveBeamdataButton = QPushButton("Save Beamdata")
        saveBeamdataButton.clicked.connect(self.saveBeamdata)
        grid.addWidget(saveBeamdataButton,10,1)

        ## Twiss
        # twiss comes as [alpha_x, beta_x, epsilon_rms_x, alpha_y, beta_y, epsilon_rms_y, alpha_z, beta_z, epsilon_rms_z]
        # alpha_x
        self.textAlpha_x = QLabel("Alpha_x:")
        grid.addWidget(self.textAlpha_x, 1, 3)
        
        self.enterAlpha_x = QLineEdit()
        grid.addWidget(self.enterAlpha_x, 1, 4)
        valueOfAlpha_x = self.enterAlpha_x.text()

        # beta_x
        self.textBeta_x = QLabel("Beta_x:")
        grid.addWidget(self.textBeta_x, 2, 3)
        
        self.enterBeta_x = QLineEdit()
        grid.addWidget(self.enterBeta_x, 2, 4)
        valueOfBeta_x = self.enterBeta_x.text()
        
        # epsilon_rms_x
        self.textEpsilon_rms_x = QLabel("Epsilon_rms_x:")
        grid.addWidget(self.textEpsilon_rms_x, 3, 3)
        
        self.enterEpsilon_rms_x = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_x, 3, 4)
        valueOfEpsilon_rms_x = self.enterEpsilon_rms_x.text()
        
        # alpha_y
        self.textAlpha_y = QLabel("Alpha_y:")
        grid.addWidget(self.textAlpha_y, 4, 3)
        
        self.enterAlpha_y = QLineEdit()
        grid.addWidget(self.enterAlpha_y, 4, 4)
        valueOfAlpha_y = self.enterAlpha_y.text()
        
        # beta_y
        self.textBeta_y = QLabel("Beta_y:")
        grid.addWidget(self.textBeta_y, 5, 3)
        
        self.enterBeta_y = QLineEdit()
        grid.addWidget(self.enterBeta_y, 5, 4)
        valueOfBeta_y = self.enterBeta_y.text()
        
        # epsilon_rms_y
        self.textEpsilon_rms_y = QLabel("Epsilon_rms_y:")
        grid.addWidget(self.textEpsilon_rms_y, 6, 3)
        
        self.enterEpsilon_rms_y = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_y, 6, 4)
        valueOfEpsilon_rms_y = self.enterEpsilon_rms_y.text()
        
        # alpha_z
        self.textAlpha_z = QLabel("Alpha_z:")
        grid.addWidget(self.textAlpha_z, 7, 3)
        
        self.enterAlpha_z = QLineEdit()
        grid.addWidget(self.enterAlpha_z, 7, 4)
        valueOfAlpha_z = self.enterAlpha_z.text()
        
        # beta_z
        self.textBeta_z = QLabel("Beta_z:")
        grid.addWidget(self.textBeta_z, 8, 3)
        
        self.enterBeta_z = QLineEdit()
        grid.addWidget(self.enterBeta_z, 8, 4)
        valueOfBeta_z = self.enterBeta_z.text()
        
        # epsilon_rms_z
        self.textEpsilon_rms_z = QLabel("Epsilon_rms_z:")
        grid.addWidget(self.textEpsilon_rms_z, 9, 3)
        
        self.enterEpsilon_rms_z = QLineEdit()
        grid.addWidget(self.enterEpsilon_rms_z, 9, 4)
        valueOfEpsilon_rms_z = self.enterEpsilon_rms_z.text()

        saveTwissButton = QPushButton("Save Twiss")
        saveTwissButton.clicked.connect(self.saveTwiss)
        grid.addWidget(saveTwissButton,10,4)

        generateMultipartButton = QPushButton("Generate Multiparticles")
        #generateMultipartButton.clicked.connect(self.saveMultipart)
        grid.addWidget(generateMultipartButton,9,5)

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
        except:
            print "Bad multipart file!"
    

class LatticeEditor(QWidget):
    '''
    Widget for editing the lattice
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)

        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        self.textLatticeEditor = QLabel("Lattice Editor")
        grid.addWidget(self.textLatticeEditor, 0, 0)

        loadButton = QPushButton("Load Lattice")
        loadButton.clicked.connect(self.loadLattice)
        grid.addWidget(loadButton,1,0)

        saveButton = QPushButton("Save Lattice")
        saveButton.clicked.connect(self.saveLattice)
        grid.addWidget(saveButton,1,1)

        self.selectedElement = "Drift"

        self.elementSelector = QComboBox(self)
        self.elementSelector.addItem("Drift")
        self.elementSelector.addItem("Dipole")
        self.elementSelector.addItem("Quadrupole")
        self.elementSelector.addItem("Sextupole")
        self.elementSelector.addItem("RF-cavity")
        self.elementSelector.addItem("Higher order element")
        grid.addWidget(self.elementSelector, 2, 0)

        self.elementSelector.activated[str].connect(self.activatedElementSelector)
        self.selectedElement = "Drift"
        
        self.textName = QLabel("Name:")
        grid.addWidget(self.textName, 3, 0)

        self.enterName = QLineEdit()
        grid.addWidget(self.enterName, 3, 1)
        name = self.enterName.text()

        self.textL = QLabel("L:")
        grid.addWidget(self.textL, 4, 0)
        
        self.enterL = QLineEdit()
        grid.addWidget(self.enterL, 4, 1)
        valueOfL = self.enterL.text()

        ## Quad and sextupole
        self.textK = QLabel("K:")
        grid.addWidget(self.textK, 5, 0)
        self.textK.hide()

        self.enterK = QLineEdit()
        grid.addWidget(self.enterK, 5, 1)
        self.enterK.hide()

        ## Dipole
        self.textRho = QLabel("Rho:")
        grid.addWidget(self.textRho, 4, 0)
        self.textRho.hide()

        self.enterRho = QLineEdit()
        grid.addWidget(self.enterRho, 4, 1)
        self.enterRho.hide()

        self.textAlpha = QLabel("Alpha:")
        grid.addWidget(self.textAlpha, 5, 0)
        self.textAlpha.hide()

        self.enterAlpha = QLineEdit()
        grid.addWidget(self.enterAlpha, 5, 1)
        self.enterAlpha.hide()

        self.textn = QLabel("n:")
        grid.addWidget(self.textn, 6, 0)
        self.textn.hide()

        self.entern = QLineEdit()
        grid.addWidget(self.entern, 6, 1)
        self.entern.hide()

        ## Sextupole
        self.textOrder = QLabel("Order:")
        grid.addWidget(self.textOrder, 6, 0)
        self.textOrder.hide()

        self.enterOrder = QLineEdit()
        grid.addWidget(self.enterOrder, 6, 1)
        self.enterOrder.hide()

        createElementButton = QPushButton("Create Element")
        createElementButton.clicked.connect(self.createElement) # here arguments should be passed
        grid.addWidget(createElementButton, 7,1)

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
        self.textOrder.hide()
        self.enterOrder.hide()

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
        elif text == "Sextupole":
            self.textK.show()
            self.enterK.show()
            self.textL.show()
            self.enterL.show()
            self.textOrder.show()
            self.enterOrder.show()

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

        if not self.enterOrder.isHidden():
            valueOfOrder = self.enterOrder.text()
            try:
                Order = int(valueOfOrder)
            except:
                print "Not a number! Order set to 5"
                Order = 5

        ## Make new elements
        if self.selectedElement == "Drift":
            self.facility.createDrift(name, L)
        elif self.selectedElement == "Dipole":
            self.facility.createDipole(name, Rho, Alpha, n)
        elif self.selectedElement == "Quadrupole":
            self.facility.createQuadrupole(name, K, L)
        elif self.selectedElement == "Sextupole":
            self.facility.createSextupole(name, K, L, Order)
        else:
            return
        # Higher order (specify order)
        # Cavity
        

        self.parent.latticeoverview.initializeGL() # update the paint lattice in overview
        self.parent.parent.widget.latticeoverview.s_pressed = 1 # prepare a zoom out
        self.parent.parent.widget.latticeoverview.d_pressed = 1
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
            self.facility.setLattice(lattice)
        except:
            print "Lattice load failed!"
        self.parent.latticeoverview.initializeGL()
        self.parent.latticeoverview.paintGL()
    
class EvalWidget(QWidget):
    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)
        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        EvalButton = QPushButton("Evaluate!")
        EvalButton.setStyleSheet("background-color:green")
        EvalButton.clicked.connect(self.evaluate)
        grid.addWidget(EvalButton,0,0)

    def evaluate(self):
        multipart = self.facility.getMultipart()
        envelope = envelopeFromMultipart(multipart) #self.facility.getEnvelope()
        twiss = self.facility.getTwiss()
        resultmultipart, resultenvelope, resulttwiss = self.facility.evaluate(multipart, envelope, twiss)
        plotEverything(multipart, twiss, resultmultipart)
        return

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
        self.layout = QHBoxLayout(self)

        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.layout.addStretch(1)

        # Lattice overview
        self.latticeoverview = LatticeOverviewWidget(facility, self)
        self.layout.addWidget(self.latticeoverview)

        self.latticeoverview.setFocus() # starts with focus
        self.latticeoverview.mousePressEvent("bla") # for some reason I need a parameter in this overloaded function (set it to bla sends 0 as bla)

        ## Editor layout
        self.editorlayout = QVBoxLayout() # no need for self as param since layout will later set this as its child

        # Beam editor
        self.beameditor = BeamEditor(self, facility)
        self.editorlayout.addWidget(self.beameditor)

        # Lattice editor
        self.latticeeditor = LatticeEditor(self, facility)
        self.editorlayout.addWidget(self.latticeeditor)

        # Evaluate
        self.evalwidget = EvalWidget(self, facility)
        self.editorlayout.addWidget(self.evalwidget)

        ## More layout stuff
        self.layout.addLayout(self.editorlayout)

        self.setLayout(self.layout)

    def resizeEvent(self, event):
        newWidth = max(self.frameGeometry().width()-720,0) # 720 is magic number for getting the correct width
        self.latticeoverview.setGeometry(4,4,newWidth, self.latticeoverview.height())
        return     

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)
        self.facility = Facility()

        self.widget = FormWidget(self, self.facility)
        self.setCentralWidget(self.widget)

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