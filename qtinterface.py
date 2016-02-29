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
from IOHandler import parseLatticeString, saveLatticeString, loadLatticeString, loadSummer2015Formatzasx, loadTwiss, loadMultipart #, loadLattice, saveLattice
from accelerator import Lattice
from numpy import array
import numpy as np
from scipy import *
import time
from facility import *
import re

class LatticeOverviewWidget(QGLWidget):
    '''
    Widget for giving an overview of the lattice
    '''
    def __init__(self, facility, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.yRotDeg = 1.0
        self.setMinimumSize(600, 500)

        self.z = 2.0

        #self.timer = QtCore.QTimer(self)
        #self.timer.setInterval(20)
        #self.trigger = pyqtSignal()
        #print str(dir(self))
        #connect(self.timer, QtCore.SIGNAL('timeout()'), self.spin)
        #self.timer.start()

        self.w_pressed = 0
        self.a_pressed = 0
        self.s_pressed = 0
        self.d_pressed = 0

        self.radius = 1
        self.camX = sin(time.time()) * self.radius;
        self.camY = cos(time.time()) * self.radius;

        #cameraPos = np.array([0.0, 0.0, 3.0])
        #self.cameraPos = np.array([self.camX, self.camY, 0.0])
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

        # is mouse pressed
        #self.mousePressed = 0
        self.yRotDeg = 0.0

    def loadLattice(self):
        self.lattice = self.facility.getLattice()
        lattticeString = self.lattice.printLattice()
        self.elements = []
        nextWillBeL = 0
        for line in lattticeString.split():
            #print "line: " + line
            if nextWillBeL:
                self.elements.append([tempword, float(line)])
                nextWillBeL = 0
            if line == "drift" or line == "dipole" or line == "quad" or line == "liealgelem" or line == "cavity":  
                tempword = line
            if line == "L:":
                nextWillBeL = 1

        print "self.elements: \n" + str(self.elements)


    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0,  150))
        #self.initGeometry()

        self.bluecubeVtxArray, self.bluecubeIdxArray, self.bluecubeClrArray = self.bluecube(self.z)
        self.redcubeVtxArray, self.redcubeIdxArray, self.redcubeClrArray = self.redcube(self.z)

        greencolor = [0, 1, 0]
        self.greencubeVtxArray, self.greencubeIdxArray, self.greencubeClrArray = self.createGeomBlock(self.z, greencolor)
        self.greencube = [self.greencubeVtxArray, self.greencubeIdxArray, self.greencubeClrArray]

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
        #print "a: " + str(a)
        b = np.array([
                    [1.0, 0.0, 0.0, -position[0]],
                    [0.0, 1.0, 0.0, -position[1]],
                    [0.0, 0.0, 1.0, -position[2]],
                    [0.0, 0.0, 0.0, 1.0]
                    ])
        #print "b: " + str(b)
        return a*b

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        #glLoadIdentity()
        #gluLookAt(0.0, 105.0, 105.0,
        #  0.0, 0.0, 0.0,
        #  0.0, 1.0, 0.0)
        #glRotatef(self.yRotDeg, 0.0, 1.0, 0.0)

        self.camX = sin(time.time()) * self.radius;
        self.camY = cos(time.time()) * self.radius;

        cameraSpeed = 0.5
        if self.w_pressed:
            self.cameraPos += self.cameraDirection*cameraSpeed
            self.w_pressed = 0
            #print "Camera moved from w"
            #print "cameraPos: " + str(self.cameraPos)
        if self.s_pressed:
            self.cameraPos -= self.cameraDirection*cameraSpeed
            self.s_pressed = 0
            #print "Camera moved from s"
            #print "cameraPos: " + str(self.cameraPos)
        if self.a_pressed:
            #print "Camera moved from a"
            tempCross = np.cross(self.cameraDirection, self.up)
            normTempCross = tempCross/np.linalg.norm(tempCross)
            self.cameraPos -= normTempCross*cameraSpeed
            self.a_pressed = 0
            #print "cameraPos: " + str(self.cameraPos)
        if self.d_pressed:
            #print "Camera moved from d"
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
        #glRotate(self.yRotDeg, 90.0, 1.0, 0.0)
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
        

    def bluecube(self, z):
        bluecubeVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, z],
                 [1.0, 0.0, z],
                 [1.0, 1.0, z],
                 [0.0, 1.0, z]])
        bluecubeIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        bluecubeClrArray = np.zeros((8,3))
        bluecubeClrArray[:,2] = 1.0
        return bluecubeVtxArray, bluecubeIdxArray, bluecubeClrArray

    def redcube(self, z):
        redcubeVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, z],
                 [1.0, 0.0, z],
                 [1.0, 1.0, z],
                 [0.0, 1.0, z]])
        redcubeIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        redcubeClrArray = np.zeros((8,3))
        redcubeClrArray[:,0] = 1.0
        return redcubeVtxArray, redcubeIdxArray, redcubeClrArray

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

    # If one clicks on the openGL panel it gets the focus
    def mousePressEvent(self, bla):
        self.setFocus()
    #    self.mousePressed = 1
        self.updateGL()

    def spin(self):
        self.yRotDeg = (self.yRotDeg  + 1) % 360.0
        self.updateGL()


class BeamEditor(QWidget):
    '''
    Widget for editing the input beam
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)
        #self.setMinimumSize(500, 500)

        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        self.textBeamEditor = QLabel("Beam Editor")
        grid.addWidget(self.textBeamEditor, 0, 0)

        loadBeamdataButton = QPushButton("Load Beamdata")
        loadBeamdataButton.clicked.connect(self.loadBeamdata)
        grid.addWidget(loadBeamdataButton,1,0)

        loadTwissButton = QPushButton("Load Twiss")
        loadTwissButton.clicked.connect(self.loadTwissWithFname)
        grid.addWidget(loadTwissButton,1,1)

        loadMultipartButton = QPushButton("Load Multiparticles")
        loadMultipartButton.clicked.connect(self.loadMultipartWithFname)
        grid.addWidget(loadMultipartButton,1,2)

    def loadTwissWithFname(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Twiss file', '')
        try:
            twiss = loadTwiss(fname[0])
            # send to parent!!!!
        except:
            print "Bad twiss file!"

    def loadMultipartWithFname(self):
        fname = QFileDialog.getOpenFileName(self, 'Open Multipart file', '')
        try:
            multipart = loadMultipart(fname[0])
            # send to parent!!!!
        except:
            print "Bad multipart file!"

    # To be implemented, takes a beamdata files and gets the data that a beamdata array should have. Implement this in IOHandler
    def loadBeamdata(self):
        return

class LatticeEditor(QWidget):
    '''
    Widget for editing the lattice
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)
        #self.setMinimumSize(500, 500)

        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        self.textLatticeEditor = QLabel("Lattice Editor")
        grid.addWidget(self.textLatticeEditor, 0, 0)

        #loadAction = QAction(QtGui.QIcon('icons/open.png'), 'Load', self)
        #loadAction.setShortcut('Ctrl+O')
        #loadAction.triggered.connect(self.parent.parent.openFile)

        loadButton = QPushButton("Load Lattice")
        loadButton.clicked.connect(self.parent.parent.openFile)
        grid.addWidget(loadButton,1,0)

        #grid.addWidget(loadAction,0,0)
        #self.toolbar = self.addToolBar('Load')
        #self.toolbar.addAction(loadAction)

        #saveAction = QAction(QtGui.QIcon('icons/save.png'), 'Save', self)
        #saveAction.setShortcut('Ctrl+S')
        #saveAction.triggered.connect(self.parent.parent.saveFile)

        saveButton = QPushButton("Save Lattice")
        saveButton.clicked.connect(self.parent.parent.saveFile)
        grid.addWidget(saveButton,1,1)

        #grid.addWidget(saveAction,0,1)
        #self.toolbar = self.addToolBar('Save')
        #self.toolbar.addAction(saveAction)

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
        

        self.textL = QLabel("L:")
        grid.addWidget(self.textL, 3, 0)

        
        self.enterL = QLineEdit()
        grid.addWidget(self.enterL, 3, 1)
        # make sure that only floats and number can be written. Pass the entered value to createDrift
        valueOfL = self.enterL.text()
        #print str(valueOfL)

        self.textK = QLabel("K:")
        grid.addWidget(self.textK, 4, 0)
        self.textK.hide()

        self.enterK = QLineEdit()
        grid.addWidget(self.enterK, 4, 1)
        self.enterK.hide()
        # make sure that only floats and number can be written. Pass the entered value to createDrift

        createElementButton = QPushButton("Create Element")
        createElementButton.clicked.connect(self.createElement) # here arguments should be passed
        grid.addWidget(createElementButton, 5,1)

    def activatedElementSelector(self, text):
        print text
        self.selectedElement = text
        if text == "Quadrupole":
            self.textK.show()
            self.enterK.show()
        else:
            self.textK.hide()
            self.enterK.hide()

    def createElement(self):
        ## Take the inputs in fields
        valueOfL = self.enterL.text()
        try:
            L = float(valueOfL)
        except:
            print "Not a number! L set to 0.0"
            L = 0.0
        #print str(L)

        if not self.enterK.isHidden():
            valueOfK = self.enterK.text()
            try:
                K = float(valueOfK)
            except:
                print "Not a number! K set to 0.0"
                K = 0.0
            #print "K is " + str(K)

        name = "new element"
        spaceChargeOn = 0
        multipart = 0
        twiss = 0
        beamdata = 0
        nbrOfSplits = 1

        if self.selectedElement == "Drift":
            self.facility.createDrift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        elif self.selectedElement == "Quadrupole":
            self.facility.createQuad(name, K, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
        else:
            return
        

        self.parent.latticeoverview.initializeGL() # update the paint lattice in overview
        self.parent.parent.widget.latticeoverview.s_pressed = 1 # prepare a zoom out
        self.parent.parent.widget.latticeoverview.d_pressed = 1
        self.parent.latticeoverview.mousePressEvent(0) # repaint by setting focus    
    
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
        return
        #self.facility.evaluate() # to be implemented

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
        self.layout = QHBoxLayout(self)


        self.layout.addStretch(1)
        #self.layout.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)

        #grid = QGridLayout()
        #self.setLayout(grid)

        self.parent = parent

        try:
            self.parent.openFile()
        except:
            print "Baaaaaaaad lattice file!"
        self.latticeoverview = LatticeOverviewWidget(facility, self)
        self.layout.addWidget(self.latticeoverview) # need to make this expand when increasing horizontal size as well
        #policy = self.latticeoverview.sizePolicy()
        #policy.setHorizontalPolicy(QSizePolicy.Maximum)
        #self.setHorizontalPolicy()
        #grid.addWidget(self.latticeoverview, 0, 0)
        self.latticeoverview.setFocus() # starts with focus
        self.latticeoverview.mousePressEvent(0) # for some reason I need a parameter in this overloaded function (set it to bla sends 0 as bla)

        #editorgrid = QGridLayout() # How can I make this work? (grid within grid) 
        self.editorlayout = QVBoxLayout(self)

        self.beameditor = BeamEditor(self, facility)
        #self.layout.addWidget(self.beameditor)
        self.editorlayout.addWidget(self.beameditor)
        #grid.addWidget(self.beameditor, 0, 1)

        self.latticeeditor = LatticeEditor(self, facility)
        #self.layout.addWidget(self.latticeeditor)
        self.editorlayout.addWidget(self.latticeeditor)
        #grid.addWidget(self.latticeeditor, 1, 1)

        self.evalwidget = EvalWidget(self, facility)
        #self.layout.addWidget(self.evalwidget)
        self.editorlayout.addWidget(self.evalwidget)
        #grid.addWidget(self.evalwidget, 2, 1)

        #grid.addWidget(editorgrid,0,1)
        self.layout.addLayout(self.editorlayout)

        self.setLayout(self.layout)

    

        

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)
        self.facility = Facility()
        self.widget = FormWidget(self, self.facility)
        self.setCentralWidget(self.widget)

        exitAction = QAction(QtGui.QIcon('icons/delete.png'), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(qApp.quit)

        #loadAction = QAction(QtGui.QIcon('icons/open.png'), 'Load', self)
        #loadAction.setShortcut('Ctrl+O')
        #loadAction.triggered.connect(self.openFile)
#
        #saveAction = QAction(QtGui.QIcon('icons/save.png'), 'Save', self)
        #saveAction.setShortcut('Ctrl+S')
        #saveAction.triggered.connect(self.saveFile)
        
        self.toolbar = self.addToolBar('Exit')
        self.toolbar.addAction(exitAction)
        #self.toolbar = self.addToolBar('Load')
        #self.toolbar.addAction(loadAction)
        #self.toolbar = self.addToolBar('Save')
        #self.toolbar.addAction(saveAction)



    def keyPressEvent(self, e):
        
        if e.key() == Qt.Key_Escape:
            self.widget.latticeoverview.spin()
        #if e.key() == "t"
            #print "t was pressed"
        #print str(e.key()) + " was pressed, which is the " + chr(e.key()) + " key"
        #if e.key() == 84:
            #print "t was pressed"
        #keyNbr = e.key()
        #lettera = "a"
        #letterakey = lettera.key()

        # keys: w = 87, a = 65, s = 83, d = 68
        if e.key() not in range(256):
            return
        if chr(e.key()) == 'W':
            self.widget.latticeoverview.w_pressed = 1
            #print "W pressed!!!!!"
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

    def openFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '')
        try:
            spaceChargeOn = 0

            datafilepart = "../data/" + "inpart1000" + ".txt"
            datafiletwiss = "../data/" + "intwiss" + ".txt"

            multipart, twiss = loadSummer2015Formatzasx(datafilepart, datafiletwiss)
            beamdata = getBeamdata()
            nbrOfSplits = 1

            loadedLatticeString = loadLatticeString(fname[0])
            loadedLattice = parseLatticeString(loadedLatticeString, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)
            print "lattice loaded"
            self.facility.setLattice(loadedLattice)
            self.widget.latticeoverview.initializeGL()
            self.widget.latticeoverview.paintGL()
        except AttributeError:
            fnameasstring = ''

    def saveFile(self):
        fname = QFileDialog.getSaveFileName(self, 'Save file', '')
        try:
            latticeToSave = self.facility.getLattice()
            saveLatticeString(fname[0],latticeToSave) # fname[0] is the path string
        except AttributeError:
            fnameasstring = ''

    


if __name__ == '__main__':
    app = QApplication(['DAT Widget Interface'])
    window = DATWidgetInterface()
    window.show()
    app.exec_()