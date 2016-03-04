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
        self.lattice = self.facility.getLattice()
        lattticeString = self.lattice.printLattice()
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
    def mousePressEvent(self):
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

        saveBeamdataButton = QPushButton("Save Beamdata")
        saveBeamdataButton.clicked.connect(self.saveBeamdata)
        grid.addWidget(saveBeamdataButton,1,0)

        saveTwissButton = QPushButton("Save Twiss")
        saveTwissButton.clicked.connect(self.saveTwiss)
        grid.addWidget(saveTwissButton,1,1)

        saveMultipartButton = QPushButton("Save Multiparticles")
        saveMultipartButton.clicked.connect(self.saveMultipart)
        grid.addWidget(saveMultipartButton,1,2)

        loadBeamdataButton = QPushButton("Load Beamdata")
        loadBeamdataButton.clicked.connect(self.loadBeamdata)
        grid.addWidget(loadBeamdataButton,2,0)

        loadTwissButton = QPushButton("Load Twiss")
        loadTwissButton.clicked.connect(self.loadTwiss)
        grid.addWidget(loadTwissButton,2,1)

        loadMultipartButton = QPushButton("Load Multiparticles")
        loadMultipartButton.clicked.connect(self.loadMultipart)
        grid.addWidget(loadMultipartButton,2,2)

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

        self.textK = QLabel("K:")
        grid.addWidget(self.textK, 5, 0)
        self.textK.hide()

        self.enterK = QLineEdit()
        grid.addWidget(self.enterK, 5, 1)
        self.enterK.hide()

        createElementButton = QPushButton("Create Element")
        createElementButton.clicked.connect(self.createElement) # here arguments should be passed
        grid.addWidget(createElementButton, 6,1)

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
        name = self.enterName.text()

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

        ## Make new elements
        if self.selectedElement == "Drift":
            self.facility.createDrift(name, L)
        elif self.selectedElement == "Quadrupole":
            self.facility.createQuad(name, K, L)
        else:
            return
        

        self.parent.latticeoverview.initializeGL() # update the paint lattice in overview
        self.parent.parent.widget.latticeoverview.s_pressed = 1 # prepare a zoom out
        self.parent.parent.widget.latticeoverview.d_pressed = 1
        self.parent.latticeoverview.mousePressEvent() # repaint by setting focus

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

        self.parent = parent

        ## Layout
        self.layout = QHBoxLayout(self)

        self.layout.setAlignment(QtCore.Qt.AlignLeft)
        self.layout.addStretch(1)

        # Lattice overview
        self.latticeoverview = LatticeOverviewWidget(facility, self)
        self.layout.addWidget(self.latticeoverview)

        self.latticeoverview.setFocus() # starts with focus
        self.latticeoverview.mousePressEvent() # for some reason I need a parameter in this overloaded function (set it to bla sends 0 as bla)

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
        newWidth = max(self.frameGeometry().width()-420,0) # 420 is magic number for getting the correct width
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