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
from IOHandler import loadLattice
from accelerator import Lattice
from numpy import array
import numpy as np
from scipy import *
import time
from facility import *

class LatticeOverviewWidget(QGLWidget):
    '''
    Widget for giving an overview of the lattice
    '''
    def __init__(self, facility, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.yRotDeg = 1.0
        self.setMinimumSize(500, 500)

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


    def loadLattice(self):
        lattticeString = self.lattice.printLattice()
        self.elements = []
        nextWillBeL = 0
        for line in lattticeString.split():
            #print "line: " + line
            if nextWillBeL:
                self.elements.append([tempword, float(line)])
                nextWillBeL = 0
            if line == "drift" or line == "dipole" or line == "quad" or line == "liealgelem":  
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

            length = elem[1]

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
        elemtype = elem[4]
        if elemtype == "quad":
            glRotate(45, 0.0, 0.0, 1.0)
        glTranslate(-0.5, -0.5, -0.5)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointerf(elem[0])
        glColorPointerf(elem[2])
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

    # If one clicks on the openGL panel it gets the focus
    def mousePressEvent(self, bla):
        self.setFocus()
        self.updateGL()






class LatticeEditor(QWidget):
    '''
    Widget for editing the lattice
    '''

    def __init__(self, parent, facility):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(500, 500)

        self.parent = parent

        self.facility = facility

        grid = QGridLayout()
        self.setLayout(grid)

        self.selectedElement = "Drift"

        self.elementSelector = QComboBox(self)
        self.elementSelector.addItem("Drift")
        self.elementSelector.addItem("Dipole")
        self.elementSelector.addItem("Quadrupole")
        self.elementSelector.addItem("Sextupole")
        self.elementSelector.addItem("RF-cavity")
        self.elementSelector.addItem("Higher order element")
        grid.addWidget(self.elementSelector, 0, 0)

        self.elementSelector.activated[str].connect(self.activatedElementSelector)
        self.selectedElement = "Drift"
        

        self.textL = QLabel("L:")
        grid.addWidget(self.textL, 1, 0)

        
        self.enterL = QLineEdit()
        grid.addWidget(self.enterL, 1, 1)
        # make sure that only floats and number can be written. Pass the entered value to createDrift
        valueOfL = self.enterL.text()
        #print str(valueOfL)

        self.textK = QLabel("K:")
        grid.addWidget(self.textK, 2, 0)
        self.textK.hide()

        self.enterK = QLineEdit()
        grid.addWidget(self.enterK, 2, 1)
        self.enterK.hide()
        # make sure that only floats and number can be written. Pass the entered value to createDrift

        createElementButton = QPushButton("Create Element")
        createElementButton.clicked.connect(self.createElement) # here arguments should be passed
        grid.addWidget(createElementButton, 3,1)

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


    #def createDrift(self, name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits):
        #self.facility.createDrift(name, L, spaceChargeOn, multipart, twiss, beamdata, nbrOfSplits)

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

        self.parent = parent

        try:
            lattice = loadLattice("data/" + "savedlattice" + ".npy") # it should be facility.importLattice(........)
        except:
            print "Baaaaaaaad lattice file!"
        self.latticeoverview = LatticeOverviewWidget(facility, self)
        self.layout.addWidget(self.latticeoverview)
        self.latticeoverview.setFocus() # starts with focus
        self.latticeoverview.mousePressEvent(0) # for some reason I need a parameter in this overloaded function (set it to bla sends 0 as bla)

        self.latticeeditor = LatticeEditor(self, facility)
        self.layout.addWidget(self.latticeeditor)

        self.setLayout(self.layout)

    

        

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)
        self.facility = Facility()
        self.widget = FormWidget(self, self.facility)
        self.setCentralWidget(self.widget)

    def keyPressEvent(self, e):
        
        if e.key() == Qt.Key_Escape:
            #print "hej"
            self.widget.rotatingcube.spin()
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


if __name__ == '__main__':
    app = QApplication(['DAT Widget Interface'])
    window = DATWidgetInterface()
    window.show()
    app.exec_()