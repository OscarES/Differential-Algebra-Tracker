import math

import sys
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import GLU
from PyQt5 import QtGui
from PyQt5.QtOpenGL import *
from PyQt5.QtWidgets import QMainWindow, QApplication, QGridLayout, QPushButton, QWidget, QHBoxLayout
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5 import QtCore
from PyQt5 import QtOpenGL
from IOHandler import loadLattice
from accelerator import Lattice
from numpy import array
import numpy as np

class LatticeOverviewWidget(QGLWidget):
    '''
    Widget for giving an overview of the lattice
    '''
    def __init__(self, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.yRotDeg = 100.0
        self.setMinimumSize(500, 500)

        #self.timer = QtCore.QTimer(self)
        #self.timer.setInterval(20)
        #self.trigger = pyqtSignal()
        #print str(dir(self))
        #connect(self.timer, QtCore.SIGNAL('timeout()'), self.spin)
        #self.timer.start()

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0,  150))
        self.initGeometry()

        glEnable(GL_DEPTH_TEST)

    def resizeGL(self, width, height):
        if height == 0: height = 1

        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glLoadIdentity()
        glTranslate(0.0, 0.0, -50.0)
        glScale(20.0, 20.0, 20.0)
        glRotate(self.yRotDeg, 0.2, 1.0, 0.3)
        glTranslate(-0.5, -0.5, -0.5)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointerf(self.cubeVtxArray)
        glColorPointerf(self.bluecubeClrArray)
        glDrawElementsui(GL_QUADS, self.cubeIdxArray)

    def initGeometry(self):
        self.cubeVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0],
                 [1.0, 0.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [0.0, 1.0, 1.0]])
        self.cubeIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        self.cubeClrArray = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 1.0 ]]

        self.bluecubeVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0],
                 [1.0, 0.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [0.0, 1.0, 1.0]])
        self.bluecubeIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        self.bluecubeClrArray = np.zeros((8,3))
        #self.bluecubeClrArray[:][0:2] = 0.0
        self.bluecubeClrArray[:,2] = 1.0
        print str(self.bluecubeClrArray)
        #self.bluecubeClrArray.fill([0.0, 0.0, 1.0])

    def spin(self):
        self.yRotDeg = (self.yRotDeg  + 100) % 360.0
        #self.parent.statusBar().showMessage('rotation %f' % self.yRotDeg)
        self.updateGL()



class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.yRotDeg = 100.0
        self.setMinimumSize(500, 500)

        #self.timer = QtCore.QTimer(self)
        #self.timer.setInterval(20)
        #self.trigger = pyqtSignal()
        #print str(dir(self))
        #connect(self.timer, QtCore.SIGNAL('timeout()'), self.spin)
        #self.timer.start()

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0,  150))
        self.initGeometry()

        glEnable(GL_DEPTH_TEST)

    def resizeGL(self, width, height):
        if height == 0: height = 1

        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glLoadIdentity()
        glTranslate(0.0, 0.0, -50.0)
        glScale(20.0, 20.0, 20.0)
        glRotate(self.yRotDeg, 0.2, 1.0, 0.3)
        glTranslate(-0.5, -0.5, -0.5)

        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        glVertexPointerf(self.cubeVtxArray)
        glColorPointerf(self.cubeClrArray)
        glDrawElementsui(GL_QUADS, self.cubeIdxArray)

    def initGeometry(self):
        self.cubeVtxArray = array(
                [[0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0],
                 [1.0, 0.0, 1.0],
                 [1.0, 1.0, 1.0],
                 [0.0, 1.0, 1.0]])
        self.cubeIdxArray = [
                0, 1, 2, 3,
                3, 2, 6, 7,
                1, 0, 4, 5,
                2, 1, 5, 6,
                0, 3, 7, 4,
                7, 6, 5, 4 ]
        self.cubeClrArray = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 1.0 ]]

    def spin(self):
        self.yRotDeg = (self.yRotDeg  + 100) % 360.0
        #self.parent.statusBar().showMessage('rotation %f' % self.yRotDeg)
        self.updateGL()




class LatticeEditor(QWidget):
    '''
    Widget for editing the lattice
    '''

    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(500, 500)

# layout manager (aranges the different widgets)
class FormWidget(QWidget):
    def __init__(self, parent):        
        #  Initialize GLUT and process user parameters
        glutInit(sys.argv);

        #  Request double buffered true color window with Z-buffer
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

        #  Enable Z-buffer depth test
        glEnable(GL_DEPTH_TEST)

        super(FormWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)

        self.rotatingcube = GLWidget(self)
        self.layout.addWidget(self.rotatingcube)

        try:
            lattice = loadLattice("data/" + "savedlattice" + ".npy")
        except:
            print "Baaaaaaaad lattice file!"
        self.latticeoverview = LatticeOverviewWidget(self)#, lattice)
        self.layout.addWidget(self.latticeoverview)

        self.latticeeditor = LatticeEditor(self)
        self.layout.addWidget(self.latticeeditor)

        self.setLayout(self.layout)

    

        

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)  
        self.widget = FormWidget(self)
        self.setCentralWidget(self.widget)

    def keyPressEvent(self, e):
        
        if e.key() == Qt.Key_Escape:
            print "hej"
            self.widget.rotatingcube.spin()
            self.widget.latticeoverview.spin()

if __name__ == '__main__':
    app = QApplication(['DAT Widget Interface'])
    window = DATWidgetInterface()
    window.show()
    app.exec_()