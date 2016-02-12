import math

from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5 import QtGui
from PyQt5.QtOpenGL import *
from PyQt5.QtWidgets import QMainWindow, QApplication, QGridLayout, QPushButton, QWidget, QHBoxLayout
from PyQt5.QtCore import Qt

class LatticeOverviewWidget(QGLWidget):
    '''
    Widget for giving an overview of the lattice
    '''

    def __init__(self, parent,lattice):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(500, 500)
        self.lattice = lattice
        #self.lattice.printLattice()

    def paintGL(self):
        '''
        Drawing routine
        '''
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()


        glBegin(GL_TRIANGLES)
        glColor(1.0, 0.0, 0.0)
        glVertex(0, 0, 0)
        glVertex(1, 0, 0)
        glVertex(0, 1, 0)
        glEnd()

        #glEnableClientState(GL_VERTEX_ARRAY)

        #glDrawArrays(GL_TRIANGLES, 0, 3)
        #glFlush()

    def resizeGL(self, w, h):
        '''
        Resize the GL window 
        '''
        
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(40.0, 1.0, 1.0, 30.0)
    
    def initializeGL(self):
        '''
        Initialize GL
        '''
        
        # set viewing projection
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glClearDepth(1.0)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(40.0, 1.0, 1.0, 30.0)

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
        super(FormWidget, self).__init__(parent)
        self.layout = QHBoxLayout(self)

        lattice = 1
        self.latticeoverview = LatticeOverviewWidget(self, lattice)
        self.layout.addWidget(self.latticeoverview)

        self.latticeeditor = LatticeEditor(self)
        self.layout.addWidget(self.latticeeditor)

        self.setLayout(self.layout)

class DATWidgetInterface(QMainWindow):
    ''' Example class for using SpiralWidget'''
    
    def __init__(self):
        QMainWindow.__init__(self)  
        widget = FormWidget(self)
        self.setCentralWidget(widget)

if __name__ == '__main__':
    app = QApplication(['DAT Widget Interface'])
    window = DATWidgetInterface()
    window.show()
    app.exec_()