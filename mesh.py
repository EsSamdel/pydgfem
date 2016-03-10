# -*-coding:Latin-1 -*

from itertools import count

import matplotlib.pyplot as plt
import numpy as np

from data import *


class Elements:
    """
        Define finite elements (only for 1D)
    """
    _ids = count(0)
    _all = []

    def __init__(self, mesure, coord, degree, faceList):
        self.id = self._ids.next()
        Elements._all.append(self)
        self.mesure = mesure
        self.coord = coord
        self.center = (coord[0] + coord[1]) * 0.5

        self.degree = 0
        self.nddl = 0
        self.setDegree(degree)
        self.localPos = []
        if self.id == 0:
            self.localPos = [0, self.nddl]
        else:
            self.localPos = [Elements._all[self.id - 1].localPos[1], Elements._all[self.id - 1].localPos[1] + self.nddl]

        self.faces = faceList

    def setDegree(self, degree):
        self.degree = degree
        self.nddl = degree + 1

    def getGeomPoint(self, x):
        return self.coord[0] + x * self.mesure


class Faces:
    """
        Define faces of elements (only for 1D)
    """
    _ids = count(0)
    _all = []

    def __init__(self, color):
        self.id = self._ids.next()
        Faces._all.append(self)
        self.color = color
        self.eltLeft = -1
        self.eltRight = -1
        self.nodes = []
        self.elements = []

    def addElement(self, element):
        self.elements.append(element)


# ---------------------------------------------------------


class Mesh:
    """
        Class to handle list of elements, faces, geometrical points, and their conectivity
    """

    def __init__(self, xOri=0.0, xEnd=1.0, nCells=10):
        self.xStart = xOri
        self.xEnd = xEnd
        self.nElements = nCells
        self.cellSize = (xEnd - xOri) / float(nCells)

        # Maillage 1D
        for iFace in range(self.nElements + 1):
            Faces(color=0)

        for iElt in range(self.nElements):
            faces = [Faces._all[iElt], Faces._all[iElt + 1]]
            coord = [self.xStart + iElt * self.cellSize, self.xStart + (iElt + 1) * self.cellSize]
            Elements(self.cellSize, coord, degree, faces)

            faces[0].eltRight = Elements._all[-1]
            faces[1].eltLeft = Elements._all[-1]

        # Definition des bords
        Faces._all[0].color = 1
        Faces._all[-1].color = 2

        self.elements = Elements._all
        self.innerFaces = [iFace for iFace in Faces._all if iFace.color == 0]
        self.boundaryFaces = [iFace for iFace in Faces._all if iFace.color != 0]

    #        for iElt in self.elements:
    #            print iElt.localPos
    #        for iFace in self.innerFaces:
    #            print iFace.id, iFace.eltLeft, iFace.eltRight
    #        for iFace in self.boundaryFaces:
    #            print iFace.id, iFace.eltLeft, iFace.eltRight

    def getNddlTot(self):
        """
            Give the total number of ddls
        """
        nddlTot = 0
        for iElt in Elements._all:
            nddlTot += iElt.nddl
        return nddlTot

    def getNNodesTot(self):
        """
            Give the total number of nodes
        """
        nNodesTot = 0
        for iElt in Elements._all:
            nNodesTot += len(iElt.coord)
        return nNodesTot

    def getCellCenter(self):
        """
            Give an array containing cell center position
        """
        x = np.zeros(self.nElements)
        for iElt in range(self.nElements):
            x[iElt] = Elements._all[iElt].center
        return x


class MeshWithBathy(Mesh):
    """
        Inherit from Mesh. Handle bathymetry data for shallow water problems
    """

    def __init__(self, xOri=0.0, xEnd=1.0, nCells=10):
        Mesh.__init__(self, xOri, xEnd, nCells)

        self.bathySize = 2 * self.nElements
        self.bathyCoord = np.zeros(self.bathySize)
        self.bathy = np.zeros(self.bathySize)

    def initBathy(self, func):
        """
            init bathymetry data given an analytical function func(x)
            :param func:
        """
        # This way the bathy is only order 2
        for i in range(self.nElements):
            self.bathyCoord[2 * i] = Elements._all[i].coord[0]
            self.bathyCoord[2 * i + 1] = Elements._all[i].coord[1]
            self.bathy[2 * i] = func(self.bathyCoord[2 * i])
            self.bathy[2 * i + 1] = func(self.bathyCoord[2 * i + 1])

    def getSlope(self, iCell):
        return self.bathy[2 * iCell + 1] - self.bathy[2 * iCell]

    def getBathy(self, iCell, xCell):
        slope = self.getSlope(iCell)
        return self.bathy[2 * iCell] + xCell * slope

    def plotBathy(self):
        #        print self.bathyCoord
        plt.plot(self.bathyCoord, self.bathy)
