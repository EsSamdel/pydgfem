# -*-coding:Latin-1 -*

from itertools import count
import numpy as np
import matplotlib.pyplot as plt


class Elements:
    _ids = count(0)
    _all = []
    
    def __init__(self, mesure, coord, faceList):
        self.id = self._ids.next()
        Elements._all.append(self)        
        self.mesure = mesure
        self.coord = coord
        self.faces = faceList
        
        self.center = (coord[0] + coord[1]) * 0.5


class Faces:
    _ids = count(0)
    _all = []    
    
    def __init__(self, color):
        self.id = self._ids.next()
        Faces._all.append(self)
        self.color = color
        self.elements = []
        self.nodes = []

    def addElement(self, element):
        self.elements.append(element)
        
        
#---------------------------------------------------------


class Mesh:
    
    def __init__(self, xOri=0.0, xEnd=1.0, nCells=10):
        self.xStart = xOri
        self.xEnd = xEnd
        self.nElements = nCells
        self.cellSize = (xEnd-xOri) / float(nCells)


        # Maillage 1D
        for iFace in range(self.nElements + 1):
            Faces(color=0)
            
        for iElt in range(self.nElements):
            faces = [Faces._all[iElt], Faces._all[iElt+1]]
            coord = [self.xStart + iElt * self.cellSize, self.xStart + (iElt+1) * self.cellSize]
            Elements(self.cellSize, coord, faces)
            for iFace in faces:
                iFace.addElement(Elements._all[-1])

        Faces._all[0].color = 1
        Faces._all[-1].color = 2

#        for iElt in Elements._all:
#            print iElt.center



class UniformMesh:

    def __init__(self, xOri=0.0, xEnd=1.0, nCells=10):
        self.xOrigine = xOri
        self.xFin = xEnd
        self.nCells = nCells
        self.tailleCell = (xEnd-xOri) / float(nCells)

        self.nNodes = self.nCells + 1
        self.nArete = self.nCells + 1
        self.nElement = self.nCells
        self.nAreteBord = 2

        self.Noeuds = np.zeros(self.nNodes)
        for i in range(self.nNodes):
            self.Noeuds[i] = self.xOrigine + i * self.tailleCell

        self.Aretes = np.zeros(self.nArete)
        for i in range(self.nArete):
            self.Aretes[i] = i

        self.Elements = np.zeros((self.nElement, 2))
        for i in range(self.nElement):
            self.Elements[i,0] = i
            self.Elements[i,1] = i+1

        self.ElementsDesAretes = np.zeros((self.nArete, 2))
        self.ElementsDesAretes[0,0] = -1
        self.ElementsDesAretes[0,1] = 0
        for i in range(self.nArete - 2):
            self.ElementsDesAretes[i+1,0] = i
            self.ElementsDesAretes[i+1,1] = i+1

        self.ElementsDesAretes[self.nArete-1,0] = self.nCells-1
        self.ElementsDesAretes[self.nArete-1,1] = -1

        self.AretesBord = np.zeros(2)
        self.AretesBord[0] = 0
        self.AretesBord[1] = self.nCells-1

        self.CouleurAreteBord = np.zeros(self.nAreteBord)
        self.CouleurAreteBord[0] = 100
        self.CouleurAreteBord[1] = 200

    def milieuxCellules(self):
        nCellules = self.Noeuds.size - 1
        x = np.zeros(nCellules)
        for i in range(nCellules):
            x[i] = 0.5 * (self.Noeuds[i] + self.Noeuds[i+1])
        return x

    def largeurElement(self, i):
        try:
            mes = abs( self.Noeuds[self.Elements[i,1]] - self.Noeuds[self.Elements[i,0]] )
            return mes
        except:
            mes = 0.0
            print('Maillage.largeurElement :: Element inconnu !!')
            return mes

#-----------------------------------------------------------------------------------

class PeriodicMesh(UniformMesh):
    """
      Periodic mesh
    """

    def __init__(self, a=0.0, b=1.0, c=10):

        UniformMesh.__init__(self, a, b, c)

        # On retire les aretes de bord
        self.nAreteBord = 0
        self.AretesBord = []
        self.CouleurAreteBord = []

        # Enlever une arete
        self.Aretes = self.Aretes[0:(self.nArete-1)]

        # Les cellules ne changent pas

        # On change ElementsDesAretes
        self.ElementsDesAretes[0,0] = self.nElement-1
        Temp = self.ElementsDesAretes[0:(self.nArete-1), :]
        self.ElementsDesAretes = []
        self.ElementsDesAretes = Temp

        self.nArete = self.nArete-1


#-----------------------------------------------------------------------------------

class UniformMeshWithBathy(UniformMesh):
    
    def __init__(self, xOri=0.0, xEnd=1.0, nCells=10):
        
        UniformMesh.__init__(self, xOri, xEnd, nCells)

        self.bathySize = self.nCells * 2
        self.bathyCoord = np.zeros(self.bathySize)        
        self.bathy = np.zeros(self.bathySize)
        
        
    def initBathy(self, func):
        """
        func : func(x) give ground profile
        """ 
        for i in range(self.nCells):
            self.bathyCoord[2*i] = self.Noeuds[i]
            self.bathyCoord[2*i+1] = self.Noeuds[i+1]
            self.bathy[2*i] = func(self.Noeuds[i])
            self.bathy[2*i + 1] = func(self.Noeuds[i+1])

        
    def plotBathy(self):
        plt.plot(self.bathyCoord, self.bathy)
#        plt.show()
        
    
    def getBathy(self, iCell, xCell):
        slope = self.bathy[2*iCell+1] - self.bathy[2*iCell]
        return self.bathy[2*iCell] + xCell*slope
        
        
    def getSlope(self, iCell):
        return self.bathy[2*iCell+1] - self.bathy[2*iCell]