# -*-coding:Latin-1 -*

import numpy as np

#---------------------------------------------------------

class UniformMesh:

    def __init__(self, a=0.0, b=1.0, c=10):
        self.xOrigine = a
        self.xFin = b
        self.nCells = c
        self.tailleCell = (b-a) / float(c)

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

#---------------------------------------------------------

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

