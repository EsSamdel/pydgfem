# -*-coding:Latin-1 -*
"""
    :: Information ::
Implementation de la méthode Galerkin discontinue pour la résolution numérique
d'équations aux dérivées partielles.
Integration en temps explicite par les méthodes de Runge Kutta et SSP
"""

# Import :
import os
path = os.path.abspath(os.path.curdir)

import numpy as np
from data import *
from mesh import *
from icbc import *
from spatialDiscretisation1d import *
from timeScheme import *
from output import Graph

#------------------------------------------
# :: CREATE MESH ::
#------------------------------------------
Mesh = PeriodicMesh(xDebut, xFin, nCells)
x = Mesh.milieuxCellules()
plot = Graph(x)

#------------------------------------------
# :: INITIALISATION ::
#------------------------------------------
fInit = ProjectionSurEspaceEF(initOndeAcoust, nvar, degre, Mesh)

# Mass matrix :
M = Masse(nvar, degre, Mesh)
Minv = np.linalg.inv(M)

# Initial solution :
U0 = np.linalg.solve(M, fInit)
y0 = ValeursAuxCentres(U0, nvar, degre, Mesh)
plot.addValues(y0)

#------------------------------------------
# :: TIME LOOP ::
#------------------------------------------
Un = U0
t = 0.0
compt = 1

while t < tFinal:

    Unext, Dt = timeSchemeDict[timeFunction](Un, Minv, degre, Mesh)
    Un = Unext

    if t + Dt > tFinal:
        Dt = tFinal - t
    t = t + Dt

    print('temps : ', t)

yFinal = ValeursAuxCentres(Un, nvar, degre, Mesh)
plot.addValues(yFinal)
plot.show()
