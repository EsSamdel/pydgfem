# -*-coding:Latin-1 -*
"""
    :: Information ::
Implementation de la m�thode Galerkin discontinue pour la r�solution num�rique
d'�quations aux d�riv�es partielles.
Integration en temps explicite par les m�thodes de Runge Kutta et SSP
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
from output import *

#------------------------------------------
# :: CREATE MESH ::
#------------------------------------------
Mesh = UniformMeshWithBathy(xDebut, xFin, nCells)
Mesh.initBathy(bathyPente)
Mesh.plotBathy()
x = Mesh.milieuxCellules()

#------------------------------------------
# :: INITIALISATION ::
#------------------------------------------
fInit = ProjectionSurEspaceEF(dambreak, nvar, degre, Mesh)

# Mass matrix :
M = Masse(nvar, degre, Mesh)
Minv = np.linalg.inv(M)

# Initial solution :
U0 = np.linalg.solve(M, fInit)
y0 = ValeursAuxCentres(U0, nvar, degre, Mesh)
plot = GraphAnim(x, y0)
#plot.addValues(y0)
#plot.draw()

#------------------------------------------
# :: TIME LOOP ::
#------------------------------------------
Un = U0
t = 0.0
compt = 1

while t < tFinal:
    compt += 1

    Unext, Dt = timeSchemeDict[timeFunction](Un, Minv, degre, Mesh)
    Un = Unext

    if t + Dt > tFinal:
        Dt = tFinal - t
    t = t + Dt

    if not compt%2:
        ytmp = ValeursAuxCentres(Un, nvar, degre, Mesh)
        plot.update_line(ytmp)

    print('temps : ', t)

yFinal = ValeursAuxCentres(Un, nvar, degre, Mesh)
plot.update_line(yFinal)
