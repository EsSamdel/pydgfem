# -*-coding:Latin-1 -*
"""
    :: Information ::
Implementation de la méthode Galerkin discontinue pour la résolution numérique
d'équations aux dérivées partielles.
Integration en temps explicite par les méthodes de Runge Kutta et SSP
"""

# Import :
import os

from mesh import *
from output import *
from timeScheme import *

path = os.path.abspath(os.path.curdir)

# ------------------------------------------
# :: CREATE MESH ::
# ------------------------------------------
mesh = MeshWithBathy(case.x1, case.x2, nCells)
mesh.initBathy(case.bathymetry)
x = mesh.getCellCenter()

# ------------------------------------------
# :: INITIALISATION ::
# ------------------------------------------
fInit = projectionOnFE(case.initialState, nvar, mesh)

# Mass matrix :
M = Masse(nvar, mesh)
Minv = np.linalg.inv(M)

# Initial solution :
U0 = np.linalg.solve(M, fInit)
y0 = averageValueOnCell(U0, nvar, mesh)

mesh.plotBathy()
plot = GraphAnim(x, y0)
# time.sleep(5)

# ------------------------------------------
# :: TIME LOOP ::
# ------------------------------------------
Un = U0
t = 0.0
compt = 1

while t < tFinal:
    compt += 1

    Unext, Dt = timeSchemeDict[timeFunction](Un, Minv, mesh)
    Un = Unext

    if t + Dt > tFinal:
        Dt = tFinal - t
    t = t + Dt

    if not compt % 2:
        ytmp = averageValueOnCell(Un, nvar, mesh)
        plot.update_line(ytmp)
    # time.sleep(1)

    print('temps : ', t)
