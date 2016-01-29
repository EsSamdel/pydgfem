# -*-coding:Latin-1 -*
"""
    :: Information ::
Implementation de la méthode Galerkin discontinue pour la résolution numérique
d'équations aux dérivées partielles.
Integration en temps explicite par les méthodes de Runge Kutta et SSP
"""

# Import :
import os, time
path = os.path.abspath(os.path.curdir)

import numpy as np
from data import *
from mesh import *
from icbc import *
from spatialDiscretisation1d import *
from timeScheme import *
from output import *

xDebut = 0.0
xFin = 1.0
nCells = 200
dx = (xFin - xDebut) / nCells

def funcBathy(x):
    z = 0.0
#    if x > 0.3 and x < 0.6:
#        z += 100*exp(0.1/((x-0.3)*(x-0.6)))
    return z
        
def funcInit(x):
    eta = 1.0
    if x > 0.1 and x < 0.4:
        eta = 1.0 + 100*exp(0.1/((x-0.1)*(x-0.4)))
    if eta <= funcBathy(x) + eps:
        eta = funcBathy(x) + 2.0*eps
    y = np.zeros(2)    
    y[0] = eta
    y[1] = 0.0
    return y 

#------------------------------------------
# :: CREATE MESH ::
#------------------------------------------
Mesh = MeshWithBathy(xDebut, xFin, nCells)
Mesh.initBathy(funcBathy)
x = Mesh.getCellCenter()

#------------------------------------------
# :: INITIALISATION ::
#------------------------------------------
fInit = projectionOnFE(funcInit, nvar, Mesh)

# Mass matrix :
M = Masse(nvar, Mesh)
Minv = np.linalg.inv(M)

# Initial solution :
U0 = np.linalg.solve(M, fInit)
y0 = averageValueOnCell(U0, nvar, Mesh)

Mesh.plotBathy()
plot = GraphAnim(x, y0)

#------------------------------------------
# :: TIME LOOP ::
#------------------------------------------
Un = U0
t = 0.0
compt = 1

while t < tFinal:
    compt += 1

    Unext, Dt = timeSchemeDict[timeFunction](Un, Minv, Mesh)
    Un = Unext

    if t + Dt > tFinal:
        Dt = tFinal - t
    t = t + Dt

    if not compt%2:
        ytmp = averageValueOnCell(Un, nvar, Mesh)
        plot.update_line(ytmp)
#        time.sleep(1)

    print('temps : ', t)



