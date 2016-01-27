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


def funcBathy(x):
    z = 0.
    if x > 0.5 and x < 0.7:
        z = 1.5
    return z
        
def funcInit(x):
    y = np.zeros(2)

    eta = 2.0
    if x > -1.0 and x < -0.5:
        eta = 3.0
    z = funcBathy(x)
    if eta < z:
        eta = z
    h = eta - z

    u = 0.0    
    
    y[0] = eta
    y[1] = h*u
    return y 
    
xDebut = -1.0
xFin = 2.0
nCells = 50

#------------------------------------------
# :: CREATE MESH ::
#------------------------------------------
Mesh = UniformMeshWithBathy(xDebut, xFin, nCells)
Mesh.initBathy(funcBathy)
x = Mesh.milieuxCellules()

#------------------------------------------
# :: INITIALISATION ::
#------------------------------------------
fInit = ProjectionSurEspaceEF(funcInit, nvar, degre, Mesh)

# Mass matrix :
M = Masse(nvar, degre, Mesh)
Minv = np.linalg.inv(M)

# Initial solution :
U0 = np.linalg.solve(M, fInit)
y0 = ValeursAuxCentres(U0, nvar, degre, Mesh)

Mesh.plotBathy()
#plt.plot(x, y0[0])
#plt.show()

plot = GraphAnim(x, y0)

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
#        time.sleep(1)

    print('temps : ', t)

yFinal = ValeursAuxCentres(Un, nvar, degre, Mesh)
plot.update_line(yFinal)


