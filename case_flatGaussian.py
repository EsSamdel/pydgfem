# -*- coding: utf-8 -*-

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
xFin = 20.0
nCells = 50

def funcBathy(x):
    z = 0.
    return z
        
def funcInit(x):
    eta = 2.4*(1.0 + exp(0.25 * ((x-10.0)*(x-10.0) + 100)))
    
    y = np.zeros(2)
    y[0] = eta
    y[1] = 0
    return y 
    
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


