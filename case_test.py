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
xFin = 10
nCells = 5
dx = (xFin - xDebut) / nCells

def funcBathy(x):
    z = x/5.0
    if x > 5.0:
        z = 2.0 - x/5.0
    return z
        
def funcInit(x):
    y = np.zeros(2)

    eta = 0.9
    if x < 1.0:
        eta += 0.5*(1.0-x)

    z = funcBathy(x)
    if eta < z:
        eta = z + eps
    
    y[0] = eta
    y[1] = 0.0
    return y 

#------------------------------------------
# :: CREATE MESH ::
#------------------------------------------
Mesh = Mesh(xDebut, xFin, nCells)
Mesh.initBathy(funcBathy)
x = Mesh.milieuxCellules()



