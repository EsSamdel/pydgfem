# -*-coding:Latin-1 -*

import numpy as np
from data import *
from spatialDiscretisation1d import *

#------------------------------------------------------------------

def eulerExplicit(Un, Minv, Mesh):
    intCell = fluxGradPhiCellIntegrator(Un, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(Un, nvar, Mesh)
    sourceTerm = APhiCellIntegrator(Un, nvar, Mesh)
    Dt = cfl * Mesh.cellSize/lambdaMax

#    print intCell
#    print '---'
#    print intBord
#    print '---'
#    print sourceTerm
#    print '---'
    
    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)
    Unext = Un + Dt*it

    return Unext, Dt

#------------------------------------------------------------------

def RK2(Un, Minv, Mesh):

    intCell = fluxGradPhiCellIntegrator(Un, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(Un, nvar, Mesh)
    sourceTerm = APhiCellIntegrator(Un, nvar, Mesh)
    Dt = cfl * Mesh.cellSize/lambdaMax

    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)
    Utild = Un + 0.5*Dt*it

    intCell = fluxGradPhiCellIntegrator(Utild, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(Utild, nvar, Mesh)
    sourceTerm = APhiCellIntegrator(Utild, nvar, Mesh)
    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)

    Unext = Un + Dt*it

    return Unext, Dt


#------------------------------------------------------------------

def SSP2(Un, Minv, Mesh):

    intCell = fluxGradPhiCellIntegrator(Un, nvar, Mesh)    
    intBord, lambdaMax = fluxPhiFaceIntegrator(Un, nvar, Mesh)
    sourceTerm = APhiCellIntegrator(Un, nvar, Mesh)
    Dt = cfl * Mesh.cellSize/lambdaMax

    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)
    U_1 = Un + Dt*it

    intCell = fluxGradPhiCellIntegrator(U_1, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(U_1, nvar, Mesh)
    sourceTerm = APhiCellIntegrator(U_1, nvar, Mesh)
    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)

    Unext = 0.5*Un + 0.5*( U_1 + Dt*it )

    return Unext, Dt

#------------------------------------------------------------------

def SSP3(Un, Minv, Mesh):

    intCell = fluxGradPhiCellIntegrator(Un, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(Un, nvar, Mesh)
    Dt = cfl * Mesh.cellSize/lambdaMax

    res = intCell - intBord
    it = np.dot(Minv, res)
    U_1 = Un + Dt*it

    intCell = fluxGradPhiCellIntegrator(U_1, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(U_1, nvar, Mesh)
    res = intCell - intBord
    it = np.dot(Minv, res)
    U_2 = 0.75*Un  +  0.25* (U_1 + Dt*it)

    intCell = fluxGradPhiCellIntegrator(U_2, nvar, Mesh)
    intBord, lambdaMax = fluxPhiFaceIntegrator(U_2, nvar, Mesh)
    res = intCell - intBord
    it = np.dot(Minv, res)

    Unext = (1./3.)*Un + (2./3.)*(U_2 + Dt*it)

    return Unext, Dt

#------------------------------------------------------------------

timeSchemeDict = {
  'eulerEplicit': eulerExplicit,
  'RK2': RK2,
  'SSP1': eulerExplicit,
  'SSP2': SSP2,
  'SSP3': SSP3
}
