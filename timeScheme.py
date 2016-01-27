# -*-coding:Latin-1 -*

import numpy as np
from data import *
from spatialDiscretisation1d import *

#------------------------------------------------------------------

def eulerExplicit(Un, Minv, degre, Mesh):
    intCell = ContributionInterne(Un, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(Un, nvar, degre, Mesh)
    sourceTerm = TermeSource(Un, nvar, degre, Mesh)
    Dt = cfl * deltax/lambdaMax

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

def RK2(Un, Minv, degre, Mesh):

    intCell = ContributionInterne(Un, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(Un, nvar, degre, Mesh)
    sourceTerm = TermeSource(Un, nvar, degre, Mesh)
    Dt = cfl * deltax/lambdaMax

    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)
    Utild = Un + 0.5*Dt*it

    intCell = ContributionInterne(Utild, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(Utild, nvar, degre, Mesh)
    sourceTerm = TermeSource(Utild, nvar, degre, Mesh)
    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)

    Unext = Un + Dt*it

    return Unext, Dt


#------------------------------------------------------------------

def SSP2(Un, Minv, degre, Mesh):

    intCell = ContributionInterne(Un, nvar, degre, Mesh)    
    intBord, lambdaMax = FluxDeBord(Un, nvar, degre, Mesh)
    sourceTerm = TermeSource(Un, nvar, degre, Mesh)
    Dt = cfl * deltax/lambdaMax

#    print intCell
#    print '---'
#    print intBord
#    print '---'
#    print sourceTerm
#    print '---'

    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)
    U_1 = Un + Dt*it

    intCell = ContributionInterne(U_1, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(U_1, nvar, degre, Mesh)
    sourceTerm = TermeSource(U_1, nvar, degre, Mesh)
    res = intCell - intBord + sourceTerm
    it = np.dot(Minv, res)

    Unext = 0.5*Un + 0.5*( U_1 + Dt*it )

    return Unext, Dt

#------------------------------------------------------------------

def SSP3(Un, Minv, degre, Mesh):

    intCell = ContributionInterne(Un, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(Un, nvar, degre, Mesh)
    Dt = cfl * deltax/lambdaMax

    res = intCell - intBord
    it = np.dot(Minv, res)
    U_1 = Un + Dt*it

    intCell = ContributionInterne(U_1, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(U_1, nvar, degre, Mesh)
    res = intCell - intBord
    it = np.dot(Minv, res)
    U_2 = 0.75*Un  +  0.25* (U_1 + Dt*it)

    intCell = ContributionInterne(U_2, nvar, degre, Mesh)
    intBord, lambdaMax = FluxDeBord(U_2, nvar, degre, Mesh)
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
