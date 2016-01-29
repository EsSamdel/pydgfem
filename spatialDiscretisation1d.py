# -*-coding:Latin-1 -*

import numpy as np
import copy
from math import *

from mesh import *
from integration1d import *
from referenceElement1d import *
from boundaryConditions import *

from fluxEuler import *
from fluxSW import *
from fluxWBSW import *


#---------------------------------------------------------

def averageValueOnCell(u, nvar, Mesh):
    """ u contient les valeurs aux degrés de liberté d'une fonction,
        projetée sur un espace EF de degre 'degree'. Cette fonction 
        extrait les valeurs aux noeuds, par moyenne des valeurs a
        gauche et à droite de chaque noeud :
        k elts
        x[qk1, qk2 ...]
    """    
    
    val = []
    for ivar in range(nvar):
        val.append([])
    

    for iElt in Mesh.elements:
        degree = iElt.degree
        nddlLoc = iElt.nddl
        omega, x, nPoints = Quadrature1d(2*degree+3)

        uLoc = u[iElt.localPos[0]*nvar: iElt.localPos[1]*nvar]
        uMoyLoc = np.zeros(nvar)

        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degree)

            # Reconstruction de la valeur de u :
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddlLoc):
                    uVal[-1] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]
            
                uMoyLoc[ivar] += uVal[ivar]
           
        for ivar in range(nvar):
            val[ivar].append(uMoyLoc[ivar] / nPoints)

    return val

#---------------------------------------------------------

def projectionOnFE(function, nvar, Mesh):
    """ Calcul la projection de la fonction 'function' sur l'espace EF
        de degre 'degre' pour un maillage 1d :
        Pj = int(f*phi[i]) = (b-a)*Sum(omega[k] * f(a + (b-a)x[i]) * phi[i])
    """
    nddltotal = Mesh.getNddlTot()
    integrale = np.zeros(nddltotal*nvar)

    for iElt in Mesh.elements:
        k = iElt.id
        a = iElt.coord[0]
        b = iElt.coord[1]
        nddlLoc = iElt.nddl
        degreeLoc = iElt.degree
        
        integralePhi = np.zeros(nddlLoc*nvar)
        omega, x, nPoints = Quadrature1d(2*degreeLoc+3)

        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degreeLoc)
            for ivar in range(nvar):
                for kddl in range(nddlLoc):
                    integralePhi[(nddlLoc)*ivar + kddl] += (b-a)*omega[kint] * function((b-a)*x[kint] + a)[ivar] * Phi[kddl]

        # On le place au bon endroit
        integrale[(nddlLoc*nvar)*k:(nddlLoc*nvar)*(k+1)] = integralePhi

    return integrale

#---------------------------------------------------------

def Masse(nvar, Mesh):
    """ Calcul de la matrice de masse EF de degre 'degree' pour un maillage 1d :
        M = int(phi[i]*phi[j]) = (b-a)*Sum(omaga[k]*phi[i]*phi[j])
        U = {q1, q2}
        Ul = {q11, q12, q21, q22}
    """
    
    nddlTot = Mesh.getNddlTot()
    M = np.zeros((nddlTot*nvar, nddlTot*nvar))

    for iElt in Mesh.elements:
        nddlLoc = iElt.nddl
        degree = iElt.degree
        massLoc = np.zeros((nddlLoc*nvar, nddlLoc*nvar))
        omega, x, nPoints = Quadrature1d(2*degree+3)
    
        for k in range(nPoints):
            Phi = FonctionBase(x[k], degree)
            for ivar in range(nvar):
                for i in range(nddlLoc):
                    for j in range(nddlLoc):
                        massLoc[(ivar*nddlLoc) + i, (ivar*nddlLoc) + j] += omega[k]*Phi[i]*Phi[j]
    
        M[(nddlLoc*nvar)*iElt.id:(nddlLoc*nvar)*(iElt.id+1), (nddlLoc*nvar)*iElt.id:(nddlLoc*nvar)*(iElt.id+1)] = iElt.mesure * massLoc
    
#        print massLoc
#        print "----------"
#        print M

    return M

#---------------------------------------------------------

def fluxGradPhiCellIntegrator(u, nvar, Mesh):
    """ 
        Compute cell integral : int_cell F . grad(Phi)
        matrix
    """

    integral = np.zeros(u.size)
    
    for iElt in Mesh.elements:
        nddlLoc = iElt.nddl
        degree = iElt.degree
        integralLoc = np.zeros(nddlLoc*nvar)
        omega, x, nPoints = Quadrature1d(2*degree+3)

        uLoc = u[iElt.localPos[0]*nvar: iElt.localPos[1]*nvar]

        for kint in range(nPoints):
            gradPhi = GradientFonctionBase(x[kint], degree)
            Phi = FonctionBase(x[kint], degree)

            # Reconstruction de la valeur de u :
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddlLoc):
                    uVal[-1] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]

            # Calcul du flux :
            z = Mesh.getBathy(iElt.id, x[kint])
#            print uVal[0], z, iCell, x[kint]
            fluxLoc = fluxWBSW(uVal, z)

            # Compute integral loc :
            for ivar in range(nvar):
                for iddl in range(nddlLoc):
                    integralLoc[(ivar*nddlLoc) + iddl] += omega[kint]*fluxLoc[ivar]*gradPhi[iddl]

        integral[(nddlLoc*nvar)*iElt.id: (nddlLoc*nvar)*(iElt.id+1)] = integralLoc

    return integral


#---------------------------------------------------------


def fluxPhiFaceIntegrator(u, nvar, Mesh):
    """ 
        Calcul la matrice des flux aux aretes
        pour un maillage 1d :
    """

    integral = np.zeros(u.size)

    for iFace in Mesh.innerFaces:
        # reconstruct left and right state :
        # left :
        eltL = iFace.eltLeft
        degree = eltL.degree
        nddlLoc = eltL.nddl

        uLeft = np.zeros(nvar)
        x = 1.0        
        Phi = FonctionBase(x, degree)
        uLoc = u[eltL.localPos[0]*nvar: eltL.localPos[1]*nvar]
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                uLeft[ivar] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]

        zL = Mesh.getBathy(eltL.id, x)

        # right :
        eltR = iFace.eltRight
        degree = eltR.degree
        nddlLoc = eltR.nddl

        uRight = np.zeros(nvar)
        x = 0.0        
        Phi = FonctionBase(x, degree)
        uLoc = u[eltR.localPos[0]*nvar: eltR.localPos[1]*nvar]
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                uRight[ivar] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]

        zR = Mesh.getBathy(eltR.id, x)

        # Flux :
        fluxLoc, lambdaMax = laxFriedrichs(uLeft, uRight, zL, zR)
                
        # integral :
        fluxR = np.zeros(nddlLoc*nvar)
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                fluxR[ivar*nddlLoc + iddl] = fluxLoc[ivar]*Phi[iddl]
        integral[(nddlLoc*nvar)*eltR.id: (nddlLoc*nvar)*(eltR.id+1)] -= fluxR
        
        x = 1.0        
        Phi = FonctionBase(x, degree)
        fluxL = np.zeros(nddlLoc*nvar)
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                fluxL[ivar*nddlLoc + iddl] = fluxLoc[ivar]*Phi[iddl]
        integral[(nddlLoc*nvar)*eltL.id: (nddlLoc*nvar)*(eltL.id+1)] += fluxL
        
    
    for iFace in Mesh.boundaryFaces:
        orient = 'left'
        x = 0.0
        uIn = np.zeros(nvar)
        uBound = np.zeros(nvar)
        
        if iFace.eltLeft == -1:
            eltIn = iFace.eltRight
            x = 0.0
        elif iFace.eltRight == -1:            
            eltIn = iFace.eltLeft
            orient = 'right'
            x = 1.0
            
        nddlLoc = eltIn.nddl
        Phi = FonctionBase(x, eltIn.degree)
        uLoc = u[eltIn.localPos[0]*nvar: eltIn.localPos[1]*nvar]
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                uIn[ivar] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]

        z = Mesh.getBathy(eltIn.id, x)
                
        if boundaryType[iFace.color] == 'wall' :
            uBound[0] = uIn[0]
            uBound[1] = -uIn[1]
        
        fluxLoc = np.zeros(nvar)

        if orient == 'left':
            fluxLoc, lambdaMax = laxFriedrichs(uBound, uIn, z, z)
        else:                
            fluxLoc, lambdaMax = laxFriedrichs(uIn, uBound, z, z)

        # integral :
        flux = np.zeros(nddlLoc*nvar)
        for ivar in range(nvar):
            for iddl in range(nddlLoc):
                flux[ivar*nddlLoc + iddl] = fluxLoc[ivar]*Phi[iddl]
        
        if orient == 'left':
            integral[(nddlLoc*nvar)*eltIn.id: (nddlLoc*nvar)*(eltIn.id+1)] -= flux
        else:
            integral[(nddlLoc*nvar)*eltIn.id: (nddlLoc*nvar)*(eltIn.id+1)] += flux

    return integral, lambdaMax


#---------------------------------------------------------


def APhiCellIntegrator(u, nvar, Mesh):
    """ 
        Compute cell integral : int_cell A . Phi
        matrix
    """
    
    integral = np.zeros(u.size)
    
    for iElt in Mesh.elements:
        nddlLoc = iElt.nddl
        degree = iElt.degree
        intLoc = np.zeros(nddlLoc*nvar)
        omega, x, nPoints = Quadrature1d(2*degree+3)

        uLoc = u[iElt.localPos[0]*nvar: iElt.localPos[1]*nvar]

        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degree)

            # Reconstruction de la valeur de u :
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddlLoc):
                    uVal[-1] += Phi[iddl]*uLoc[(nddlLoc*ivar)+iddl]    
    
            # Source term :
            slope = Mesh.getSlope(iElt.id)
            A = sourceTerm(uVal, slope)
            
            for ivar in range(nvar):
                for iddl in range(nddlLoc):
                    intLoc[(ivar*nddlLoc) + iddl] += omega[kint]*A[ivar]*Phi[iddl]

        integral[(nddlLoc*nvar)*iElt.id: (nddlLoc*nvar)*(iElt.id+1)] = intLoc

    return integral