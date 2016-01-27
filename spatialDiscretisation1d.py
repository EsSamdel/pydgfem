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

def NddlTotal(degre, Mesh):
    """ Nombre total de ddl pour un EF de degre 'degre', sur un maillage 1D
    """
    # TODO a  metre dans mesh
    nddl = (degre + 1) * Mesh.nCells     # Discontinu
    return nddl

#---------------------------------------------------------

def ValeursAuxNoeuds(x, nvar, degre, Mesh):
    """ x contient les valeurs aux degrés de liberté d'une fonction,
        projetée sur un espace EF de degre 'degre'. Cette fonction
        en extrait les valeurs aux noeuds, par moyenne des valeurs a
        gauche et �  droite de chaque noeud :
        k elts
        x[qk1, qk2 ...]
    """
    nddllocal = NddlLocal(degre)

    y = []
    for ivar in range(nvar):
        y.append(np.zeros(Mesh.nSommets))
        y[ivar][0] = x[ivar]

        for i in range(Mesh.nSommets-2):
            y[ivar][i+1] = 0.5 * (x[NumerotationGlobale(nddllocal-1, i, degre)+ivar] + x[NumerotationGlobale(0, i+1, degre)+ivar])

        y[ivar][Mesh.nSommets - 1] = x[NumerotationGlobale(nddllocal-1, Mesh.nCells-1, degre)+ivar]

    return y

#---------------------------------------------------------

def ValeursAuxCentres(x, nvar, degre, Mesh):
    """ x contient les valeurs aux degrés de liberté d'une fonction,
        projetée sur un espace EF de degre 'degre'. Cette fonction
        en extrait la moyenne dans la maille :
        k elts
        x[qk1, qk2 ...]
    """
    nddllocal = NddlLocal(degre)

    y = []
    for ivar in range(nvar):
        y.append(np.zeros(Mesh.nCells))

        for icell in range(Mesh.nCells):
            # Extraction des valeurs locales a l'element
            uloc = x[(nddllocal*nvar)*icell: (nddllocal*nvar)*(icell+1)]

            valeurcentre = 0.
            omega, xi, n = Quadrature1d(degre+1)

            for iquad in range(n):
                xchap = xi[iquad]
                omegai = omega[iquad]
                phi = FonctionBase(xchap, degre)

                for iddl in range(nddllocal):
                    valeurcentre += omegai * uloc[(nddllocal*ivar) + iddl] * phi[iddl]

            y[-1][icell] = valeurcentre

    return y

#---------------------------------------------------------

def NumerotationGlobale(iddl, icell, degre):
    """ Numero global du ddl 'iddl' de la cellule 'icell',
        pour un EF de degre 'degre'
    """
    iddlGlobal = icell * (degre + 1) + iddl
    return iddlGlobal

#---------------------------------------------------------

def GetUlocElement(u, icell, nvar, degre):
    nddllocal = NddlLocal(degre)
    uloc = u[(nddllocal*nvar)*icell: (nddllocal*nvar)*(icell+1)]
    return uloc

#---------------------------------------------------------

def ProjectionSurEspaceEF(function, nvar, degre, Mesh):
    """ Calcul la projection de la fonction 'function' l'espace EF
        de degre 'degre' pour un maillage 1d :
        Pj = int(f*phi[i]) = (b-a)*Sum(omega[k] * f(a + (b-a)x[i]) * phi[i])
    """
    nddltotal = NddlTotal(degre, Mesh)
    integrale = np.zeros(nddltotal*nvar)

    for k in range(Mesh.nElement):
        a = Mesh.Noeuds[Mesh.Elements[k,0]]
        b = Mesh.Noeuds[Mesh.Elements[k,1]]

        # Integrale de (f,phi) pour tous les phi, ddl dans une cellule
        nddllocal = NddlLocal(degre)
        integralePhi = np.zeros(nddllocal*nvar)
        omega, x, nPoints = Quadrature1d(2*degre+3)

        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degre)
            for ivar in range(nvar):
                for kddl in range(nddllocal):
                    integralePhi[(nddllocal)*ivar + kddl] += (b-a)*omega[kint] * function((b-a)*x[kint] + a)[ivar] * Phi[kddl]

        # On le place au bon endroit
        integrale[(nddllocal*nvar)*k:(nddllocal*nvar)*(k+1)] = integralePhi

    return integrale

#---------------------------------------------------------

def Masse(nvar, degre, Mesh):
    """ Calcul de la matrice de masse EF d'degre 'degre'
        pour un maillage 1d :
        M = int(phi[i]*phi[j]) = (b-a)*Sum(omaga[k]*phi[i]*phi[j])
        U = {q1, q2}
        Ul = {q11, q12, q21, q22}
    """
    # Definition de la matrice de masse sur l'element de reference
    nddllocal = NddlLocal(degre)
    nddltotal = NddlTotal(degre, Mesh)

    MassLoc = np.zeros((nddllocal*nvar, nddllocal*nvar))
    omega, x, nPoints = Quadrature1d(2*degre + 1)

    for k in range(nPoints):
        Phi = FonctionBase(x[k], degre)
        for ivar in range(nvar):
            for i in range(nddllocal):
                for j in range(nddllocal):
                    MassLoc[(ivar*nddllocal) + i, (ivar*nddllocal) + j] += omega[k]*Phi[i]*Phi[j]

    # Calcul de la matrice de masse globale
    M = np.zeros((nddltotal*nvar, nddltotal*nvar))
    for k in range(Mesh.nElement):
        mes = Mesh.largeurElement(k)
        M[(nddllocal*nvar)*k:(nddllocal*nvar)*(k+1), (nddllocal*nvar)*k:(nddllocal*nvar)*(k+1)] = mes * MassLoc

#    print MassLoc
#    print "----------"
#    print M

    return M

#---------------------------------------------------------

def ContributionInterne(u, nvar, degre, Mesh):
    """ Calcul de la contribution interne aux cellules
    """
    nddllocal = NddlLocal(degre)
    contribInterne = np.zeros(u.size)

    for iCell in range(Mesh.nElement):
        AdvLoc = np.zeros(nddllocal*nvar)
        omega, x, nPoints = Quadrature1d(2*degre+3)

        uLoc = GetUlocElement(u, iCell, nvar, degre)

        for kint in range(nPoints):
            gradPhi = GradientFonctionBase(x[kint], degre)
            Phi = FonctionBase(x[kint], degre)

            # Reconstruction de la valeur de u :
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddllocal):
                    uVal[-1] += Phi[iddl]*uLoc[(nddllocal*ivar)+iddl]

            # Calcul du flux :
            z = Mesh.getBathy(iCell, x[kint])
#            print uVal[0], z, iCell, x[kint]
            fluxLoc = fluxWBSW(uVal, z)

            # Calcul de l'adv local :
            for ivar in range(nvar):
                for iddl in range(nddllocal):
                    AdvLoc[(ivar*nddllocal) + iddl] += omega[kint]*fluxLoc[ivar]*gradPhi[iddl]

        contribInterne[(nddllocal*nvar)*iCell: (nddllocal*nvar)*(iCell+1)] = AdvLoc

    return contribInterne

#---------------------------------------------------------

def FluxDeBord(u, nvar, degre, Mesh):
    """ Calcul la matrice des flux aux aretes
        pour un maillage 1d :
    """
    nddllocal = NddlLocal(degre)
    fluxBord = np.zeros(u.size)

    for k in range(Mesh.nArete):
        
        eltL = Mesh.ElementsDesAretes[k,0]
        leftCell = eltL
        uLeft = np.zeros(nvar)
        leftBound = False
        if eltL == -1: # bord
            eltL = Mesh.ElementsDesAretes[-1,0]
            leftBound = True
            leftCell = 0
        
        x = 1.0
        Phi = FonctionBase(x, degre)
        uLoc = GetUlocElement(u, eltL, nvar, degre)
        for ivar in range(nvar):
            for iddl in range(nddllocal):
                uLeft[ivar] += Phi[iddl]*uLoc[(nddllocal*ivar)+iddl]

        eltR = Mesh.ElementsDesAretes[k,1]
        rightCell = eltR
        uRight = np.zeros(nvar)        
        rightBound = False
        if eltR == -1: # bord
            eltR = Mesh.ElementsDesAretes[0,1]
            rightBound = True
            rightCell = Mesh.ElementsDesAretes[k,0]

        x = 0.0
        Phi = FonctionBase(x, degre)
        uLoc = GetUlocElement(u, eltR, nvar, degre)
        for ivar in range(nvar):
            for iddl in range(nddllocal):
                uRight[ivar] += Phi[iddl]*uLoc[(nddllocal*ivar)+iddl]


        zL = Mesh.getBathy(leftCell, 1.0)
        zR = Mesh.getBathy(rightCell, 0.0)
        
        if leftBound:                        
#            uLeft = WBSWEwallBoundary(uRight, zL)
            uLeft[0] = uRight[0]
            uLeft[1] = -uRight[1]
        if rightBound:
#            uRight = WBSWEwallBoundary(uLeft, zR)
            uRight[0] = uLeft[0]
            uRight[1] = -uLeft[1]

        # Flux numerique
        fluxLoc, lambdaMax = laxFriedrichs(uLeft, uRight, zL, zR)
#        if fluxLoc[0] or fluxLoc[1] != 0:
#        print fluxLoc
#        print k, zL, zR, uLeft, uRight, eltL, leftCell, eltR, rightCell
#        print ''
        
        # Placement de la valeur a gauche :
        if Mesh.ElementsDesAretes[k,0] != -1:
            iCell = Mesh.ElementsDesAretes[k,0]
            x = 1.
            phi = FonctionBase(x, degre)
            fluxL = np.zeros(nddllocal*nvar)
            for ivar in range(nvar):
                for iddl in range(nddllocal):
                    fluxL[ivar*nddllocal + iddl] = fluxLoc[ivar]*phi[iddl]
            fluxBord[(nddllocal*nvar)*iCell: (nddllocal*nvar)*(iCell+1)] += fluxL

        # Placement de la valeur a droite :
        if Mesh.ElementsDesAretes[k,1] != -1:
            iCell = Mesh.ElementsDesAretes[k,1]
            x = 0.
            phi = FonctionBase(x, degre)
            fluxR = np.zeros(nddllocal*nvar)
            for ivar in range(nvar):
                for iddl in range(nddllocal):
                    fluxR[ivar*nddllocal + iddl] = fluxLoc[ivar]*phi[iddl]
            fluxBord[(nddllocal*nvar)*iCell: (nddllocal*nvar)*(iCell+1)] -= fluxR

    return fluxBord, lambdaMax

#---------------------------------------------------------

def TermeSource(u, nvar, degre, Mesh):
    """ Calcul de la contribution interne aux cellules du second membre
    """
    nddllocal = NddlLocal(degre)
    sourceTerme = np.zeros(u.size)

    for iCell in range(Mesh.nElement):
        stLoc = np.zeros(nddllocal*nvar)
        omega, x, nPoints = Quadrature1d(2*degre+3)

        uLoc = GetUlocElement(u, iCell, nvar, degre)

        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degre)

            # Reconstruction de la valeur de u :
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddllocal):
                    uVal[-1] += Phi[iddl]*uLoc[(nddllocal*ivar)+iddl]

            # Calcul du flux :
            slope = Mesh.getSlope(iCell)
            fluxLoc = sourceTerm(uVal, slope)

            # Calcul de l'adv local :
            for ivar in range(nvar):
                for iddl in range(nddllocal):
                    stLoc[(ivar*nddllocal) + iddl] += omega[kint]*fluxLoc[ivar]*Phi[iddl]
                    
        sourceTerme[(nddllocal*nvar)*iCell: (nddllocal*nvar)*(iCell+1)] = stLoc

    return sourceTerme


#-----------------------------------------------------------------------------------

def assertPositivity(u, nvar, degre, Mesh):
    
    for iCell in range(Mesh.nElement):
        omega, x, nPoints = Quadrature1d(2*degre+3)    
        uLoc = GetUlocElement(u, iCell, nvar, degre)

        zL = Mesh.getBathy(iCell, 0.)
        zR = Mesh.getBathy(iCell, 1.)
        zmax = max(zL, zR)
    
        for kint in range(nPoints):
            Phi = FonctionBase(x[kint], degre)    
            uVal = []
            for ivar in range(nvar):
                uVal.append(0.)
                for iddl in range(nddllocal):
                    uVal[-1] += Phi[iddl]*uLoc[(nddllocal*ivar)+iddl]
                    
            if (uVal[0]-zmax) < eps:
                uVal[0] = zmax + eps