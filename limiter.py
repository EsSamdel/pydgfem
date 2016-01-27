# -*- coding: utf-8 -*-


import numpy as np
from math import *

from mesh import *
from integration1d import *
from referenceElement1d import *
from boundaryConditions import *
from spatialDiscretisation1d import *



def assertPositivity():
    pass



#---------------------------------------------------------

def limiter(u, degre, Mesh):
    """ Function limiter for TVD:
    """

    if degre >= 1:
        uLimiter = np.zeros(u.size)
        nddllocal = NddlLocal(degre)

        for k in range(Mesh.nElement):

            uLimiter[NumerotationGlobale(0, k, degre)] = u[NumerotationGlobale(0, k, degre)]

            # Elements voisins :
            voisinGauche = k-1
            voisinDroite = k+1
            if voisinGauche < 0:
                voisinGauche = Mesh.nElement - 1
            if voisinDroite > Mesh.nElement - 1:
                voisinDroite = 0

            uMoy = u[NumerotationGlobale(0, k, degre)]
            uMoyGauche = u[NumerotationGlobale(0, voisinGauche, degre)]
            uMoyDroite = u[NumerotationGlobale(0, voisinDroite, degre)]

            deltaGauche = uMoy - uMoyGauche
            deltaDroite = uMoyDroite - uMoy
            u_1 = u[NumerotationGlobale(1, k, degre)]

            val = minmod(u_1, deltaGauche, deltaDroite)
            uLimiter[NumerotationGlobale(1, k, degre)] = val

    else:
        uLimiter = u

    #print('u ', u)
    #print('limiter : ', uLimiter)

    return uLimiter

#---------------------------------------------------------

def minmod(a, b, c):
    """ Function minmod :
    """

    if a >= 0. and b >= 0. and c >= 0.:
        val = min(a, b, c)
    elif a <= 0. and b <= 0. and c <= 0.:
        val = max(a, b, c)
    else:
        val = 0.

    return val
