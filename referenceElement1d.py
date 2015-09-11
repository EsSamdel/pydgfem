# -*-coding:Latin-1 -*

import numpy as np
from math import sqrt

#---------------------------------------------------------

def NddlLocal(degre):
    """ Donne le nombre de ddl local a un element fini
    """
    n = degre + 1
    return n

#---------------------------------------------------------

def FonctionBaseLagrange(x, degre):
    """Calcule les fonctions de base au point x pour un degré donné
    """

    phi = np.zeros(NddlLocal(degre))

    if degre == 0:
        phi[0] = 1.

    elif degre == 1:
        phi[0] = 1.-x
        phi[1] = x

    elif degre == 2:
        phi[0] = 2*(x-1)*(x-0.5)
        phi[1] = 4*x*(1.-x)
        phi[2] = 2*x*(x-0.5)

    elif degre == 3 :
        phi[0] = 9/2*(1/3-x)*(2/3-x)*(1-x)
        phi[1] = 27/2*x*(2/3-x)*(1-x)
        phi[2] = 27/2*x*(x-1/3)*(1-x)
        phi[3] = 9/2*x*(x-1/3)*(x-2/3)

    elif degre == 4:
        phi[0] = 3/32*(1/4-x)*(1/2-x)*(3/4-x)*(1-x)
        phi[1] = 3/128*x*(1/2-x)*(3/4-x)*(1-x)
        phi[2] = 1/64*x*(x-1/4)*(3/4-x)*(1-x)
        phi[3] = 3/128*x*(x-1/4)*(x-1/2)*(1-x)
        phi[4] = 3/32*x*(x-1/4)*(x-1/2)*(x-3/4)

    else:
        print('degre = ', degre)
        print('Les fonctions de base en 1d pour ce degre ne sont pas implementees.')

    return phi

#---------------------------------------------------------

def GradientFonctionBaseLagrange(x,degre):
    """Calcule les gradients des fonctions de base au point x pour un degré donné
    """

    gradPhi = np.zeros(NddlLocal(degre))

    if degre == 0:
        gradPhi[0] = 0.

    elif degre == 1:
        gradPhi[0] = -1.
        gradPhi[1] = 1.

    elif degre == 2:
        gradPhi[0] = 2*(x-1) + 2*(x-0.5)
        gradPhi[1] = -4*x + 4*(1.-x)
        gradPhi[2] = 2*x + 2*(x-0.5)

    elif degre == 3:
        gradPhi[0] = -9./2.*(2./3.-x)*(1.-x)-9./2.*(1./3.-x)*(1.-x)-9./2.*(1./3.-x)*(2./3.-x)
        gradPhi[1] = 27./2.*(2./3.-x)*(1.-x)-27./2.*x*(1.-x)-27./2.*x*(2./3.-x)
        gradPhi[2] = 27./2.*(x-1./3.)*(1.-x)+27./2.*x*(1.-x)-27./2.*x*(x-1./3.)
        gradPhi[3] = 9./2.*(x-1./3.)*(x-2./3.)+9./2.*x*(x-2./3.)+9./2.*x*(x-1./3.)

    elif degre == 4:
        gradPhi[0] = 3./32.*(1./4.-x)*(1./2.-x)*(3./4.-x)*(1.-x)
        gradPhi[1] = 3./128.*x*(1./2.-x)*(3./4.-x)*(1.-x)
        gradPhi[2] = 1./64.*x*(x-1./4.)*(3./4.-x)*(1.-x)
        gradPhi[3] = 3./128.*x*(x-1./4.)*(x-1./2.)*(1.-x)
        gradPhi[4] = 3./32.*x*(x-1./4.)*(x-1./2.)*(x-3./4.)

    else:
        print('degre = ', degre)
        print('Les gradients des fonctions de base en 1d pour ce degre ne sont pas implementees.')

    return gradPhi

#---------------------------------------------------------

def FonctionBaseLegendre(x, degre):
    """Calcule les fonctions de base au point x pour un degré donné
    """

    phi = np.zeros(NddlLocal(degre))

    if degre >= 0:
        phi[0] = 1.

    if degre >= 1:
        phi[1] = (x-0.5)

    if degre >= 2:
        phi[2] = sqrt(5.)*(6.*x*x-6.*x+1.)

    if degre >= 3:
        phi[3] = sqrt(7.)*(20.*x*x*x-30.*x*x+12.*x-1.)

    if degre >= 4:
        print('degre = ', degre)
        print('Les fonctions de base en 1d pour ce degre ne sont pas implementees.')

    return phi

#---------------------------------------------------------

def GradientFonctionBaseLegendre(x, degre):
    """Calcule les gradients des fonctions de base au point x pour un degré donné
    """

    gradphi = np.zeros(NddlLocal(degre))

    if degre >= 0:
        gradphi[0] = 0.

    if degre >= 1:
        gradphi[1] = 1.


    if degre >= 2:
        gradphi[2] = sqrt(5.)*(12.*x-6.)


    if degre >= 3:
        gradphi[3] = sqrt(7.)*(60.*x*x-60.*x+12.)


    if degre >= 4:
        print('degre = ', degre)
        print('Les gradients des fonctions de base en 1d pour ce degre ne sont pas implementees.')

    return gradphi

#---------------------------------------------------------

FonctionBase = FonctionBaseLegendre
GradientFonctionBase = GradientFonctionBaseLegendre
