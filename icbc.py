# -*-coding:Latin-1 -*

from data import *
import numpy as np
from math import sin, cos, exp, sqrt

PI = 3.14159265359
mode = 2.


def initSimpleWave(x):
    y = np.zeros(2)

    h = 1.0

    if (x > 1) and (x < 1.4):
      h = 1.0 + 1*exp(0.1/((x-1.)*(x-1.4)))
    else:
      h = 1.0

    u = sqrt(g*h)


    y[0] = h
    y[1] = h * u
    return y


def dambreak(x):
    y = np.zeros(2)

    h = 1.0

    if (x > 1) and (x < 3):
      h = 2
    else:
      h = 1

    u = 0.0


    y[0] = h
    y[1] = h * u
    return y
    
#------------------------------------------------------

def bathyPente(x):
    return x/2.0
    
def testWellBalanced(x):
    y = np.zeros(2)

    eta = 2.0    
    z = bathyPente(x)
    if eta < z:
        eta = z
    h = eta - z

    u = 0.0    
    
    y[0] = eta
    y[1] = h*u
    return y


#------------------------------------------------------


def initOndeEntro(x):
    y = np.zeros(3)

    P = 1.0
    u = 1.0
    rho = 1.0

    if (x > 0.1) and (x < 0.3):
      rho = 1. + 100*exp(0.1/((x-0.1)*(x-0.3)))
    else:
      rho = 1.

    y[0] = rho
    y[1] = rho * u
    y[2] = P/(0.4) + 0.5*rho* u*u
    return y

def initOndeAcoust(x):
    y = np.zeros(3)

    rhoRef = 1.2046
    uRef = 0.030886
    pRef = 101300.0

    sigma = 0.02
    c = sqrt(1.4*pRef / rhoRef)

    dp = 200 * exp(-(x -0.2)*(x -0.2) / (2.0*sigma*sigma))
    drho = dp / (c*c)
    dv = dp / (rhoRef*c)

    p = pRef + dp
    u = uRef + dv
    rho = rhoRef + drho

    y[0] = rho
    y[1] = rho * u
    y[2] = p/(0.4) + 0.5*rho* u*u

    return y