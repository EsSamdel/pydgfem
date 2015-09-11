# -*-coding:Latin-1 -*

from data import *
import numpy as np
from math import sin, cos, exp, sqrt

PI = 3.14159265359
mode = 2.

# Le second membre :
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