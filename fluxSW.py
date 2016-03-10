import numpy as np
from math import sqrt
from data import *

gamma = 1.4
G8 = gamma - 1.0


def fluxSW(u):
    F = np.zeros(2)
    F[0] = u[1]
    F[1] = u[1] * u[1] / u[0] + 0.5 * g * u[0] * u[0]
    return F


def laxFriedrichs(uL, uR):
    """ Compute Lax-Friedrichs flux
    :param uR:
    :param uL:
    """
    cl = sqrt(g * uL[0])
    ul = uL[1] / uL[0]
    eigenL = [ul - cl, ul, ul + cl]

    cr = sqrt(g * uR[0])
    ur = uR[1] / uR[0]
    eigenR = [ur - cr, ur, ur + cr]

    lambdaMax = 0.
    for eigen in eigenL:
        lambdaMax = max(lambdaMax, abs(eigen))
    for eigen in eigenR:
        lambdaMax = max(lambdaMax, abs(eigen))

    C1 = 0.5 * (fluxSW(uL) + fluxSW(uR))
    C2 = 0.5 * lambdaMax * (uL - uR)

    Flux = C1 + C2

    return Flux, lambdaMax


def HLL(uL, uR):
    cL = sqrt(g * uL[0])
    vL = uL[1] / uL[0]
    eigenL = [vL - cL, vL, vL + cL]

    cR = sqrt(g * uR[0])
    vR = uR[1] / uR[0]
    eigenR = [vR - cR, vR, vR + cR]

    lambdaMax = 0.
    for eigen in eigenL:
        lambdaMax = max(lambdaMax, abs(eigen))
    for eigen in eigenR:
        lambdaMax = max(lambdaMax, abs(eigen))

    hStar = 0.5 * (uL[0] + uR[0]) - 0.25 * (vR - vL) * (uL[0] + uR[0]) / (cL + cR)

    qL = 1.
    if hStar > uL[0]:
        qL = sqrt(0.5 * (hStar + uL[0]) * hStar)
    qR = 1.
    if hStar > uR[0]:
        qR = sqrt(0.5 * (hStar + uR[0]) * hStar)

    SL = vL - cL * qL
    SR = vR + cR * qR

    if SL >= 0:
        Flux = fluxSW(uL)
    elif SR >= 0:
        Flux = (SR * fluxSW(uL) - SL * fluxSW(uR) + SR * SL * (uR - uL)) / (SR - SL)
    else:
        Flux = fluxSW(uR)

    return Flux, lambdaMax


def HLLC(uL, uR):
    cL = sqrt(g * uL[0])
    vL = uL[1] / uL[0]
    eigenL = [vL - cL, vL, vL + cL]

    cR = sqrt(g * uR[0])
    vR = uR[1] / uR[0]
    eigenR = [vR - cR, vR, vR + cR]

    lambdaMax = 0.
    for eigen in eigenL:
        lambdaMax = max(lambdaMax, abs(eigen))
    for eigen in eigenR:
        lambdaMax = max(lambdaMax, abs(eigen))

    hStar = 0.5 * (uL[0] + uR[0]) - 0.25 * (vR - vL) * (uL[0] + uR[0]) / (cL + cR)

    qL = 1.
    if hStar > uL[0]:
        qL = sqrt(0.5 * (hStar + uL[0]) * hStar)
    qR = 1.
    if hStar > uR[0]:
        qR = sqrt(0.5 * (hStar + uR[0]) * hStar)

    SL = vL - cL * qL
    SR = vR + cR * qR
    Sstar = (SL * uR[0] * (vR - SR) - SR * uL[0] * (vL - SL)) / (uR[0] * (vR - SR) - uL[0] * (vL - SL))

    uLstar = np.zeros(nvar)
    coefL = uL[0] * ((SL - vL) / (SL - Sstar))
    uLstar[0] = coefL
    uLstar[1] = coefL * Sstar

    uRstar = np.zeros(nvar)
    coefR = uR[0] * ((SR - vR) / (SR - Sstar))
    uRstar[0] = coefR
    uRstar[1] = coefR * Sstar

    if SL >= 0:
        Flux = fluxSW(uL)
    elif SL <= 0 <= Sstar:
        Flux = fluxSW(uL) + SL * (uLstar - uL)
    elif Sstar <= 0 <= SR:
        Flux = fluxSW(uR) + SR * (uRstar - uR)
    else:
        Flux = fluxSW(uR)

    return Flux, lambdaMax
