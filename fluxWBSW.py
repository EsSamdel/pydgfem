# -*- coding: utf-8 -*-

from math import sqrt

import numpy as np

from data import *

gamma = 1.4
G8 = gamma - 1.0


def fluxWBSW(u, z):
    h = u[0] - z
    F = np.zeros(2)
    F[0] = u[1]
    F[1] = u[1] * u[1] / h + 0.5 * g * (u[0] * u[0] - 2 * u[0] * z)
    return F


def sourceTerm(u, slope):
    S = np.zeros(2)
    S[0] = 0
    S[1] = -g * u[0] * slope
    return S


def computeWBState(uL, uR, zL, zR):
    zStar = max(zL, zR)
    zTild = max(0., zStar - uL[0])
    hTildL = max(0., uL[0] - zStar)
    hTildR = max(0., uR[0] - zStar)
    etaTildL = hTildL + zTild
    etaTildR = hTildR + zTild

    wL = np.zeros(2)
    wL[0] = etaTildL
    wL[1] = hTildL / (etaTildL - zL) * uL[1]

    wR = np.zeros(2)
    wR[0] = etaTildR
    wR[1] = hTildR / (etaTildR - zR) * uR[1]

    return wL, wR


def laxFriedrichs(uL, uR, zL, zR):
    """ Compute Lax-Friedrichs flux
    :param uL:
    :param uR:
    :param zL:
    :param zR:
    """
    epsilon = 0.001

    hL = uL[0] - zL
    #    if hL < 0.0:
    #        print 'error left : ', hL, uL[0], zL
    if hL < epsilon:
        uL[0] = zL + epsilon
        uL[1] = 0

    hR = uR[0] - zR
    #    if hR < 0.0:
    #        print 'error right : ', hR, uR[0], zR
    if hR < epsilon:
        uR[0] = zR + epsilon
        uR[1] = 0

    zStar = max(zL, zR)
    zTild = zStar - max(0., zStar - uL[0])
    hTildL = max(0., uL[0] - zStar)
    hTildR = max(0., uR[0] - zStar)
    etaTildL = hTildL + zTild
    etaTildR = hTildR + zTild

    #    print uL[0], zL, uR[0], zR, zStar, zTild, hTildL, etaTildL, hTildR, etaTildR

    eigenL = abs(uL[1] / (uL[0] - zL)) + sqrt(g * (uL[0] - zL))
    eigenR = abs(uR[1] / (uR[0] - zR)) + sqrt(g * (uR[0] - zR))
    lambdaMax = max(eigenL, eigenR)

    wL = np.zeros(2)
    wL[0] = etaTildL
    wL[1] = hTildL / (uL[0] - zL) * uL[1]

    wR = np.zeros(2)
    wR[0] = etaTildR
    wR[1] = hTildR / (uR[0] - zR) * uR[1]

    C1 = 0.5 * (fluxWBSW(wL, zTild) + fluxWBSW(wR, zTild))
    C2 = 0.5 * lambdaMax * (wL - wR)

    Flux = C1 + C2

    Flux[1] += g * etaTildL * (zTild - zL)

    #    if zL != 0.0:
    #        print "State = ", uL[0], uR[0], "Flux = ", Flux

    return Flux, lambdaMax
