from math import sqrt

import numpy as np

gamma = 1.4
G8 = gamma - 1.0


def fluxEuler(u):
    P = (gamma - 1.0) * (u[2] - 0.5 * u[1] * u[1] / u[0])
    F = np.zeros(3)
    F[0] = u[1]
    F[1] = u[1] * u[1] / u[0] + P
    F[2] = u[1] * (u[2] + P) / u[0]
    return F


def laxFriedrichs(uL, uR):
    """ Compute Lax-Friedrichs flux
    :param uR: Right state
    :param uL: Left state
    """
    pl = (gamma - 1.0) * (uL[2] - 0.5 * uL[1] * uL[1] / uL[0])
    al = sqrt(gamma * pl / uL[0])
    ul = uL[1] / uL[0]
    eigenL = [ul - al, ul, ul + al]

    pr = (gamma - 1.0) * (uR[2] - 0.5 * uR[1] * uR[1] / uR[0])
    ar = sqrt(gamma * pr / uR[0])
    ur = uR[1] / uR[0]
    eigenR = [ur - ar, ur, ur + ar]

    lambdaMax = 0.
    for eigen in eigenL:
        lambdaMax = max(lambdaMax, abs(eigen))
    for eigen in eigenR:
        lambdaMax = max(lambdaMax, abs(eigen))

    C1 = 0.5 * (fluxEuLer(uL) + fluxEuLer(uR))
    C2 = 0.5 * lambdaMax * (uL - uR)

    Flux = C1 + C2

    return Flux, lambdaMax


def roe(uL, uR):
    """ Compute Roe-Pike flux with entropy fix
        Confers : E.F Toro chap 11
    :param uR: Right state
    :param uL: Left state
    """
    pl = (gamma - 1.0) * (uL[2] - 0.5 * uL[1] * uL[1] / uL[0])
    al = sqrt(gamma * pl / uL[0])
    ul = uL[1] / uL[0]
    Hl = (uL[2] + pl) / uL[0]
    eigenLeft = [ul - al, ul, ul + al]

    pr = (gamma - 1.0) * (uR[2] - 0.5 * uR[1] * uR[1] / uR[0])
    ar = sqrt(gamma * pr / uR[0])
    ur = uR[1] / uR[0]
    Hr = (uR[2] + pr) / uR[0]
    eigenRight = [ur - ar, ur, ur + ar]

    lambdaMax = 0.
    for eigen in eigenLeft:
        lambdaMax = max(lambdaMax, abs(eigen))
    for eigen in eigenRight:
        lambdaMax = max(lambdaMax, abs(eigen))

    # Compute component of average vector :
    d_av = sqrt(uL[0] * uR[0])
    u_av = (sqrt(uL[0]) * ul + sqrt(uR[0]) * ur) / (sqrt(uL[0]) + sqrt(uR[0]))
    H_av = (sqrt(uL[0]) * Hl + sqrt(uR[0]) * Hr) / (sqrt(uL[0]) + sqrt(uR[0]))
    a_av = (G8 * (H_av - 0.5 * u_av * u_av)) ** 0.5

    # Compute average wave strength
    Delta_d = uR[0] - uL[0]
    Delta_u = ur - ul
    Delta_p = pr - pl

    alpha = np.zeros(4)
    alpha[0] = (0.5 * (1 / a_av ** 2) * (Delta_p - d_av * a_av * Delta_u))
    alpha[1] = (Delta_d - Delta_p / a_av ** 2)
    alpha[2] = (0.5 * (1 / a_av ** 2) * (Delta_p + d_av * a_av * Delta_u))
    alpha[3] = d_av

    # Compute wave speed :
    ws = np.zeros(4)
    ws[0] = (abs(u_av - a_av))
    ws[1] = (abs(u_av))
    ws[2] = (abs(u_av + a_av))
    ws[3] = (abs(u_av))

    # Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
    dws = 1.
    if ws[0] < dws:
        ws[0] = 0.5 * (ws[0] * ws[0] / dws + dws)
    if ws[2] < dws:
        ws[2] = 0.5 * (ws[2] * ws[2] / dws + dws)

    # EigenVector :
    R = [np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)]
    R[0][0] = 1.
    R[0][1] = u_av - a_av
    R[0][2] = H_av - a_av * u_av

    R[1][0] = 1.0
    R[1][1] = u_av
    R[1][2] = 0.5 * u_av * u_av

    R[2][0] = 1.0
    R[2][1] = u_av + a_av
    R[2][2] = H_av + u_av * a_av

    # du = ur - ul
    R[3][0] = 0.
    R[3][1] = 0.
    R[3][2] = 0.

    # Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
    diss = np.zeros(3)
    for i in range(4):
        diss += ws[i] * alpha[i] * R[i]

    # Compute Flux :
    Flux = 0.5 * (fluxEuler(uL) + fluxEuler(uR) - diss)

    return Flux, lambdaMax
