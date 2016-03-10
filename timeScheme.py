# -*-coding:Latin-1 -*

from spatialDiscretisation1d import *


def eulerExplicit(Un, Minv, mesh):
    residual, lambdaMax = computeResidual(Un, mesh)
    Dt = cfl * mesh.cellSize / lambdaMax

    it = np.dot(Minv, residual)
    Unext = Un + Dt * it

    return Unext, Dt


def RK2(Un, Minv, mesh):
    residual, lambdaMax = computeResidual(Un, mesh)
    Dt = cfl * mesh.cellSize / lambdaMax
    it = np.dot(Minv, residual)
    Utild = Un + 0.5 * Dt * it

    residual, lambdaMax = computeResidual(Utild, mesh)
    it = np.dot(Minv, residual)
    Unext = Un + Dt * it

    return Unext, Dt


def SSP2(Un, Minv, mesh):
    residual, lambdaMax = computeResidual(Un, mesh)
    Dt = cfl * mesh.cellSize / lambdaMax
    it = np.dot(Minv, residual)
    U_1 = Un + Dt * it

    residual, lambdaMax = computeResidual(U_1, mesh)
    it = np.dot(Minv, residual)
    Unext = 0.5 * Un + 0.5 * (U_1 + Dt * it)

    return Unext, Dt


def SSP3(Un, Minv, mesh):
    residual, lambdaMax = computeResidual(Un, mesh)
    Dt = cfl * mesh.cellSize / lambdaMax
    it = np.dot(Minv, residual)
    U_1 = Un + Dt * it

    residual, lambdaMax = computeResidual(U_1, mesh)
    it = np.dot(Minv, residual)
    U_2 = 0.75 * Un + 0.25 * (U_1 + Dt * it)

    residual, lambdaMax = computeResidual(U_2, mesh)
    it = np.dot(Minv, residual)
    Unext = (1. / 3.) * Un + (2. / 3.) * (U_2 + Dt * it)

    return Unext, Dt


# ------------------------------------------------------------------

timeSchemeDict = {
    'eulerEplicit': eulerExplicit,
    'RK2': RK2,
    'SSP1': eulerExplicit,
    'SSP2': SSP2,
    'SSP3': SSP3
}
