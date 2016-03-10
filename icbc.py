# -*-coding:Latin-1 -*

from math import exp

import numpy as np

from data import *

PI = 3.14159265359
mode = 2.


class ICBC:
    """
        Initial conditions and boundary conditions
    """

    def factory(caseName):
        if caseName == 'test':
            return Case_Test()
        if caseName == 'test2':
            return Case_Test2()
        if caseName == 'pente':
            return Case_Pente()
        if caseName == 'pente2':
            return Case_Pente2()
        if caseName == 'bump':
            return Case_bump()
        if caseName == 'choc':
            return Case_choc()
        if caseName == 'island':
            return Case_island()
        assert 0, " Test name unknown : " + caseName

    factory = staticmethod(factory)

    def __init__(self):
        self.x1 = 0.0
        self.x2 = 1.0

    def bathymetry(self, x):
        return 0.0

    def initialState(self, x):
        y = np.zeros(nvar)
        return y

    def assertLevelPositivity(self, val, x):
        if val <= self.bathymetry(x) + eps:
            return self.bathymetry(x) + 2.0 * eps
        else:
            return val

    def inflow(self):
        y = np.zeros(nvar)
        return y

    def outflow(self):
        y = np.zeros(nvar)
        return y

    def boundaryType(self, color):
        return 'wall'


class Case_Test(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 1.0

    def bathymetry(self, x):
        return 0.0

    def initialState(self, x):
        eta = 1.0
        if 0.1 < x < 0.4:
            eta = 1.0 + 100 * exp(0.1 / ((x - 0.1) * (x - 0.4)))

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = 0.0
        return y


class Case_Test2(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 1.0

    def bathymetry(self, x):
        return 0.0

    def initialState(self, x):
        eta = x
        #        if x > 0.1 and x < 0.4:
        #            eta = 1.0 + 100*exp(0.1/((x-0.1)*(x-0.4)))

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = 0.0
        return y


class Case_Pente(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 10.0

    def bathymetry(self, x):
        z = 0
        if x > 5.0:
            z = -1.0 + x / 5.0
        return z

    def initialState(self, x):
        eta = 0.8
        #        if x > 1.0 and x < 3.0:
        #            eta += 1*exp(2.0/((x-1.0)*(x-3.0)))

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = 0.0
        return y


class Case_Pente2(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 10.0

    def bathymetry(self, x):
        z = x / 10.0
        return z

    def initialState(self, x):
        eta = 1.05
        if 1.0 < x < 3.0:
            eta += 1 * exp(2.0 / ((x - 1.0) * (x - 3.0)))

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = 0.0
        return y


class Case_bump(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 10.0
        self.v0 = 1.0
        self.eta0 = 2.0

    def bathymetry(self, x):
        z = 0.
        if 3.0 < x < 7.0:
            z += 2.0 * exp(4.0 / ((x - 3.0) * (x - 7.0)))
        return z

    def initialState(self, x):
        eta = self.eta0
        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = eta * self.v0
        return y

    def inflow(self):
        y = np.zeros(nvar)
        y[0] = self.eta0
        y[1] = self.eta0 * self.v0
        return y

    def outflow(self):
        y = np.zeros(nvar)
        y[0] = self.eta0
        y[1] = self.eta0 * self.v0
        return y

    def boundaryType(self, color):
        if color == 1:
            return 'dirichlet'
        else:
            return 'dirichlet'


class Case_choc(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 10.0

    def bathymetry(self, x):
        z = 0.
        return z

    def initialState(self, x):
        eta = 2.0
        u = 1.0
        if x > 5.0:
            eta = 1.0
            u = 2.0

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = eta * u
        return y

    def inflow(self):
        y = np.zeros(nvar)
        y[0] = 2.0
        y[1] = 2.0
        return y

    def outflow(self):
        y = np.zeros(nvar)
        y[0] = 1.0
        y[1] = 2.0
        return y

    def boundaryType(self, color):
        if color == 1:
            return 'dirichlet'
        else:
            return 'dirichlet'


class Case_island(ICBC):
    def __init__(self):
        ICBC.__init__(self)
        self.x1 = 0.0
        self.x2 = 10.0

    def bathymetry(self, x):
        z = 0.
        if 3.0 < x < 7.0:
            z += 1 * exp(4.0 / ((x - 3.0) * (x - 7.0)))
        return z

    def initialState(self, x):
        eta = 0.35
        if 2.0 < x < 3.0:
            eta = 0.5

        eta = self.assertLevelPositivity(eta, x)

        y = np.zeros(2)
        y[0] = eta
        y[1] = 0.0
        return y


# -----------------------------------------------------------------------------------
case = ICBC.factory(testName)
dx = (case.x2 - case.x1) / nCells
