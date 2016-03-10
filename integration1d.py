# -*-coding:Latin-1 -*

import numpy as np
from math import sqrt


def Quadrature1d(degre):
    """ Renvoie les coefficient de la formule de quadrature
        de Gauss (points et poids) en fonction du degre du
        shéma (nb ddl).
        :param degre:
    """

    if degre == 0 or degre == 1:
        # Gauss, 1 point
        n = 1
        omegai = np.zeros(n)
        xi = np.zeros(n)

        xi[0] = 0.5
        omegai[0] = 1.

    elif degre == 2 or degre == 3:
        # Gauss, 2 points
        n = 2
        omegai = np.zeros(n)
        xi = np.zeros(n)

        omegai[0] = 0.5
        omegai[1] = 0.5
        xi[0] = 0.21132486540518711775
        xi[1] = 0.78867513459481288225

    elif degre == 4 or degre == 5:
        # Gauss, 3 points
        n = 3
        omegai = np.zeros(n)
        xi = np.zeros(n)

        omegai[0] = 0.277777777778
        omegai[1] = 0.444444444444
        omegai[2] = 0.277777777778
        xi[0] = 0.11270166537925831148
        xi[1] = 0.5
        xi[2] = 0.88729833462074168852

    elif degre == 6 or degre == 7:
        # Gauss, 4 points
        n = 4
        omegai = np.zeros(n)
        xi = np.zeros(n)

        xi[0] = 0.069431844203
        xi[1] = 0.330009478208
        xi[2] = 0.669990521792
        xi[3] = 0.930568155797
        omegai[0] = 0.173927422569
        omegai[1] = 0.326072577431
        omegai[2] = 0.326072577431
        omegai[3] = 0.173927422569

    elif degre == 8 or degre == 9:
        # Gauss, 5 points
        n = 5
        omegai = np.zeros(n)
        xi = np.zeros(n)

        xi[0] = 0.0469100770307
        xi[1] = 0.230765344947
        xi[2] = 0.5
        xi[3] = 0.769234655053
        xi[4] = 0.953089922969
        omegai[0] = 0.118463442528
        omegai[1] = 0.23931433525
        omegai[2] = 0.284444444444
        omegai[3] = 0.23931433525
        omegai[4] = 0.118463442528

    elif degre == 10 or degre == 11:
        # Gauss, 6 points
        n = 6
        omegai = np.zeros(n)
        xi = np.zeros(n)

        xi[0] = .03376524289842398608
        xi[1] = .16939530676686774318
        xi[2] = .38069040695840154568
        xi[3] = .61930959304159845432
        xi[4] = .83060469323313225682
        xi[5] = .96623475710157601392

        d = 1. / sqrt(0.0073380204222450993933)
        omegai[0] = d * 0.0073380204222450993933
        omegai[1] = d * 0.015451823343095832149
        omegai[2] = d * 0.020041279329451654676
        omegai[3] = d * 0.020041279329451654676
        omegai[4] = d * 0.015451823343095832149
        omegai[5] = d * 0.0073380204222450993933

    else:
        print('degre = ', degre)
        print('Aucune formule de quadrature 1d, de ce degre n est implementee.')

    return omegai, xi, n
