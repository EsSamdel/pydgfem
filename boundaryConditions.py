# -*- coding: utf-8 -*-

from data import *
from math import sqrt
import numpy as np

def WBSWEwallBoundary(uIn, z):
    hIn = uIn[0] - z
    vIn = uIn[1]/hIn
    hB = (vIn*vIn - 4*vIn*sqrt(g*hIn) + 4*g*hIn) / (4*g)
    
    uB = np.zeros(2)
    uB[0] = hB + z
    uB[1] = 0.
    
    return uB