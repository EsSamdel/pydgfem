# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


class Graph():

    def __init__(self, x):
        self.x = x
        self.fig, self.tab = plt.subplots(2, 2)

    def addValues(self, y):
        self.tab[0, 0].plot(self.x, y[0])
        self.tab[0, 0].set_title('rho')
        self.tab[0, 1].plot(self.x, y[1]/y[0])
        self.tab[0, 1].set_title('u')
        p = 0.4*(y[2]-0.5*y[1]*y[1]/y[0])
        self.tab[1, 0].plot(self.x, p)
        self.tab[1, 0].set_title('P')
        self.tab[1, 1].plot(self.x, y[2]/y[0])
        self.tab[1, 1].set_title('E')

    def show(self):
        plt.show()