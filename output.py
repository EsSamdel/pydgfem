# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt


class GraphAnim():
    def __init__(self, x, y):
        plt.ion()
        self.hl, = plt.plot(x, y[0])
#        self.hl.axes.set_ylim([0,3])

    def update_line(self, y):
        self.hl.set_ydata(y[0])
        plt.draw()

class Graph1():

    def __init__(self, x):
        self.x = x
        self.fig, self.tab = plt.subplots(1, 1)

    def addValues(self, y):
        self.tab.plot(self.x, y[0])
        self.tab.set_ylim([0,2])
        self.tab.set_title('h')

    def updateData(self, y):
        self.fig.set_ydata(y[0])
        
    def draw(self):
        plt.draw()

    def show(self):
        plt.show()

    def ion(self):
        plt.ion()
        

class Graph4():

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