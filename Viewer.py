"""
===========================
Matplotlib based visualizer
===========================

This file creates and updates the flow field visualizer
"""

import matplotlib.pyplot as plt
import numpy as np


class FlowContours:

    def __init__(self, initial_field, x, y, title):

        self.x = x
        self.y = y
        self.title = title

        self.fig = plt.figure()
        self.im = plt.contourf(x, y, initial_field, 15)
        self.ax = self.fig.gca()

        self.cb = plt.colorbar(self.im)
        self.format_fig(initial_field)

        plt.show(block=False)
        plt.tick_params(labelleft='off', labelbottom='off')
        plt.xticks(x)
        plt.yticks(y)

    def update(self, field):

        self.ax.cla()
        self.im = self.ax.contourf(self.x, self.y, field, 15)
        self.format_fig(field)

    def format_fig(self, field):

        ticks = np.linspace(field.min(), field.max(), num=14, endpoint=True)
        self.cb.boundaries = ticks
        self.cb.set_ticks(ticks)
        self.cb.draw_all()
        plt.title(self.title)
        plt.grid()


def keep_open():
    plt.show(block=True)
