"""
===========================
Matplotlib based visualizer
===========================

This file creates and updates the flow field visualizer
"""

import matplotlib.pyplot as plt
import numpy as np


class FlowContours:

    def __init__(self, initial_field, x, y, vol_bounds, title):

        self.x = x
        self.y = y

        self.bounds = vol_bounds

        self.title = title

        self.fig = plt.figure()
        self.im = plt.contourf(x, y, initial_field.transpose(), 15)
        self.ax = self.fig.gca()
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)

        self.cb = plt.colorbar(self.im)

        self.p_a = 1
        self.p_c = 'w'
        self.grid = False

        self.format_fig(initial_field)

        plt.show(block=False)
        plt.tick_params(labelleft='off', labelbottom='off', color='white')
        plt.ylim([vol_bounds[1], vol_bounds[3]])
        plt.xlim([vol_bounds[0], vol_bounds[2]])

    def update(self, field):

        self.ax.cla()
        self.im = self.ax.contourf(self.x, self.y, field.transpose(), 15)
        self.format_fig(field)

    def format_fig(self, field):

        if field.min() == 0. and field.max() == 0:
            ticks = np.linspace(-10, 10, num=14, endpoint=True)
        else:
            ticks = np.linspace(field.min(), field.max(), num=14, endpoint=True)

        self.cb.boundaries = ticks
        self.cb.set_ticks(ticks)
        self.cb.draw_all()

        self.fig.suptitle(self.title)

        if self.grid:

            self.ax.grid()
            self.ax.set_xticks(self.x)
            self.ax.set_yticks(self.y)

        else:
            self.ax.set_ylim([self.bounds[1], self.bounds[3]])
            self.ax.set_xlim([self.bounds[0], self.bounds[2]])


def keep_open():
    plt.show(block=True)
