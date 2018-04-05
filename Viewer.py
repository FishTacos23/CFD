"""
===========================
Matplotlib based visualizer
===========================

This file creates and updates the flow field visualizer
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


class FlowContours:

    def __init__(self, initial_field, x, y, vol_bounds, title):

        self.x = x
        self.y = y

        self.bounds = vol_bounds

        self.title = title

        self.fig = plt.figure()
        self.im = plt.contourf(x, y, initial_field, 15)
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

    def update(self, field):

        self.ax.cla()
        self.im = self.ax.contourf(self.x, self.y, field, 15)
        self.format_fig(field)

    def format_fig(self, field):

        if field.min() == 0. and field.max() == 0:
            ticks = np.linspace(-10, 10, num=14, endpoint=True)
        else:
            ticks = np.linspace(field.min(), field.max(), num=14, endpoint=True)

        self.cb.boundaries = ticks
        self.cb.set_ticks(ticks)
        self.cb.draw_all()

        self.ax.add_patch(patches.Rectangle((self.x.min(), self.y.min()), self.bounds[0]-self.x.min(),
                                            self.y.max()-self.y.min(), alpha=self.p_a, color=self.p_c))

        self.ax.add_patch(patches.Rectangle((self.x.max(), self.y.min()), self.bounds[2]-self.x.max(),
                                            self.y.max()-self.y.min(), alpha=self.p_a, color=self.p_c))

        self.ax.add_patch(patches.Rectangle((self.bounds[0], self.y.min()), self.bounds[2]-self.bounds[0],
                                            self.bounds[1]-self.y.min(), alpha=self.p_a, color=self.p_c))

        self.ax.add_patch(patches.Rectangle((self.bounds[0], self.y.max()), self.bounds[2]-self.bounds[0],
                                            self.bounds[3]-self.y.max(), alpha=self.p_a, color=self.p_c))

        self.fig.suptitle(self.title)

        if self.grid:
            self.ax.grid()
            self.ax.set_xticks(self.x)
            self.ax.set_yticks(self.y)


def keep_open():
    plt.show(block=True)
