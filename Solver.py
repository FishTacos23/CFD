"""
===========================
Solver for 2D CFD Simulation
===========================

This file uses the control volume method with practice a with the simple scheme to solve a 2d flow field
"""
import numpy as np
import math


class Solution:

    def __init__(self, ip, iu, iv, bc, p_view=None, u_view=None, v_view=None):

        """

        This class solves the Flow Field Simulations

        :param ip: initial pressures
        :type ip: ndarray of pressure field
        :param iu: initial x velocities
        :type iu: ndarray of x velocity field
        :param iv: initial y velocities
        :type iv: ndarray of y velocity field
        :param bc: boundary conditions
        :type bc: dict
        :param p_view: viewer class for pressure
        :type p_view: class
        :param u_view: viewer class for x velocity
        :type u_view: class
        :param v_view: viewer class for y velocity
        :type v_view: class
        """

        self.p_view = p_view
        self.u_view = u_view
        self.v_view = v_view

        self.bc = bc

        self.u = iu
        self.v = iv
        self.p = ip

        self.solve_message = 'SOLUTION COMPLETE'

        self.solve()

    def solve(self):

        for i in xrange(10):

            val = np.asarray([np.asarray([j for j in xrange(self.u.shape[1])]) for _ in xrange(self.u.shape[0])])

            self.u += np.sin(val*float(i)*2.*math.pi/(self.u.shape[1]*10.))
            self.v += np.cos(val*float(i)*2.*math.pi/(self.u.shape[1]*10.))
            self.p += self.u+self.v

            self.v_view.update(self.v)
            self.p_view.update(self.p)
            self.u_view.update(self.u)

        print self.solve_message
