"""
===========================
Solver for 2D CFD Simulation
===========================

This file uses the control volume method with practice a with the simple scheme to solve a 2d flow field
"""


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
        self.failed_message = 'SOLUTION FAILED'

        self.solve()

    def solve(self):
        print self.solve_message
