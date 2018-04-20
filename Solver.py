"""
===========================
Solver for 2D CFD Simulation
===========================

This file uses the control volume method with practice a with the simple scheme to solve a 2d flow field
"""
import numpy as np
import math
from sys import stdout


class Solution:

    def __init__(self, ip, iu, iv, lx, bx, ly, by, bc, cc, au, av, ap, rho, mu, p_view=None, u_view=None, v_view=None):

        """

        This class solves the Flow Field Simulations

        :param ip: initial pressures
        :type ip: ndarray
        :param iu: initial x velocities
        :type iu: ndarray
        :param iv: initial y velocities
        :type iv: ndarray
        :param lx: i node locations
        :type lx: ndarray
        :param bx: I node locations
        :type bx: ndarray
        :param ly: j node locations
        :type ly: ndarray
        :param by: J node locations
        :type by: ndarray
        :param bc: boundary conditions
        :type bc: dict
        :param cc: convergence criterion
        :type cc: float
        :param au: x velocity relaxation parameter
        :type au: float
        :param av: x velocity relaxation parameter
        :type av: float
        :param ap: pressure relaxation parameter
        :type ap: float
        :param rho: fluid density
        :type rho: float
        :param mu: fluid viscosity
        :type mu: float
        :param p_view: viewer class for pressure
        :type p_view: class
        :param u_view: viewer class for x velocity
        :type u_view: class
        :param v_view: viewer class for y velocity
        :type v_view: class
        """

        # Viewers
        self.p_view = p_view
        self.u_view = u_view
        self.v_view = v_view

        # CV Properties
        self.bc = bc

        self.x = lx
        self.X = bx
        self.Y = by
        self.y = ly

        self.ni = self.x.size
        self.nj = self.y.size
        self.n_max = self.ni*self.nj

        # Fluid Properties
        self.rho = rho
        self.mu = mu

        # Solver Properties
        self.convergence = cc
        self.au = au
        self.av = av
        self.ap = ap

        # Star Values
        self.u_s = iu
        self.v_s = iv
        self.p_s = ip

        # Initialize Values
        self.x_a_mat = np.zeros((self.n_max, self.n_max))
        self.y_a_mat = np.zeros((self.n_max, self.n_max))
        self.p_a_mat = np.zeros((self.n_max, self.n_max))
        self.x_b_mat = np.zeros(self.n_max)
        self.y_b_mat = np.zeros(self.n_max)
        self.p_b_mat = np.zeros(self.n_max)

        self.u_n_s = np.zeros(iu.shape)
        self.v_n_s = np.zeros(iv.shape)
        self.p_p = np.zeros(ip.shape)

        self.u_n = np.zeros(iu.shape)
        self.v_n = np.zeros(iv.shape)
        self.p_n = np.zeros(ip.shape)

        print('Number of Nodes = ' + str(self.n_max))

        self.solve()

        print('\nSolution Complete')

    def solve(self):

        max_diff = 1.
        it_count = 0

        while max_diff > self.convergence:

            # Solve Iteration
            self.x_mom()
            self.conserve_mass()
            self.y_mom()
            self.pressure()
            self.correct_values()

            # Calculate Max Difference
            max_diff = self.residuals()
            # max_diff = max((np.abs(self.u_s-self.u_n)/np.abs(self.u_n).max()).max(),
            #                (np.abs(self.v_s-self.v_n)/np.abs(self.v_n).max()).max(),
            #                (np.abs(self.p_s-self.p_n)/np.abs(self.p_n).max()).max())

            # Update Values
            self.u_s = self.u_n
            self.v_s = self.v_n
            self.p_s = self.p_n

            self.x_a_mat = np.zeros((self.u_s.size, self.u_s.size))
            self.y_a_mat = np.zeros((self.u_s.size, self.u_s.size))
            self.p_a_mat = np.zeros((self.u_s.size, self.u_s.size))

            # Plot Solution
            self.v_view.update(self.v_s)
            self.p_view.update(self.p_s)
            self.u_view.update(self.u_s)

            it_count += 1

            stdout.write("\r" + "Iteration # " + str(it_count))

    def x_mom(self):

        self.x_b_mat = np.zeros(self.n_max)

        # assemble coefficients
        for n_num in xrange(1, self.n_max+1):

            i, j = self.num_to_ij(n_num)

            # ADJUST FOR INDEXING
            n_id = n_num - 1

            if i == 1 or j == 1 or j == self.nj:
                self.x_a_mat[n_id][n_id] = 1.
            elif i == 2:
                self.x_b_mat[n_id] = self.bc['inlet']['value']
                self.x_a_mat[n_id][n_id] = 1.
            elif i == self.ni:
                self.x_a_mat[n_id][n_id] = 1.
                self.x_a_mat[n_id][n_id-1] = -1.
            elif self.nj/2 - 2 < j < self.nj/2 + 2 and self.ni/2 - 2 < i < self.ni/2 + 3:
                self.x_a_mat[n_id][n_id] = 1.
            else:

                # Initialize Values
                a_e, a_w, a_n, a_s = [0]*4
                a_p = 0.

                # ADJUST FOR INDEXING
                i_id = i - 1
                j_id = j - 1

                n_e_id = n_id + 1
                fe = self.rho * (self.u_s[i_id + 1][j_id] + self.u_s[i_id][j_id]) / 2.
                de = self.mu / (self.x[i_id + 1] - self.x[i_id])
                a_e = (de + max(-fe, 0.)) * (self.y[j_id + 1] - self.y[j_id])
                self.x_a_mat[n_id][n_e_id] = -a_e

                n_w_id = n_id - 1
                fw = self.rho * (self.u_s[i_id][j_id] + self.u_s[i_id - 1][j_id]) / 2.
                dw = self.mu / (self.x[i_id] - self.x[i_id - 1])
                a_w = (dw + max(fw, 0.)) * (self.y[j_id + 1] - self.y[j_id])
                self.x_a_mat[n_id][n_w_id] = -a_w

                if j == self.nj/2 + 2 and self.ni/2-2 < i < self.ni/2 + 3:
                    if self.ni/2-1 < i < self.ni/2 + 2:
                        fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                        a_p += self.mu * (self.X[i_id] - self.X[i_id - 1]) / (self.Y[j_id] - self.y[j_id])
                    else:
                        fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                        a_p += self.mu * (self.X[i_id] - self.X[i_id - 1]) / ((self.Y[j_id] - self.y[j_id])*2.)
                elif j > 2:
                    n_s_id = n_id - self.ni
                    fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                    ds = self.mu / (self.Y[j_id] - self.Y[j_id - 1])
                    a_s = (ds + max(fs, 0.)) * (self.X[i_id] - self.X[i_id - 1])
                    self.x_a_mat[n_id][n_s_id] = -a_s
                else:
                    fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                    a_p += self.mu * (self.X[i_id] - self.X[i_id-1]) / (self.Y[j_id] - self.y[j_id])

                if j == self.nj/2 - 2 and self.ni/2-2 < i < self.ni/2 + 3:
                    if self.ni/2-1 < i < self.ni/2 + 2:
                        fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                        a_p += self.mu * (self.X[i_id] - self.X[i_id - 1]) / (self.y[j_id + 1] - self.Y[j_id])
                    else:
                        fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                        a_p += self.mu * (self.X[i_id] - self.X[i_id - 1]) / ((self.y[j_id + 1] - self.Y[j_id])*2.)
                elif j < self.nj - 1:
                    n_n_id = n_id + self.ni
                    fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                    dn = self.mu / (self.Y[j_id + 1] - self.Y[j_id])
                    a_n = (dn + max(-fn, 0.)) * (self.X[i_id] - self.X[i_id - 1])
                    self.x_a_mat[n_id][n_n_id] = -a_n
                else:
                    fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                    a_p += self.mu * (self.X[i_id] - self.X[i_id-1]) / (self.y[j_id+1] - self.Y[j_id])

                a_p += a_e + a_w + a_n + a_s+(fe-fw)*(self.y[j_id+1]-self.y[j_id])+(fn-fs)*(self.X[i_id]-self.X[i_id-1])

                self.x_a_mat[n_id][n_id] = a_p / self.au

                self.x_b_mat[n_id] = (self.p_s[i_id-1][j_id]-self.p_s[i_id][j_id]) * (self.y[j_id+1] - self.y[j_id]) + \
                                     ((1.-self.au)*a_p/self.au) * self.u_s[i_id][j_id]

        # solve matrix
        a_mat = np.asmatrix(self.x_a_mat)
        b_mat = np.asmatrix(self.x_b_mat)

        self.u_n_s = np.asarray(a_mat.I*b_mat.transpose()).reshape((self.ni, self.nj), order='F')

    def y_mom(self):

        self.y_b_mat = np.zeros(self.v_s.size)

        # assemble coefficients
        for n_num in xrange(1, self.n_max + 1):

            i, j = self.num_to_ij(n_num)

            # ADJUST FOR INDEXING
            n_id = n_num - 1

            if i == 1 or j < 3 or j == self.nj:
                self.y_a_mat[n_id][n_id] = 1.
            elif i == self.ni:
                self.y_a_mat[n_id][n_id] = 1.
                self.y_a_mat[n_id][n_id - 1] = -1.
            elif self.ni/2 - 2 < i < self.ni + 2 and self.nj/2 - 3 < j < self.nj + 2:
                self.y_a_mat[n_id][n_id] = 1.
            else:

                a_e, a_w, a_n, a_s = [0] * 4
                a_p = 0.

                # ADJUST FOR INDEXING
                i_id = i - 1
                j_id = j - 1

                if self.nj/2 - 2 < j < self.nj + 3 and i == self.ni/2 - 2:
                    if j == self.nj/2 or j == self.nj/2 + 1:
                        a_p += self.mu * (self.Y[j_id] - self.Y[j_id - 1]) / (self.X[i_id] - self.x[i_id])
                        fe = self.rho * (self.u_n_s[i_id + 1][j_id] + self.u_n_s[i_id + 1][j_id - 1]) / 2.
                    else:
                        a_p += self.mu * (self.Y[j_id] - self.Y[j_id - 1]) / ((self.X[i_id] - self.x[i_id])*2.)
                        fe = self.rho * (self.u_n_s[i_id + 1][j_id] + self.u_n_s[i_id + 1][j_id - 1]) / 2.
                else:
                    n_e_id = n_id + 1
                    fe = self.rho * (self.u_n_s[i_id+1][j_id] + self.u_n_s[i_id+1][j_id-1]) / 2.
                    de = self.mu / (self.X[i_id+1] - self.X[i_id])
                    a_e = (de + max(-fe, 0.)) * (self.Y[j_id] - self.Y[j_id-1])
                    self.y_a_mat[n_id][n_e_id] = -a_e

                if self.nj/2 - 2 < j < self.nj + 3 and i == self.ni/2 + 2:
                    if j == self.nj/2 or j == self.nj/2 + 1:
                        a_p += self.mu * (self.Y[j_id] - self.Y[j_id - 1]) / (self.X[i_id] - self.x[i_id])
                        fw = self.rho * (self.u_n_s[i_id][j_id] + self.u_n_s[i_id][j_id-1]) / 2.
                    else:
                        a_p += self.mu * (self.Y[j_id] - self.Y[j_id - 1]) / ((self.X[i_id] - self.x[i_id])*2.)
                        fw = self.rho * (self.u_n_s[i_id][j_id] + self.u_n_s[i_id][j_id-1]) / 2.
                else:
                    n_w_id = n_id - 1
                    fw = self.rho * (self.u_n_s[i_id][j_id] + self.u_n_s[i_id][j_id-1]) / 2.
                    dw = self.mu / (self.X[i_id] - self.X[i_id - 1])
                    a_w = (dw + max(fw, 0.)) * (self.Y[j_id] - self.Y[j_id-1])
                    self.y_a_mat[n_id][n_w_id] = -a_w

                n_s_id = n_id - self.ni
                fs = self.rho * (self.v_s[i_id][j_id-1] + self.v_s[i_id][j_id]) / 2.
                ds = self.mu / (self.y[j_id] - self.y[j_id-1])
                a_s = (ds + max(fs, 0.)) * (self.x[i_id+1] - self.x[i_id])
                self.y_a_mat[n_id][n_s_id] = -a_s

                n_n_id = n_id + self.ni
                fn = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id][j_id+1]) / 2.
                dn = self.mu / (self.y[j_id+1] - self.y[j_id])
                a_n = (dn + max(-fn, 0.)) * (self.x[i_id+1] - self.x[i_id])
                self.y_a_mat[n_id][n_n_id] = -a_n

                a_p += a_e + a_w + a_n + a_s+(fe-fw)*(self.Y[j_id]-self.Y[j_id-1])+(fn-fs)*(self.x[i_id+1]-self.x[i_id])

                self.y_a_mat[n_id][n_id] = a_p / self.av

                self.y_b_mat[n_id] = (self.p_s[i_id][j_id-1]-self.p_s[i_id][j_id])*(self.x[i_id+1] - self.x[i_id]) + \
                                     ((1. - self.av) * a_p / self.av) * self.v_s[i_id][j_id]

        # solve matrix
        a_mat = np.asmatrix(self.y_a_mat)
        b_mat = np.asmatrix(self.y_b_mat)

        self.v_n_s = np.asarray(a_mat.I * b_mat.transpose()).reshape((self.ni, self.nj), order='F')

    def pressure(self):

        self.p_b_mat = np.zeros(self.p_s.size)

        # assemble coefficients
        for n_num in xrange(1, self.n_max + 1):

            i, j = self.num_to_ij(n_num)

            # ADJUST FOR INDEXING
            n_id = n_num - 1

            if i == 1 or i == self.ni or j == 1 or j == self.nj:
                self.p_a_mat[n_id][n_id] = 1.

            elif self.nj/2 - 2 < j < self.nj/2 + 2 and self.ni/2 - 2 < i < self.ni/2 + 2:
                self.p_a_mat[n_id][n_id] = 1.

            elif i == self.ni - 1:

                self.p_a_mat[n_id][n_id] = 1.
                self.p_a_mat[n_id][n_id+1] = -1.

            else:

                # ADJUST FOR INDEXING
                i_id = i - 1
                j_id = j - 1

                # EAST
                n_e_id = n_id + 1
                area_e = self.y[j_id+1] - self.y[j_id]

                if self.x_a_mat[n_e_id][n_e_id] == 0:
                    raise ValueError('Diverged')

                de = area_e / self.x_a_mat[n_e_id][n_e_id]
                a_e = self.rho * de * area_e

                if i == self.ni-1:
                    a_e *= self.au
                elif i < self.ni-2:
                    self.p_a_mat[n_id][n_e_id] = -a_e

                self.p_b_mat[n_id] -= self.rho * area_e * self.u_n_s[i_id + 1][j_id]

                # WEST
                n_w_id = n_id - 1
                area_w = self.y[j_id+1] - self.y[j_id]
                dw = area_w / self.x_a_mat[n_id][n_id]
                a_w = self.rho * dw * area_w

                if i == 2:
                    a_w *= self.au
                else:
                    self.p_a_mat[n_id][n_w_id] = -a_w

                self.p_b_mat[n_id] += self.rho * area_w * self.u_n_s[i_id][j_id]

                # SOUTH
                n_s_id = n_id - self.ni
                area_s = self.x[i_id+1] - self.x[i_id]
                ds = area_s / self.y_a_mat[n_id][n_id]
                a_s = self.rho * ds * area_s

                if j == 2:
                    a_s *= self.av
                else:
                    self.p_a_mat[n_id][n_s_id] = -a_s

                self.p_b_mat[n_id] += self.rho * area_s * self.v_n_s[i_id][j_id]

                # NORTH
                n_n_id = n_id + self.ni
                area_n = self.x[i_id+1] - self.x[i_id]
                dn = area_n / self.y_a_mat[n_n_id][n_n_id]
                a_n = self.rho * dn * area_n

                if j == self.nj - 1:
                    a_n *= self.av
                else:
                    self.p_a_mat[n_id][n_n_id] = -a_n

                self.p_b_mat[n_id] -= self.rho * area_n * self.v_n_s[i_id][j_id + 1]

                a_p = a_e + a_w + a_n + a_s

                self.p_a_mat[n_id][n_id] = a_p

        # solve matrix
        a_mat = np.asmatrix(self.p_a_mat)
        b_mat = np.asmatrix(self.p_b_mat)

        self.p_p = np.asarray(a_mat.I * b_mat.transpose()).reshape((self.ni, self.nj), order='F')

    def correct_values(self):

        # get corrections
        for n_num in xrange(1, self.n_max + 1):

            i, j = self.num_to_ij(n_num)

            n_id = n_num - 1
            i_id = i - 1
            j_id = j - 1

            if i > 1 and 1 < j < self.nj:
                self.u_n[i_id][j_id] = (self.u_n_s[i_id][j_id] + (self.y[j_id+1] - self.y[j_id]) *
                                        (self.p_p[i_id-1][j_id] - self.p_p[i_id][j_id]) / self.x_a_mat[n_id][n_id])

            if 1 < i < self.ni and 1 < j < self.nj:
                self.v_n[i_id][j_id] = (self.v_n_s[i_id][j_id] + (self.x[i_id+1] - self.x[i_id]) *
                                        (self.p_p[i_id][j_id-1] - self.p_p[i_id][j_id]) / self.y_a_mat[n_id][n_id])

        # adjust values
        self.p_n = self.p_s + self.ap*self.p_p

    def conserve_mass(self):

        # sum flow along inlet
        m_out = sum(self.u_n_s[-1])
        m_in = sum(self.u_n_s[1])

        self.u_n_s[-1] = self.u_n_s[-2]*(m_in/m_out)

    def num_to_ij(self, num):

        j = int(math.ceil(float(num)/float(self.ni)))
        i = num % self.ni

        if i == 0:
            i = self.ni

        return i, j

    def residuals(self):

        u_vec = np.matrix(self.u_n.reshape(self.u_n.size, order='F'))
        v_vec = np.matrix(self.v_n.reshape(self.v_n.size, order='F'))

        u_residual = np.abs(np.array(np.matrix(self.x_a_mat) * u_vec.T).reshape(self.x_b_mat.shape) - self.x_b_mat)
        v_residual = np.abs(np.array(np.matrix(self.y_a_mat) * v_vec.T).reshape(self.y_b_mat.shape) - self.y_b_mat)

        return max(u_residual.max(), v_residual.max())
