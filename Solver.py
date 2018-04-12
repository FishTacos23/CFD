"""
===========================
Solver for 2D CFD Simulation
===========================

This file uses the control volume method with practice a with the simple scheme to solve a 2d flow field
"""
import numpy as np
import math


class Solution:

    def __init__(self, ip, iu, iv, un, vn, pn, bc, cc, au, av, ap, rho, mu, p_view=None, u_view=None, v_view=None):

        """

        This class solves the Flow Field Simulations

        :param ip: initial pressures
        :type ip: ndarray
        :param iu: initial x velocities
        :type iu: ndarray
        :param iv: initial y velocities
        :type iv: ndarray
        :param un: x velocity node locations
        :type un: ndarray
        :param vn: y velocity node locations
        :type vn: ndarray
        :param pn: pressure grid locations
        :type pn: ndarray
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

        self.x = un[0]
        self.X = vn[0]
        self.Y = un[1]
        self.y = vn[1]

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
        self.x_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.y_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.p_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.p_d = np.zeros((self.u_s.size, self.u_s.size))
        self.b_mat = np.zeros(self.u_s.size)

        self.u_n_s = np.empty(iu.shape)
        self.v_n_s = np.empty(iv.shape)
        self.p_p = np.empty(ip.shape)

        self.u_n = np.empty(iu.shape)
        self.v_n = np.empty(iv.shape)
        self.p_n = np.empty(ip.shape)

        self.solve()

        print('Solution Complete')

    def solve(self):

        max_diff = 1.

        while max_diff > self.convergence:

            # Solve Iteration
            self.x_mom()
            self.conserve_mass()
            self.y_mom()
            self.pressure()
            self.correct_values()

            # Calculate Max Difference
            max_diff = max(np.max(np.abs(self.u_s-self.u_n)), np.max(np.abs(self.v_s-self.v_n)),
                           np.max(np.abs(self.p_s-self.p_n)))

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

    def x_mom(self):

        self.b_mat = np.zeros(self.u_s.size)

        # assemble coefficients
        for n_num in xrange(1, self.n_max+1):

            i, j = self.num_to_ij(n_num)

            # ADJUST FOR INDEXING
            n_id = n_num - 1

            if i == 1 or j == 1 or j == self.nj:
                self.x_a_mat[n_id][n_id] = 1.
            elif i == 2:
                self.b_mat[n_id] = self.bc['inlet']['value']
                self.x_a_mat[n_id][n_id] = 1.
            elif i == self.ni:
                self.x_a_mat[n_id][n_id] = 1.
                self.x_a_mat[n_id][n_id-1] = -1.
            else:

                # Initialize Values
                a_e, a_w, a_n, a_s = [0]*4
                fe, fw, fn, fs = [0]*4
                a_p = 0.

                # ADJUST FOR INDEXING
                i_id = i - 1
                j_id = j - 1

                if i < self.ni:
                    n_e_id = n_id + 1
                    fe = self.rho * (self.u_s[i_id + 1][j_id] + self.u_s[i_id][j_id]) / 2.
                    de = self.mu / (self.x[i_id + 1] - self.x[i_id])
                    a_e = (de + max(-fe, 0.)) * (self.y[j_id + 1] - self.y[j_id])
                    self.x_a_mat[n_id][n_e_id] = -a_e
                if i > 1:
                    n_w_id = n_id - 1
                    fw = self.rho * (self.u_s[i_id][j_id] + self.u_s[i_id - 1][j_id]) / 2.
                    dw = self.mu / (self.x[i_id] - self.x[i_id - 1])
                    a_w = (dw + max(fw, 0.)) * (self.y[j_id + 1] - self.y[j_id])
                    self.x_a_mat[n_id][n_w_id] = -a_w
                if j > 2:
                    n_s_id = n_id - self.ni
                    fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                    ds = self.mu / (self.Y[j_id] - self.Y[j_id - 1])
                    a_s = (ds + max(fs, 0.)) * (self.X[i_id] - self.X[i_id - 1])
                    self.x_a_mat[n_id][n_s_id] = -a_s
                elif j == 2:
                    fs = self.rho * (self.v_s[i_id][j_id] + self.v_s[i_id - 1][j_id]) / 2.
                    a_p += self.mu * (self.X[i_id] - self.X[i_id-1]) / (self.Y[j_id] - self.y[j_id])
                if j < self.nj - 1:
                    n_n_id = n_id + self.ni
                    fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                    dn = self.mu / (self.Y[j_id + 1] - self.Y[j_id])
                    a_n = (dn + max(-fn, 0.)) * (self.X[i_id] - self.X[i_id - 1])
                    self.x_a_mat[n_id][n_n_id] = -a_n
                elif j == self.nj - 1:
                    fn = self.rho * (self.v_s[i_id][j_id + 1] + self.v_s[i_id - 1][j_id + 1]) / 2.
                    a_p += self.mu * (self.X[i_id] - self.X[i_id-1]) / (self.y[j_id+1] - self.Y[j_id])

                a_p += a_e + a_w + a_n + a_s+(fe-fw)*(self.y[j_id+1]-self.y[j_id])+(fn-fs)*(self.X[i_id]-self.X[i_id-1])

                self.x_a_mat[n_id][n_id] = a_p / self.au

                self.b_mat[n_id] = (self.p_s[i_id-1][j_id]-self.p_s[i_id][j_id]) * (self.y[j_id+1] - self.y[j_id]) + \
                                   ((1.-self.au)*a_p/self.au) * self.u_s[i_id][j_id]

        # solve matrix
        a_mat = np.asmatrix(self.x_a_mat)
        b_mat = np.asmatrix(self.b_mat)

        self.u_n_s = np.asarray(a_mat.I*b_mat.transpose()).reshape((self.ni, self.nj), order='F')

    def y_mom(self):

        self.b_mat = np.zeros(self.v_s.size)

        # assemble coefficients
        for i in xrange(self.v_s.shape[1]):
            for j in xrange(self.v_s.shape[0]):

                n_num = j*self.u_s.shape[1] + i

                n_e = -1
                n_w = -1
                n_n = -1
                n_s = -1

                if i < self.u_s.shape[1] - 1:  # If not on the right face
                    n_e = n_num + 1
                if i > 0:  # If not on the left face
                    n_w = n_num - 1
                if j > 0:  # If not on bottom face
                    n_s = n_num - self.u_s.shape[1]
                if j < self.u_s.shape[0] - 1:  # If not on top face
                    n_n = n_num + self.u_s.shape[1]

                if i == 0 or j == 0 or n_n < 0 or n_e < 0:
                    self.y_a_mat[n_num][n_num] = 0.
                else:

                    a_e = 0
                    a_w = 0
                    a_n = 0
                    a_s = 0

                    fe = 0
                    fw = 0
                    fn = 0
                    fs = 0

                    if n_e >= 0:
                        fe = self.rho * (self.u_n_s[j][i + 1] + self.u_n_s[j-1][i + 1]) / 2.
                        de = self.mu / (self.X[i + 1] - self.X[i])
                        a_e = (de + max(-fe, 0.)) * (self.Y[j] - self.Y[j - 1])
                        self.y_a_mat[n_num][n_e] = a_e
                    if n_w >= 0:
                        fw = self.rho * (self.u_n_s[j][i] + self.u_n_s[j - 1][i]) / 2.
                        dw = self.mu / (self.X[i] - self.X[i - 1])
                        a_w = (dw + max(fw, 0.)) * (self.Y[j] - self.Y[j - 1])
                        self.y_a_mat[n_num][n_w] = a_w
                    if n_n >= 0:
                        fn = self.rho * (self.v_s[j][i] + self.v_s[j + 1][i]) / 2.
                        dn = self.mu / (self.y[j + 1] - self.y[j])
                        a_n = (dn + max(-fn, 0.)) * (self.x[i + 1] - self.x[i])
                        self.y_a_mat[n_num][n_n] = a_n
                    if n_s >= 0:
                        fs = self.rho * (self.v_s[j - 1][i] + self.v_s[j][i]) / 2.
                        ds = self.mu / (self.y[j] - self.y[j - 1])
                        a_s = (ds + max(-fs, 0.)) * (self.x[i + 1] - self.x[i])
                        self.y_a_mat[n_num][n_s] = a_s

                    a_p = a_e + a_w + a_n + a_s + (fe - fw) * (
                            self.Y[j] - self.Y[j - 1]) + (fn - fs) * (self.x[i + 1] - self.x[i])

                    self.y_a_mat[n_num][n_num] = a_p

        # solve matrix

        del_list = [i for i in xrange(self.v_s.shape[1]*2)]
        for i in xrange(self.v_s.size-1, self.v_s.size - self.v_s.shape[1], -1):
            del_list.append(i)
        for i in xrange(0, self.v_s.size, self.v_s.shape[1]):
            del_list.append(i)
        for i in xrange(self.v_s.shape[1]-1, self.v_s.size, self.v_s.shape[1]):
            del_list.append(i)

        a_mat = np.delete(self.y_a_mat, del_list, 0)
        a_mat = np.delete(a_mat, del_list, 1)
        b_mat = np.delete(self.b_mat, del_list)

        a_mat = np.asmatrix(a_mat)
        b_mat = np.asmatrix(b_mat)

        vns = np.asarray(a_mat.I * b_mat.transpose())

        k = 0
        # assemble coefficients
        for j in xrange(self.v_s.shape[0]):  # rows
            for i in xrange(self.v_s.shape[1]):  # columns
                n_num = j * self.v_s.shape[1] + i

                if n_num in del_list:
                    self.v_n_s[j][i] = 0
                else:
                    self.v_n_s[j][i] = vns[k]
                    k += 1

    def pressure(self):

        self.b_mat = np.zeros(self.p_s.size)

        # assemble coefficients
        for i in xrange(self.p_s.shape[1]):
            for j in xrange(self.p_s.shape[0]):

                n_num = j * self.u_s.shape[1] + i

                n_e = -1
                n_w = -1
                n_n = -1
                n_s = -1

                if i < self.u_s.shape[1] - 1:  # If not on the right face
                    n_e = n_num + 1
                if i > 0:  # If not on the left face
                    n_w = n_num - 1
                if j > 0:  # If not on bottom face
                    n_s = n_num - self.u_s.shape[1]
                if j < self.u_s.shape[0] - 1:  # If not on top face
                    n_n = n_num + self.u_s.shape[1]

                if i == 0 or j == 0 or n_n < 0 or n_e < 0:
                    self.p_a_mat[n_num][n_num] = 0.
                else:

                    a_e = 0
                    a_w = 0
                    a_n = 0
                    a_s = 0

                    if n_e >= 0:
                        area_e = self.y[j + 1] - self.y[j]
                        de = area_e*self.au/self.x_a_mat[n_e][n_e]  # TODO this is giving a divide by zero warning
                        a_e = self.rho*de*area_e
                        self.p_a_mat[n_num][n_e] = a_e
                        self.p_d[n_num][n_e] = de
                        self.b_mat[n_num] -= self.rho*area_e*self.u_n_s[j][i+1]
                    if n_w >= 0:
                        area_w = self.y[j + 1] - self.y[j]
                        dw = area_w * self.au / self.x_a_mat[n_num][n_num]
                        a_w = self.rho * dw * area_w
                        self.p_a_mat[n_num][n_w] = a_w
                        self.p_d[n_num][n_w] = dw
                        self.b_mat[n_num] += self.rho * area_w * self.u_n_s[j][i]
                    if n_n >= 0:
                        area_n = self.x[i + 1] - self.x[i]
                        dn = area_n * self.av / self.y_a_mat[n_n][n_n]
                        a_n = self.rho * dn * area_n
                        self.p_a_mat[n_num][n_n] = a_n
                        self.p_d[n_num][n_n] = dn
                        self.b_mat[n_num] -= self.rho * area_n * self.v_n_s[j+1][i]
                    if n_s >= 0:
                        area_s = self.x[i + 1] - self.x[i]
                        ds = area_s * self.av / self.y_a_mat[n_num][n_s]
                        a_s = self.rho * ds * area_s
                        self.p_a_mat[n_num][n_s] = a_s
                        self.p_d[n_num][n_s] = ds
                        self.b_mat[n_num] += self.rho * area_s * self.v_n_s[j][i]
                    a_p = a_e + a_w + a_n + a_s

                    self.p_a_mat[n_num][n_num] = a_p

        # solve matrix

        del_list = [i for i in xrange(self.p_s.shape[1])]
        for i in xrange(self.p_s.size - 1, self.p_s.size - self.p_s.shape[1], -1):
            del_list.append(i)
        for i in xrange(0, self.p_s.size, self.p_s.shape[1]):
            del_list.append(i)
        for i in xrange(self.p_s.shape[1] - 1, self.p_s.size, self.p_s.shape[1]):
            del_list.append(i)

        a_mat = np.delete(self.p_a_mat, del_list, 0)
        a_mat = np.delete(a_mat, del_list, 1)
        b_mat = np.delete(self.b_mat, del_list)

        a_mat = np.asmatrix(a_mat)
        b_mat = np.asmatrix(b_mat)

        pns = np.asarray(a_mat.I * b_mat.transpose())

        k = 0
        # assemble coefficients
        for j in xrange(self.v_s.shape[0]):  # rows
            for i in xrange(self.v_s.shape[1]):  # columns
                n_num = j * self.v_s.shape[1] + i

                if n_num in del_list:
                    self.p_p[j][i] = 0
                else:
                    self.p_p[j][i] = pns[k]
                    k += 1

    def correct_values(self):

        # get corrections
        for i in xrange(self.p_s.shape[1]):
            for j in xrange(self.p_s.shape[0]):
                if i == 0:
                    pass
                else:
                    self.u_n[j][i] = self.u_n_s[j][i] + self.p_d[j][i]*(self.p_p[j][i-1] - self.p_p[j][i])
                    self.v_n[j][i] = self.v_n_s[j][i] + self.p_d[j][i]*(self.p_p[j-1][i] - self.p_p[j][i])

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
