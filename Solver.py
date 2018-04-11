"""
===========================
Solver for 2D CFD Simulation
===========================

This file uses the control volume method with practice a with the simple scheme to solve a 2d flow field
"""
import numpy as np


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

        # Flow Properties
        self.bc = bc

        self.rho = rho
        self.mu = mu

        self.convergence = cc
        self.au = au
        self.av = av
        self.ap = ap

        # grid locations
        self.li = un[0]
        self.bj = un[1]
        self.bi = vn[0]
        self.lj = vn[1]

        # Star Values
        self.u_s = iu
        self.v_s = iv
        self.p_s = ip

        # Initialize Matrix Values
        self.x_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.y_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.p_a_mat = np.zeros((self.u_s.size, self.u_s.size))
        self.p_d = np.zeros((self.u_s.size, self.u_s.size))
        self.b_mat = np.zeros(self.u_s.size)

        # New Star Values
        self.u_n_s = np.empty(iu.shape)
        self.v_n_s = np.empty(iv.shape)
        self.p_p = np.empty(ip.shape)

        # New Values
        self.u_n = np.empty(iu.shape)
        self.v_n = np.empty(iv.shape)
        self.p_n = np.empty(ip.shape)

        self.solve_message = 'SOLUTION COMPLETE'

        self.solve()

    def solve(self):

        max_diff = 1.

        while max_diff > self.convergence:

            self.x_mom()
            self.conserve_mass()
            self.y_mom()
            self.pressure()
            self.correct_values()

            max_diff = max(np.max(np.abs(self.u_s-self.u_n)), np.max(np.abs(self.v_s-self.v_n)),
                           np.max(np.abs(self.p_s-self.p_n)), np.max(np.abs(self.u_n[:][-1]-self.u_n[:][-2])))

            self.u_s = self.u_n
            self.v_s = self.v_n
            self.p_s = self.p_n

            self.v_view.update(self.v_s)
            self.p_view.update(self.p_s)
            self.u_view.update(self.u_s)

            self.x_a_mat = np.zeros((self.u_s.size, self.u_s.size))
            self.y_a_mat = np.zeros((self.u_s.size, self.u_s.size))
            self.p_a_mat = np.zeros((self.u_s.size, self.u_s.size))
            self.b_mat = np.zeros(self.u_s.size)

        print self.solve_message

    def x_mom(self):

        self.b_mat = np.zeros(self.u_s.size)

        # assemble coefficients
        for i in xrange(self.u_s.shape[1]):  # columns
            for j in xrange(self.u_s.shape[0]):  # rows

                n_num = j*self.u_s.shape[1] + i

                if i == 0 or j == 0 or j == self.u_s.shape[0] - 1:
                    pass

                elif i == 1:

                    self.b_mat[n_num] = self.bc['inlet']['value']
                    self.x_a_mat[n_num][n_num] = 1.
                else:

                    n_e = -1
                    n_w = -1
                    n_n = -1
                    n_s = -1

                    a_e = 0
                    a_w = 0
                    a_n = 0
                    a_s = 0

                    fe = 0
                    fw = 0
                    fn = 0
                    fs = 0

                    a_p = 0.

                    if i < self.u_s.shape[1]-1:  # If not on the right face
                        n_e = n_num+1
                    if i > 0:  # If not on the left face
                        n_w = n_num-1
                    if j > 0:  # If not on bottom face
                        n_s = n_num - self.u_s.shape[1]
                    if j < self.u_s.shape[0] - 1:  # If not on top face
                        n_n = n_num + self.u_s.shape[1]

                    if n_e >= 0:
                        fe = self.rho * (self.u_s[j][i + 1] + self.u_s[j][i]) / 2.
                        de = self.mu / (self.li[i + 1] - self.li[i])
                        a_e = (de + max(-fe, 0.)) * (self.lj[j + 1] - self.lj[j])
                        self.x_a_mat[n_num][n_e] = a_e
                    if n_w >= 0:
                        fw = self.rho * (self.u_s[j][i] + self.u_s[j][i - 1]) / 2.
                        dw = self.mu / (self.li[i] - self.li[i - 1])
                        a_w = (dw + max(fw, 0.)) * (self.lj[j + 1] - self.lj[j])
                        self.x_a_mat[n_num][n_w] = a_w
                    if n_n >= 0:
                        fn = self.rho * (self.v_s[j + 1][i] + self.v_s[j+1][i - 1]) / 2.
                        dn = self.mu / (self.bj[j + 1] - self.bj[j])
                        a_n = (dn + max(-fn, 0.)) * (self.bi[i] - self.bi[i-1])
                        self.x_a_mat[n_num][n_n] = a_n
                    else:
                        a_p = self.mu * (self.bi[i+1]-self.bi[i]) / (self.bj[j]-self.lj[j])
                    if n_s >= 0:
                        fs = self.rho * (self.v_s[j][i] + self.v_s[j][i - 1]) / 2.
                        ds = self.mu / (self.bj[j] - self.bj[j - 1])
                        a_s = (ds + max(-fs, 0.)) * (self.bi[i] - self.bi[i-1])
                        self.x_a_mat[n_num][n_s] = a_s
                    else:
                        a_p = self.mu * (self.bi[i + 1] - self.bi[i]) / (self.lj[j + 1] - self.bj[j])

                    a_p += a_e+a_w+a_n+a_s+(fe-fw)*(
                            self.lj[j+1]-self.lj[j])+(fn-fs)*(self.bi[i]-self.bi[i-1])

                    self.x_a_mat[n_num][n_num] = a_p

        # solve matrix

        del_list = [i for i in xrange(self.u_s.shape[1])]
        for i in xrange(self.u_s.size, self.u_s.size - self.u_s.shape[1], -1):
            del_list.append(i)
        for i in xrange(0, self.u_s.size, self.u_s.shape[1]):
            del_list.append(i)

        a_mat = np.delete(self.x_a_mat, del_list, 0)
        a_mat = np.delete(a_mat, del_list, 1)
        b_mat = np.delete(self.b_mat, del_list)

        a_mat = np.asmatrix(a_mat)
        b_mat = np.asmatrix(b_mat)

        uns = np.asarray(a_mat.I*b_mat.transpose())

        k = 0
        # assemble coefficients
        for j in xrange(self.u_s.shape[0]):  # rows
            for i in xrange(self.u_s.shape[1]):  # columns
                n_num = j * self.u_s.shape[1] + i

                if n_num in del_list:
                    self.u_n_s[j][i] = 0
                else:
                    self.u_n_s[j][i] = uns[k]
                    k += 1

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
                        de = self.mu / (self.bi[i] - self.bi[i - 1])
                        a_e = (de + max(-fe, 0.)) * (self.bj[j] - self.bj[j - 1])
                        self.y_a_mat[n_num][n_e] = a_e
                    if n_w >= 0:
                        fw = self.rho * (self.u_n_s[j][i] + self.u_n_s[j - 1][i]) / 2.
                        dw = self.mu / (self.bi[i + 1] - self.bi[i])
                        a_w = (dw + max(fw, 0.)) * (self.bj[j] - self.bj[j - 1])
                        self.y_a_mat[n_num][n_w] = a_w
                    if n_n >= 0:
                        fn = self.rho * (self.v_s[j][i] + self.v_s[j + 1][i]) / 2.
                        dn = self.mu / (self.lj[j + 1] - self.lj[j])
                        a_n = (dn + max(-fn, 0.)) * (self.li[i + 1] - self.li[i])
                        self.y_a_mat[n_num][n_n] = a_n
                    if n_s >= 0:
                        fs = self.rho * (self.v_s[j - 1][i] + self.v_s[j][i]) / 2.
                        ds = self.mu / (self.lj[j] - self.lj[j - 1])
                        a_s = (ds + max(-fs, 0.)) * (self.li[i + 1] - self.li[i])
                        self.y_a_mat[n_num][n_s] = a_s

                    a_p = a_e + a_w + a_n + a_s + (fe - fw) * (
                            self.bj[j] - self.bj[j-1]) + (fn - fs) * (self.li[i + 1] - self.li[i])

                    self.y_a_mat[n_num][n_num] = a_p

        # solve matrix

        del_list = [i for i in xrange(self.v_s.shape[1]*2)]
        for i in xrange(self.v_s.size, self.v_s.size - self.v_s.shape[1], -1):
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

                if i == 0 or j == 0 or n_n == 0:
                    self.p_a_mat[n_num][n_num] = 0.
                else:

                    a_e = 0
                    a_w = 0
                    a_n = 0
                    a_s = 0

                    if n_e >= 0:
                        area_e = self.lj[j+1] - self.lj[j]
                        de = area_e*self.au/self.x_a_mat[n_num][n_e]
                        a_e = self.rho*de*area_e
                        self.p_a_mat[n_num][n_e] = a_e
                        self.p_d[n_num][n_e] = de
                        self.b_mat[n_num] -= self.rho*area_e*self.u_n_s[j][i+1]
                    if n_w >= 0:
                        area_w = self.lj[j + 1] - self.lj[j]
                        dw = area_w * self.au / self.x_a_mat[n_num][n_num]
                        a_w = self.rho * dw * area_w
                        self.p_a_mat[n_num][n_w] = a_w
                        self.p_d[n_num][n_w] = dw
                        self.b_mat[n_num] += self.rho * area_w * self.u_n_s[j][i]
                    if n_n >= 0:
                        area_n = self.li[i + 1] - self.li[1][i]
                        dn = area_n * self.av / self.y_a_mat[n_num][n_n]
                        a_n = self.rho * dn * area_n
                        self.p_a_mat[n_num][n_n] = a_n
                        self.p_d[n_num][n_n] = dn
                        self.b_mat[n_num] -= self.rho * area_n * self.v_n_s[j+1][i]
                    if n_s >= 0:
                        area_s = self.li[i + 1] - self.li[i]
                        ds = area_s * self.av / self.y_a_mat[n_num][n_s]
                        a_s = self.rho * ds * area_s
                        self.p_a_mat[n_num][n_s] = a_s
                        self.p_d[n_num][n_s] = ds
                        self.b_mat[n_num] += self.rho * area_s * self.v_n_s[j][i]
                    a_p = a_e + a_w + a_n + a_s

                    self.p_a_mat[n_num][n_num] = a_p

        # solve matrix
        a_mat = np.asmatrix(self.p_a_mat)
        b_mat = np.asmatrix(self.b_mat)

        self.p_p = np.asarray(a_mat.I * b_mat).reshape(self.p_s.shape)

    def correct_values(self):

        # get corrections
        for i in xrange(self.p_s.shape[1]):
            for j in xrange(self.p_s.shape[0]):
                if i == 0:
                    pass
                else:
                    self.u_n[i][j] = self.u_n_s + self.p_d[i][j]*(self.p_p[i-1][j] - self.p_p[i][j])
                    self.v_n[i][j] = self.v_n_s + self.p_d[i][j]*(self.p_p[i][j-1] - self.p_p[i][j])

        # adjust values
        self.p_n = self.p_s + self.ap*self.p_p

    def conserve_mass(self):

        # sum flow along inlet
        m_out = sum(self.u_n_s[:, -1])
        m_in = sum(self.u_n_s[:, 1])

        self.u_n_s[:, -1] = self.u_n_s[:, -2]*(m_in/m_out)
