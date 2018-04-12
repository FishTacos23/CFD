"""
===========================
Preprocessing for 2D CFD Simulation
===========================

This file sets up the geometry and mesh of the analysis
"""

import numpy as np


def gen_mesh(height, length, gs_x, gs_y):

    """

    This function uses the geometry and grid spacing to lay out nodes in the mesh for using practice A.
    It returns the x and y location of these nodes as a numpy array

    :param height: y distance
    :type height: float
    :param length: x distance
    :type length: float
    :param gs_x: length of the square grid
    :type gs_x: float
    :param gs_y: length of the square grid
    :type gs_y: float

    :return: position of nodes in mesh
    """

    ux = np.linspace(-gs_x, length, int(length / gs_x)+2, endpoint=True)
    uy = np.linspace(-gs_y/2., height+gs_y/2., int(height / gs_y)+2, endpoint=True)

    vx = np.linspace(-gs_x/2., length+gs_x/2., int(length / gs_x)+2, endpoint=True)
    vy = np.linspace(-gs_y, height, int(height / gs_y) + 2, endpoint=True)

    px = np.linspace(-gs_x/2., length+gs_x/2., int(length / gs_x)+2, endpoint=True)
    py = np.linspace(-gs_y/2., height+gs_y/2., int(height / gs_y)+2, endpoint=True)

    return [ux, uy], [vx, vy], [px, py]


def initial_con(x_num, y_num, u=0, v=0, p=0):

    """

    This function creates an x_num by y_num set of initial values for velocity in x, velocity in y, and pressure
    It returns these values as a numpy arrays

    :param x_num: number of nodes in x
    :type x_num: int
    :param y_num: x number of nodes in y
    :type y_num: int
    :param u: starting constant value for velocity field in x
    :type u: float, function
    :param v: starting constant value for velocity field in y
    :type v: float, function
    :param p: starting constant value for pressure field
    :type p: float, function

    :return: initial u velocities, v velocities, and pressures
    """

    if callable(u):
        u_start = u()
    else:
        u_start = u

    if callable(v):
        v_start = v()
    else:
        v_start = v

    if callable(p):
        p_start = p()
    else:
        p_start = p

    pressure = np.zeros((x_num, y_num), dtype=float) + p_start
    u_vel = np.zeros((x_num, y_num), dtype=float) + u_start
    v_vel = np.zeros((x_num, y_num), dtype=float) + v_start

    return u_vel, v_vel, pressure


def boundary_cond(inlet, outlet, top, bottom, values):

    """

    :param inlet: type of inlet condition
    :type inlet: str
    :param outlet: type of outlet condition
    :type outlet: str
    :param top: type of top condition
    :type top: str
    :param bottom: type of bottom condition
    :type bottom: str
    :param values: list of values for corresponding boundary conditions, if no values necessary use None
    :type values: list

    :return: dictionary of boundary condition dictionaries
    """

    boundary_conditions = {'inlet': None, 'outlet': None, 'top': None, 'bottom': None}

    if inlet == 'velocity':
        boundary_conditions['inlet'] = {'type': inlet, 'value': values[0]}
    if outlet == 'velocity gradient':
        boundary_conditions['outlet'] = {'type': outlet, 'value': values[1]}
    if top == 'no slip':
        boundary_conditions['top'] = {'type': top, 'value': values[2]}
    if bottom == 'no slip':
        boundary_conditions['bottom'] = {'type': bottom, 'value': values[3]}

    return boundary_conditions
