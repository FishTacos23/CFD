"""
===========================
Controller for 2D CFD Simulation
===========================

This file is how the user interacts with the solver and visualizer
"""

import Viewer
import Preprocessing
import Solver


def project_solution(height, length, gs, u_in):

    # Get Mesh
    u_grid, v_grid, p_grid = Preprocessing.gen_mesh(height, length, gs)

    # Get Initial Values
    u, v, p = Preprocessing.initial_con(u_grid[0].shape[0], u_grid[1].shape[0], 0, 0, 0)

    # Get Boundary Conditions
    bc = Preprocessing.boundary_cond('velocity', 'velocity gradient', 'no slip', 'no slip', [u_in, None, None, None])

    # Create viewers
    p_viewer = Viewer.FlowContours(p, p_grid[0], p_grid[1], [0, 0, length, height], 'Pressure')
    x_v_viewer = Viewer.FlowContours(u, u_grid[0], u_grid[1], [0, 0, length, height], 'X Velocity')
    y_v_viewer = Viewer.FlowContours(v, v_grid[0], v_grid[1], [0, 0, length, height], 'Y Velocity')

    # p_viewer.grid = True
    # x_v_viewer.grid = True
    # y_v_viewer.grid = True

    cc = .1
    rho = 1000.
    mu = .001

    au = .7
    av = .7
    ap = .3

    Solver.Solution(p, u, v, u_grid, v_grid, p_grid, bc, cc, au, av, ap, rho, mu, p_viewer, x_v_viewer, y_v_viewer)
    Viewer.keep_open()


if __name__ == '__main__':

    volume_height = 0.01
    volume_length = 0.05
    grid_spacing = 0.0025
    u_inlet = 0.001

    project_solution(volume_height, volume_length, grid_spacing, u_inlet)
