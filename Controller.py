"""
===========================
Controller for 2D CFD Simulation
===========================

This file is how the user interacts with the solver and visualizer
"""

import Viewer
import Preprocessing
import Solver


def project_solution(height, length, gs_x, gs_y, u_in):

    # Get Mesh
    u_grid, v_grid, p_grid = Preprocessing.gen_mesh(height, length, gs_x, gs_y)

    # Get Initial Values
    u, v, p = Preprocessing.initial_con(u_grid[0].shape[0], u_grid[1].shape[0], u=0.001, v=0.0001, p=0.001)

    # Get Boundary Conditions
    bc = Preprocessing.boundary_cond('velocity', 'velocity gradient', 'no slip', 'no slip', [u_in, None, None, None])

    # Create viewers
    p_viewer = Viewer.FlowContours(p, p_grid[0], p_grid[1], [0, 0, length, height], 'Pressure')
    x_v_viewer = Viewer.FlowContours(u, u_grid[0], u_grid[1], [0, 0, length, height], 'X Velocity')
    y_v_viewer = Viewer.FlowContours(v, v_grid[0], v_grid[1], [0, 0, length, height], 'Y Velocity')

    cc = .000001
    rho = 1000.
    mu = .001

    au = .5
    av = .5
    ap = .5

    Solver.Solution(p, u, v, u_grid, v_grid, p_grid, bc, cc, au, av, ap, rho, mu, p_viewer, x_v_viewer, y_v_viewer)
    Viewer.keep_open()


if __name__ == '__main__':

    volume_height = 0.01
    volume_length = 0.05
    grid_spacing_x = 0.0125
    grid_spacing_y = 0.0025
    u_inlet = 0.001

    project_solution(volume_height, volume_length, grid_spacing_x, grid_spacing_y, u_inlet)
