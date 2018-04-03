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
    x, y = Preprocessing.gen_mesh(height, length, gs)

    # Get initial values
    u, v, p = Preprocessing.initial_con(x.shape[0], y.shape[0], 0, 0, 0)
    bc = Preprocessing.boundary_cond('velocity', 'velocity gradient', 'no slip', 'no slip', [u_in, None, None, None])

    # Create viewers
    pressure_viewer = Viewer.FlowContours(p, x, y, 'Pressure')
    x_velocity_viewer = Viewer.FlowContours(u, x, y, 'X Velocity')
    y_velocity_viewer = Viewer.FlowContours(v, x, y, 'Y Velocity')

    Solver.Solution(p, u, v, bc, pressure_viewer, x_velocity_viewer, y_velocity_viewer)


if __name__ == '__main__':

    volume_height = 0.01
    volume_length = 0.05
    grid_spacing = 0.0025
    u_inlet = 0.001

    project_solution(volume_height, volume_length, grid_spacing, u_inlet)
