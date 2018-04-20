"""
===========================
Controller for 2D CFD Simulation
===========================

This file is how the user interacts with the solver and visualizer
"""

import Viewer
import Preprocessing
import Solver


def project_solution(height, length, gs_x, gs_y, u_in, cc, au, av, ap, rho, mu):

    # Get Mesh
    u_grid, v_grid, p_grid = Preprocessing.gen_mesh(height, length, gs_x, gs_y)

    # Get Initial Values
    u, v, p = Preprocessing.initial_con(u_grid[0].shape[0], u_grid[1].shape[0], u=0., v=0., p=0.)

    # Get Boundary Conditions
    bc = Preprocessing.boundary_cond('velocity', 'velocity gradient', 'no slip', 'no slip', [u_in, None, None, None])

    # Create viewers
    p_viewer = Viewer.FlowContours(p, p_grid[0], p_grid[1], [0, 0, length, height], 'Pressure')
    x_v_viewer = Viewer.FlowContours(u, u_grid[0], u_grid[1], [0, 0, length, height], 'X Velocity')
    y_v_viewer = Viewer.FlowContours(v, v_grid[0], v_grid[1], [0, 0, length, height], 'Y Velocity')

    s = Solver.Solution(p, u, v, u_grid[0], v_grid[0], v_grid[1], u_grid[1], bc, cc, au, av, ap, rho, mu, p_viewer,
                        x_v_viewer, y_v_viewer)

    p = s.p_n[s.ni/2:-1, s.nj/2]
    dp = p[-1]-p[0]

    print 'dp/dx = ' + str(dp/(gs_x*(len(p)-1)))
    print 'max U = ' + str(s.u_n[s.ni/2:-1, s.nj/2].max())
    print 'max V = ' + str(s.v_n[s.ni/2:-1, s.nj/2].max())

    Viewer.keep_open()


if __name__ == '__main__':

    volume_height = 0.01
    volume_length = 0.05
    grid_spacing_x = 0.003
    grid_spacing_y = 0.0006
    u_inlet = 0.001
    criteria = .00000001
    density = 1000.
    viscosity = .001

    p_relax = .5
    u_relax = .5
    v_relax = .5

    project_solution(volume_height, volume_length, grid_spacing_x, grid_spacing_y, u_inlet, criteria, u_relax, v_relax,
                     p_relax, density, viscosity)
