#!/usr/bin/env python 

"""
    This is a problem for the Navier--Stokes equation, 

    dt u + (u . nabla) u - Delta u + grad p = f
    nabla . u = 0  
    with no-slip boundary condition boundary

"""
from fenics import *
import numpy as np
tol = 1e-10

import sys
sys.path.insert(0,"..") 
from problems.director import *
from problems.base_problem import base_problem
from problems.ics import *


class no_slip_box(base_problem):

    def __init__(self, parameters):
        self.name = "no_slip_box"
        self.dim = parameters["dim"]
        self.dt = parameters["dt"]
        self.n = parameters["meshsize"]
        self.mu = parameters["mu"]
        self.rey = parameters["rey"]
        self.psi_freq = parameters["psi_freq"]
        self.psi_amplitude = parameters["bc_psi_amplitude"]
        self.e_orth = parameters["e_orth"]
        self.e_para = parameters["e_para"]
        self.lamb_para = parameters["lamb_para"]
        self.lamb_orth = parameters["lamb_orth"]
        self.B = parameters["B"]
        self.F = parameters["F"]

        self.mesh = RectangleMesh(Point(0,0),Point(1.5,1),int(1.5 * self.n),self.n)

        #Spaces for Chorin projection
        self.V = VectorFunctionSpace(self.mesh, "P", 2)
        self.Q = FunctionSpace(self.mesh, "P", 1)

        #Define the initial conditions for Chorin projection
        self.u_i = ics(parameters["u_i"], self.V).ic
        self.p_i = ics(parameters["p_i"], self.Q).ic

        #Define initial condtions for the Nernst--Planck--Poisson system
        self.psi_i = ics(parameters["psi_i"], self.Q).ic
        self.n_plus_i = ics(parameters["n_plus_i"], self.Q).ic
        self.n_minus_i = ics(parameters["n_minus_i"], self.Q).ic

        #Define the director
        self.d = project(create_director(parameters["director"]), self.V)

        self.function_names = ['velocity', 'pressure', 'electric_potential', 'n_plus', 'n_minus', 'director']
        self.functions = [self.u_i, self.p_i, self.psi_i, self.n_plus_i, self.n_minus_i, self.d]

        #Define the boundary
        self.boundary = 'near(x[0], 0) || near(x[0], 1.5) || near(x[1],0) || near(x[1],1)'
        self.left = 'near(x[0],0)'
        self.right = 'near(x[0], 1.5)'
        self.top_n_bottom = 'near(x[1], 0) || near(x[1], 1)'

        #Define Dirichlet boundary conditions for u
        self.bc_u = DirichletBC(self.V, Constant((0,0)), self.boundary)

        #Define the boundary conditions for the electric field psi
        self.tau = parameters["tau"]
        self.xi = Expression('x[1] <= tol || x[1] >= 1 - tol ? 0 : a * sin(2*pi*n*t + pi/2) * (1 - (2/1.5) * x[0])', a = self.psi_amplitude, n = self.psi_freq, t = 0, pi = np.pi, degree = 4, tol = tol)
        self.xi_time_derivative = Expression('x[1] <= tol || x[1] >= 1 - tol ? 0 : a * 2*pi*n * cos(2*pi*n*t + pi/2) * (1 - (2/1.5) * x[0])', a = self.psi_amplitude, n = self.psi_freq, t = 0, pi = np.pi, degree = 4, tol = tol)

    def time_update(self, t):
        if self.psi_freq > 0:
            self.xi.t += self.dt
            self.xi_time_derivative.t += self.dt
