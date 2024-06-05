#!/usr/bin/env python 

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

from problems.no_slip_box import *
from solvers.chorin_projection import *


def cfl_number(u, dim, dt, n):
    max_u = norm(u.vector(), 'linf')
    cfl = dim * max_u * (dt/n)
    return cfl

def init_problem(name, parameters):
    if name == "no_slip_box":
        return no_slip_box(parameters)

def init_ns_solver(name, problem):
    if name == "chorin_projection":
        return chorin_projection(problem)

def calculate_error(n_plus, n_minus, psi, u, n_plus_, n_minus_, psi_, u_):

    #Calculte the L-infinity error of the solution to the tentative estimate for the charges
    error_plus = norm(n_plus.vector() - n_plus_.vector(), 'linf')
    error_minus = norm(n_minus.vector() - n_minus_.vector(), 'linf')

    #Calculte the H^1 error of the solution to the tentative estimate for the electric field psi
    error_psi = errornorm(psi, psi_, 'H10')

    #Calculate the L^2 error of the solution to the tentative estimate for the velocity u
    error_u = errornorm(u, u_, 'L2')

    #Calculte the error
    error = error_u + error_psi + error_plus + error_minus

    return error