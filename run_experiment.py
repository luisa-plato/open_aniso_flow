from fenics import *
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10,10)
import pandas as pd
from tqdm import tqdm

import sys
sys.path.insert(0,"..")

from problems.no_slip_box import *
from solvers.chorin_projection import *
from solvers.poisson import * 
from solvers.nernst_planck import *
from common.save_data import *
from common.extras import *
import telegram_send

set_log_level(50)
tol = 1E-14

def run_experiment(parameters):

    #Set the final time and the time-step size
    T = parameters["T"]
    dt = parameters["dt"]
    num_steps = int(T/dt)
    dim = parameters["dim"]
    freq = parameters["save_freq"]

    e_para = parameters["e_para"]

    #Set the tolerance and the maximum number of iterations, when the fixed point solver is selected
    #If sequential is selected finish the computation after one iteration
    if parameters["coupling"] == "sequential":
        max_iter = 1
        fp_tol = 0
    elif parameters["coupling"] == "fp_solver":
        max_iter = parameters["max_iter"]
        fp_tol = parameters["fp_tol"]

    #Fix the meshsize
    n = parameters["meshsize"]

    #Set the initial time to zero
    t = 0

    #Initialize the given problem
    problem = init_problem(parameters["problem"], parameters)

    #Initialization of the solvers
    ns_solver = init_ns_solver(parameters["ns_solver"], problem)
    p_solver = poisson(problem)
    np_solver = nernst_planck(problem)

    #Write parameter to a text file
    path = create_dic(problem.name)
    write_param_file(path, parameters)

    #Create vtk files
    vtk_files = create_vtk(path, problem.function_names)

    #Write the initial data to the created vtk_files
    write_to_vtk(vtk_files, problem.functions, t)

    #Define tentative solution functions
    u_ = Function(problem.V)
    p_ = Function(problem.Q)
    n_plus_ = Function(problem.Q)
    n_minus_ = Function(problem.Q)
    psi_ = Function(problem.Q)

    #Calculate initial electric potential
    problem.psi_i.assign(p_solver.solve(problem.B * (problem.n_plus_i - problem.n_minus_i), d = problem.d))

    #Create empty lists to plot the cfl number over time
    cfl_numbers = [0]
    times = [t]

    #Start the time-stepping loop
    for i in tqdm(range(num_steps)):

        # Update current time
        t += dt
        times.append(t)

        #Update the boundary condtions if necessary
        problem.time_update(t)

        #Set default value for the error of the fixed point solver to one
        error = 1

        #Set the index of the iteration of the fixed point solver to 0
        iter = 0

        #Set tentative solution to solution at previous time step
        u_.assign(problem.u_i)
        p_.assign(problem.p_i)
        n_plus_.assign(problem.n_plus_i)
        n_minus_.assign(problem.n_minus_i)
        psi_.assign(problem.psi_i)

        #Start loop for the fixed point solver
        while iter < max_iter and error > fp_tol:

            iter += 1

            #Solve Navier--Stokes equation
            ns_force = - problem.B * (n_plus_ - n_minus_) * grad(psi_)
            #    - grad(problem.d) * (e_para * dot(problem.d, grad(psi_)) * grad(psi_))
            #ns_force = e_para * dot(grad(psi_), problem.d) * dot(grad(psi_), grad(problem.d))
            u, p = ns_solver.solve(psi = psi_, d = problem.d, f = ns_force)

            #Calculate the L^2 error of the solution to the tentative estimate for the velocity u  
            error_u = errornorm(u, u_, 'L2')

            #Update tentative solution
            u_.assign(u)
            p_.assign(p)

            #Solve Poisson equation
            psi = p_solver.solve(problem.B * (n_plus_ - n_minus_), d = problem.d)

            #Calculte the H^1 error of the solution to the tentative estimate for the electric field psi
            error_psi = errornorm(psi, psi_, 'H10')

            #Update tentative solution
            psi_.assign(psi)

            #Solve Nernst--Planck equation
            n_plus, n_minus = np_solver.solve(u = u_, psi = psi_, d = problem.d)

            #Calculte the L-infinity error of the solution to the tentative estimate for the charges
            error_plus = norm(n_plus.vector() - n_plus_.vector(), 'linf')
            error_minus = norm(n_minus.vector() - n_minus_.vector(), 'linf')

            #Calculate the error for the fixed point iteration
            error = error_u + error_psi + error_plus + error_minus
            #print(error)

            #Update tentative solution
            n_plus_.assign(n_plus)
            n_minus_.assign(n_minus)
            

        #Set solution at next time step to tentative solution
        problem.n_plus_i.assign(n_plus_)
        problem.n_minus_i.assign(n_minus_)
        problem.u_i.assign(u_)
        problem.p_i.assign(p_)
        problem.psi_i.assign(psi_)

        #Save solutions to vtk file every freq time step
        if (i/freq).is_integer():
            write_to_vtk(vtk_files, problem.functions, t)

        #Calculate the cfl number
        cfl_numbers.append(cfl_number(u, dim, dt, n))

    #Save the velocity plot
    save_plot(path, 'velocity', u)

    #Save the electric field plot
    save_plot(path, 'electric_field', grad(psi))

    #Save the CFL-number plot
    save_plot(path, 'cfl', cfl_numbers, x_axis=times)
    
    #Send message via telegram bot
    telegram_send.send(messages=["Hurray!! The simmulation is done!"])

    with open(path + "/plots/velocity.png", "rb") as f:
        telegram_send.send(images=[f])

    return path