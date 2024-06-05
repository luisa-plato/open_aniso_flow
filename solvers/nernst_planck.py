from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

import sys
sys.path.insert(0,"..")

class nernst_planck:

    def __init__(self, problem):
        self.Q = problem.Q
        self.n_plus_i = problem.n_plus_i
        self.n_minus_i = problem.n_minus_i
        self.dt = problem.dt
        self.mu = problem.mu
        self.e_orth = problem.e_orth
        self.e_para = problem.e_para
        self.lamb_para = problem.lamb_para
        self.lamb_orth = problem.lamb_orth
        self.F = problem.F

    def solve(self, u, psi, d = Constant((0,0))):

        #Define trial- and testfunctions
        self.n_plus = TrialFunction(self.Q)
        self.n_minus = TrialFunction(self.Q)
        self.q_plus = TestFunction(self.Q)
        self.q_minus = TestFunction(self.Q)

        #Define linearized variational forms for the charges
        a_plus = self.n_plus * self.q_plus * dx\
            + 0.5 * self.dt * self.mu * self.lamb_orth * dot(grad(self.n_plus),grad(self.q_plus)) * dx\
            + 0.5 * self.dt * self.mu * self.lamb_para * dot(grad(self.n_plus),d) * dot(d, grad(self.q_plus)) * dx\
            + self.dt * self.F * self.lamb_orth * self.n_plus * dot(grad(psi), grad(self.q_plus)) * dx\
            + self.dt * self.F * self.lamb_para * self.n_plus * dot(grad(psi), d) * dot(d, grad(self.q_plus)) * dx\
            - self.dt * self.n_plus * dot(u, grad(self.q_plus)) * dx

        a_minus = self.n_minus * self.q_minus * dx\
	        + 0.5 * self.dt * self.mu * self.lamb_orth * dot(grad(self.n_minus),grad(self.q_minus)) * dx\
            + 0.5 * self.dt * self.mu * self.lamb_para * dot(grad(self.n_minus),d) * dot(d, grad(self.q_minus)) * dx\
	        - self.dt * self.F * self.lamb_orth * self.n_minus * dot(grad(psi), grad(self.q_minus)) * dx\
            - self.dt * self.F * self.lamb_para * self.n_minus * dot(grad(psi), d) * dot(d, grad(self.q_minus)) * dx\
            - self.dt  * self.n_minus * dot(u, grad(self.q_minus)) * dx

        L_plus = self.n_plus_i * self.q_plus * dx\
            - 0.5 * self.dt * self.mu * self.lamb_orth * dot(grad(self.n_plus_i),grad(self.q_plus)) * dx\
            - 0.5 * self.dt * self.mu * self.lamb_para * dot(grad(self.n_plus_i),d) * dot(d, grad(self.q_plus)) * dx

        L_minus = self.n_minus_i * self.q_minus * dx\
            - 0.5 * self.dt * self.mu * self.lamb_orth * dot(grad(self.n_minus_i),grad(self.q_minus)) * dx\
            - 0.5 * self.dt * self.mu * self.lamb_para * dot(grad(self.n_minus_i),d) * dot(d, grad(self.q_minus)) * dx

        self.n_plus = Function(self.Q)
        self.n_minus = Function(self.Q)

        solve(a_plus == L_plus, self.n_plus)
        solve(a_minus == L_minus, self.n_minus)
        return (self.n_plus, self.n_minus)
    

