from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

import sys
sys.path.insert(0,"..")

class chorin_projection:

    def __init__(self, problem):
        self.V = problem.V
        self.Q = problem.Q
        self.bc_u = problem.bc_u
        self.dt = problem.dt
        self.u_i = problem.u_i
        self.p_i = problem.p_i
        self.rey = problem.rey
        self.e_para = problem.e_para

    def solve(self, psi, d, f = Constant((0,0))):

        #Define trial- and testfunctions
        self.u = TrialFunction(self.V)
        self.p = TrialFunction(self.Q)
        self.v = TestFunction(self.V)
        self.q = TestFunction(self.Q)

        u_ = Function(self.V)
        p_ = Function(self.Q)

        #Define variational form for the tentative velocity
        F_tent = dot(self.u - self.u_i, self.v) * dx\
            + self.dt * (1./self.rey) * inner(grad(0.5 * (self.u + self.u_i)), grad(self.v)) * dx\
            + self.dt * dot(dot(self.u, nabla_grad(self.u_i)),  self.v) * dx\
            - self.dt * dot(f, self.v) * dx
            #+ self.dt * dot(grad(psi), grad(self.v) * grad(psi)) * dx\
            #+ self.dt * self.e_para * dot(grad(psi), d) * dot(grad(psi), (grad(self.v) * d)) * dx
        
        a_tent, L_tent = lhs(F_tent), rhs(F_tent)
        
        #Define variational form for the pressure update
        a_press = dot(grad(self.p), grad(self.q))*dx
        L_press = - (1./self.dt) * div(u_) *  self.q * dx

        #Define variational form for the velocity
        a = dot(self.u, self.v) * dx
        L = dot(u_, self.v) * dx - self.dt * dot(grad(p_), self.v) * dx

        self.u = Function(self.V)
        self.p = Function(self.Q)

        solve(a_tent == L_tent, self.u, self.bc_u)
        u_.assign(self.u)
        
        solve(a_press == L_press, self.p)
        p_.assign(self.p)

        lb = assemble(p_ * dx)
        vec = np.subtract(np.array(p_.vector().get_local()), lb)
        self.p.vector()[:] = vec
        
        solve(a == L, self.u, self.bc_u)
        
        return (self.u, self.p)
    
