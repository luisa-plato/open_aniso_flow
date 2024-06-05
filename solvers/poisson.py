from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

import sys
sys.path.insert(0,"..")

class poisson:

    def __init__(self, problem):
        self.Q = problem.Q
        self.psi_i = problem.psi_i
        self.mu = problem.mu
        self.e_orth = problem.e_orth
        self.e_para = problem.e_para
        self.tau = problem.tau
        self.xi = problem.xi

    def solve(self, f, d = Constant((0,0))):

        #Define trial- and testfunctions
        self.psi = TrialFunction(self.Q)
        self.phi = TestFunction(self.Q)

        #Define the variatonal form
        a = self.e_orth * dot(grad(self.psi), grad(self.phi)) * dx\
            + self.e_para * dot(d,grad(self.psi)) * dot(d, grad(self.phi)) * dx\
            + self.tau * self.psi * self.phi * ds
        
        L = f * self.phi * dx\
            + self.xi * self.phi * ds

        self.psi = Function(self.Q)

        solve(a == L, self.psi)
        return (self.psi)