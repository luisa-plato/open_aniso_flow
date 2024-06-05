
import argparse
import numpy as np

import sys
sys.path.insert(0,"..")
from run_experiment import *

parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--dim", metavar="dimension", type=int, default=2)
parser.add_argument("--meshsize", type=int, default = 32)
parser.add_argument("--dt", type=float, default = 0.01)
parser.add_argument("--T", type=float, default = 1.0)
parser.add_argument("--problem", type=str, default="no_slip_box")
parser.add_argument("--coupling", type=str, default="sequential", help="choose between a sequential coupling or a fixed point solver")
parser.add_argument("--ns_solver", type=str, default="chorin_projection")
parser.add_argument("--mu", type=float, default = 0.004, help="specify the viscosity for the charges")
parser.add_argument("--rey", type=float, default= 0.001, help="specify the Reynolds number for the velocity field")
parser.add_argument("--max_iter", type=int, default="6")
parser.add_argument("--fp_tol", type=float, default=1.0e-6)
parser.add_argument("--save_freq", type=int, default=2)
parser.add_argument("--director", type=str, default = "wind_around_egg")
parser.add_argument("--u_i", type=str, default = "zero")
parser.add_argument("--p_i", type=str, default = "zero")
parser.add_argument("--psi_i", type=str, default = "zero")
parser.add_argument("--n_plus_i", type=str, default = "constant")
parser.add_argument("--n_minus_i", type=str, default = "constant")
parser.add_argument("--bc_psi_amplitude", type=float, default = 1)
parser.add_argument("--psi_freq", type=float, default = 5)
parser.add_argument("--e_orth", type=float, default = 1)
parser.add_argument("--e_para", type=float, default = 0.2)
parser.add_argument("--lamb_para", type=float, default = 0.42)
parser.add_argument("--lamb_orth", type=float, default = 1)
parser.add_argument("--tau", type=float, default = 1)
parser.add_argument("--B", type=float, default= 34.75)
parser.add_argument("--F", type=float, default= 2.8)


args = parser.parse_args()
parameters = vars(args)

path = run_experiment(parameters)
print(path)