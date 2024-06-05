#!/usr/bin/env python 

class base_problem:
    
    def __init__(self, name, dt, n, functions, bcs, ics, spaces, t = 0):
        self.name = name
        self.dt = dt
        self.t = t
        self.n = n
        self.dim = dim
        self.functions = functions
        self.bcs = bcs
        self.ics
        self.spaces
    
    def time_update(self, t):
        self.t = t