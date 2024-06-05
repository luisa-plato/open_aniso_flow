#!/usr/bin/env python

import numpy as np
from fenics import *

tol = 1E-14

class ics:
    
    def __init__(self, name, space, value = 0.5):
        self.name = name
        self.space = space
        if name == "zero":
            self.ic = self.zero()
        elif name == "vortex":
            self.ic = self.vortex()
        elif name == "left_gaussian":
            self.ic = self.left_gaussian()
        elif name == "right_gaussian":
            self.ic = self.right_gaussian()
        elif name == "constant":
            self.ic = self.constant(value)
        elif name == "left_half":
            self.ic = self.left_half()
        elif name == "right_half":
            self.ic = self.right_half()
        elif name == "right_slab":
            self.ic = self.right_slab()
        elif name == "small_square":
            self.ic = self.small_square()

    def zero(self):
        d = self.space._ufl_element.reference_value_size()
        if d == 1:
            return project(Constant(0), self.space)
        else:
            return project(Constant(tuple(np.zeros(d))), self.space)

    def vortex(self):
        return project(Expression((' sin(pi * x[0]) * cos(pi * x[1])', '-sin(pi * x[1]) * cos(pi * x[0])'), degree = 2, pi = np.pi), self.space)

    def left_gaussian(self):
        return project(Expression('10*exp(-40*pow(x[0]-0.25,2)-40*pow(x[1]-0.5,2))', degree = 1, tol=tol), self.space)

    def right_gaussian(self):
        return project(Expression('10*exp(-40*pow(x[0]-1.25,2)-40*pow(x[1]-0.5,2))', degree = 1, tol=tol), self.space)

    def constant(self, value):
        return project(Constant(value), self.space)

    def left_half(self):
        return project(Expression('x[0] < 0.75 + tol? 1 : tol', degree = 2, tol = tol), self.space)
    
    def right_half(self):
        return project(Expression('x[0] > 0.75 + tol? 1 : tol', degree = 2, tol = tol), self.space)
    
    def right_slab(self):
        return project(Expression('x[0] > 1 + tol && (x[1] > 0.4 - tol && x[1] < 0.6 + tol)? 1 : tol', degree = 2, tol = tol), self.space)
    
    def small_square(self):
        return project(Expression('(x[0] > 0.5 - tol && x[0] < 1 + tol) && (x[1] > 0.3 - tol && x[1] < 0.7 + tol)? 1 : tol', degree = 2, tol = tol), self.space)

