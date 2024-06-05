#Different director fields as expressions

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

set_log_level(50)
tol = 1E-14

def create_director(name):
	if name == "source_and_saddle":
		return source_and_saddle()
	elif name == "zero":
		return Constant((0,0))
	elif name == "s_curve":
		return s_curve()
	elif name == "curly_bracket":
		return curly_bracket()
	elif name == "epsilon":
		return epsilon()
	elif name == "wind_around_egg":
		return wind_around_egg()
	elif name == "broom":
		return broom()
	elif name == "one_circle":
		return one_circle()
	elif name == "multiple_circle":
		return multiple_circle()
	elif name == "matrix_circle":
		return matrix_circle()
	elif name == "one_source":
		return one_source()
	else:
		print("Raise an error!")

#----------------------------------------------------- flat land -----------------------------------------------------------------
class director(UserExpression):

	def eval(self, values, x):
		values[0] = 1/sqrt(2)
		values[1] = 1/sqrt(2)

	def value_shape(self):
		return (2,)

class source_and_saddle(UserExpression):

	def eval(self, values, x):
		if ((x[0] == 0.146477 or x[0] == 0.853553) and x[1] == 0.5):
			values[0] = 0
			values[1] = 0
		else:
			d = [4*(x[0]-0.5)**2 + 4*(x[1]-0.5)**2-0.5, 2*(x[1]-0.5)]
			norm = sqrt(d[0]**2 + d[1]**2)
			values[0] = d[0]/norm
			values[1] = d[1]/norm

	def value_shape(self):
		return (2,)

class curly_bracket(UserExpression):

	def eval(self, values, x):
		T = np.arcsin(sin(4*pi*x[1]))
		alpha = T
		values[0] = cos(alpha)
		values[1] = sin(alpha)

	def value_shape(self):
		return (2,)

class s_curve(UserExpression):

	def eval(self, values, x):
		T = np.arcsin(sin(4*pi*x[1]))
		alpha = pi/2 + T
		values[0] = cos(alpha)
		values[1] = sin(alpha)

	def value_shape(self):
		return (2,)

class epsilon(UserExpression):

	def eval(self, values, x):
		T = np.arcsin(sin(4*pi*x[1]))
		alpha = pi*(1-x[1]*4)
		values[0] = cos(alpha)
		values[1] = sin(alpha)

	def value_shape(self):
		return (2,)

class broom(UserExpression):

	def eval(self, values, x):
		alpha = 0.5 * np.arctan2(x[1]-0.5,x[0]-0.2) - 0.5 * np.arctan2(x[1]-0.5,x[0]-0.8) + pi/2
		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)

class wind_around_egg(UserExpression):

	def eval(self, values, x):
		alpha = - 0.5 * np.arctan2(x[1]-0.5,x[0]-0.7) + 0.5 * np.arctan2(x[1]-0.5,x[0]-0.8)
		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)

class one_circle(UserExpression):

	def eval(self, values, x):
		alpha = - 0.5 * np.arctan2(x[1]-0.5,x[0]-0.65) + np.arctan2(x[1]-0.5,x[0]-0.75) - 0.5 * np.arctan2(x[1]-0.5,x[0]-0.85) 
		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)

class multiple_circle(UserExpression):

	def eval(self, values, x):
		d = 0.1
		alpha = 0
		for j in range(3):
			y = np.sqrt(3)*j*d + 0.2
			alpha += - 0.5 * np.arctan2(x[1]-y,x[0]) + np.arctan2(x[1]-y,x[0]-0.2) - 0.5 * np.arctan2(x[1]-y,x[0]-0.4)\
			- 0.5 * np.arctan2(x[1]-y,x[0]-0.6) + np.arctan2(x[1]-y,x[0]-0.8) - 0.5 * np.arctan2(x[1]-y,x[0]-1)
		for i in range(4):
			y = np.sqrt(3)*i*d + (1-np.sqrt(3)/2)*0.2
			alpha += - 0.5 * np.arctan2(x[1]-y,x[0]+0.3) + np.arctan2(x[1]-y,x[0]+0.1) - 0.5 * np.arctan2(x[1]-y,x[0]-0.1)\
				- 0.5 * np.arctan2(x[1]-y,x[0]-0.3) + np.arctan2(x[1]-y,x[0]-0.5) - 0.5 * np.arctan2(x[1]-y,x[0]-0.7)\
				- 0.5 * np.arctan2(x[1]-y,x[0]-0.9) + np.arctan2(x[1]-y,x[0]-1.1) - 0.5 * np.arctan2(x[1]-y,x[0]-1.3)

		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)

class matrix_circle(UserExpression):

	def __init__(self, m=4, **kwargs):
		self.m = m
		self.n = m
		self.d = 1/(m+1)
		super().__init__()

	def eval(self, values, x):
		d = self.d
		alpha = 0

		m = self.m
		n = self.n
		for i in range(n+1):
			for j in range(m+1):
				#print(i,j)
				d_m = -np.sqrt(3) * i * d
				d_n = -3 * j * d
				alpha += np.arctan2(x[1]+d_m, x[0]+d_n)\
					- 0.5 * np.arctan2(x[1]+d_m, x[0]+d_n+d)\
					- 0.5 * np.arctan2(x[1]+d_m, x[0]+d_n-d)

		for i in range(n+2):
			for j in range(m+2):
				d_p = -0.5 * np.sqrt(3) * (2*i+1) * d
				d_q = -0.5 * 3 * (2*j+1) * d
				alpha += np.arctan2(x[1]+d_p, x[0]+d_q)\
					- 0.5 * np.arctan2(x[1]+d_p, x[0]+d_q+d)\
					- 0.5 * np.arctan2(x[1]+d_p, x[0]+d_q-d)

		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)

class one_source(UserExpression):

	def eval(self, values, x):
		alpha = 0.5 * np.arctan2(x[1]-0.5,x[0] - 0.65) - np.arctan2(x[1]-0.5,x[0]-0.75) + 0.5 * np.arctan2(x[1]-0.5,x[0]-0.85) 
		values[0] = cos(alpha)
		values[1] = sin(alpha) 

	def value_shape(self):
		return (2,)