""" This file contains functions for numerical differentation and integration """

import numpy as np
import sympy as sp

def first_derivative(function, x, dx):

# Using the derivation of the first derivative found in the lecture notes (dx = h)

	fdx = (function(x + dx) - function(x - dx)) / (2*dx)
	return fdx

def test_first_derivative():
	
# Using imported library sympy to symbolically calculate the first derivative.

	x = sp.Symbol('x')
	func = fun(x)
	return sp.diff(func)
	
def second_derivative( function, x, dx):

# Using the derivation for the second derivative found in the lecture notes (dx = h). 

	f2dx = (function(x + dx) + function(x - dx) - 2*function(x)) / (dx**2)
	return f2dx

def test_second_derivative():

# Using imported library sympy to symbolically calculate the second derivative.
	x = sp.Symbol('x')
	func = fun(x)
	return sp.diff(sp.diff(func))
	
def riemann_int(fun, x):

# Calculating left-side riemann sum with uniform grid. Delta_x is constant

	delta_x = abs(x[1] - x[0])
	rie_sum = delta_x * np.sum(fun)
	return rie_sum

def test_analytical_int():

# Testing the numeric integrals functions by comparing them to the analytical integral.

	a = 1
	b = 2
	x = sp.Symbol('x')
	f = fun(x)
	rie_sum = sp.integrate(f, (x, a, b))
	return rie_sum

def trapezoid_int(fun, x):

# Calculating the integral using trapezoidal rule. Uniform grid -> delta_x is constant

	delta_x = abs(x[1] - x[0])
	fun_right = fun[1:]
	fun_left = fun[:-1]
	trapz_int = 1/2*delta_x*np.sum(fun_right + fun_left)
	return trapz_int

def simpson_int(f, a, b, N):

# Calculating the integral using simpson rule. Uniform grid -> delta_x is constant. Different formulas for even and odd amount of subintervals.

# The implementation of the uneven amount of subintervals formula lead into a bad approximation.

	if N % 2 == 1:
		I_even = simpson_int(f, a, b, N-1)
		x = np.linspace(a, b, N)
		f = fun(x)
		h = abs(x[1] - x[0])
		delta_I = h/12*(-f[N-3]+8*f[N-2]+5*f[N-1])
		simp_int = delta_I + I_even
		return simp_int
		
	x = np.linspace(a, b, N+1)
	f = fun(x)
	h = abs(x[1] - x[0])
	simp_int = h/3*np.sum( f[0:-1:2] + 4*f[1::2] + f[2::2])
	return simp_int



def fun(x):

# Defining the function which derivatives and integrals we wish to calculate

	return 3*x**2	



def main():

# Calculating first and second derivative. Numerical approximation is more accurate when dx is small. x indicates the point in which we want to calculate the function value.

	dx = 0.1
	x = 3
	print()
	print("Numerical first derivative:", first_derivative(fun, x, dx), ", x = ", x)
	print("Analytical first derivative:", test_first_derivative())
	print()
	print("Numerical second derivative:", second_derivative(fun, x, dx), ", x = ", x)
	print("Analytical second derivative:", test_second_derivative())
	print()
	
# Calculating the numerical integrals. Numpy linspace to create a vector.
	
	x = np.linspace(1, 2, 100)
	f = fun(x)
	print("Numerical Riemann sum:", riemann_int(f, x))
	print("Numerical trapezoidal sum", trapezoid_int(f, x))
	
	a = 1
	b = 2
	N = 100
	print("Numerical Simpson integral", simpson_int(f, a, b, N))
	print("Analytical integral:", test_analytical_int())

if __name__=="__main__":

	main()
