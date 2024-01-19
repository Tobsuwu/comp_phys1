""" This program calculates the steepest descent of the function """

import numpy as np

def f(x, N):
    # Create an N-dimensional one element vector, where x is the initial guess for the solution
    g = np.zeros((N, 1))
    for i in range(len(g)):
        g[[i, 1]] = x**2
    return g

def derivative_1D(g, x, dx):
     # Calculates the 1D derivative of the given row of N-dimensional vector

     return (g(x + dx) - g(x - dx)) / (2 * dx)

def gradient_ND(g, x, N, dx):
    G = np.zeros((N, 1))
    for i in range(len(g)):
        G[i] = derivative_1D(g[i], x, dx)
    return G


def steepest_descent(a, x, g, N):
    # Calculates the steepest descent by using the equations (31) and (32) of the lecture slides. If x is a
    # constant vector, its derivative is 0, so the gamma = a.

    x_n = x
    gamma = a
    x_new = x_n - gamma*gradient_ND(g, x, N, dx=0.01)
    return x_new



def main():
    # Making the initial guess. Dimensions of guess need to match the dimension given to function f.
    x = np.array([[1], [1], [1], [1], [1], [1], [1], [1], [1], [1]])
    print(steepest_descent(0.001, x, f(x, 10), 10))

if __name__ == '__main__':
    main()