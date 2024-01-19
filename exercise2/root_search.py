""" This program searches the root of the given function. """

import numpy as np

def f(x):
    # Defines the function

    return x**2-4

def linear_grid(a, b, num_points):
    # Returns a linear grid

    return np.linspace(a, b, num_points)

def custom_grid(r_0, h, dim):
    # Returns a custom grid with given requirements

    r = [0]*dim
    for i in range(1, dim):
        r[i] = r_0*(np.exp(i*h)-1)
    return r

def secant(f, a_0, b_0, N):
    # Approximates the root by using the secant method. Requires a starting guess for the interval (grid start
    # and end points). The starting guess has to have different sign function values in endpoints.
    # Parameter N determines the max amount of iterations.

    if f(a_0)*f(b_0) >= 0:
        # Tests the starting guess
        print("Try different guess :)")
        return None
    a_n = a_0
    b_n = b_0
    for i in range(1, N+1):
        # Calculates the root using the secant algorithm.

        x_n = a_n - (b_n - a_n)*f(a_n)/(f(b_n)-f(a_n))
        if f(a_n)*f(x_n) < 0:
            b_n = x_n
        elif f(b_n)*f(x_n) < 0:
            a_n = x_n
        elif f(x_n) == 0:
            return x_n
    return a_n - (b_n - a_n)*f(a_n)/(f(b_n)-f(a_n))


def main():
    r_0 = 10**(-5)
    r_max = 100
    dim = 100
    h = np.log(r_max/r_0+1)/(dim-1)
    x = linear_grid(-1, 3, 100)
    root = secant(f, x[0], x[-1], 5)
    print(root)


if __name__=="__main__":

    main()