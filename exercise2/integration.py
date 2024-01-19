""" This program calculates integrals using the trapz, simpson, quad and nquad integrations
from the scipy package """

import numpy as np
import scipy.integrate as si

def trapezoid(function, x):
    # Calculates the integral of the given function using the trapezoidal rule.

    return si.trapz(function, x)

def simpson(function, x):
    # Calculates the integral of the given function using the simpsons rule.

    return si.simpson(function, x, even="avg")

def function_a1(x):
    # Defines the function which integral we wish to calculate.
    return x**2*np.exp(-2*x)*np.sin(2*x)

def function_a2(x):
    # Defines the function which integral we wish to calculate.
    return np.sin(x)/x

def function_a3(x):
    # Defines the function which integral we wish to calculate.
    return np.exp(np.sin(x*x*x))

def function_b(x, y):
    # Defines the function which integral we wish to calculate.
    return x*np.exp(-np.sqrt(x*x+y*y))

def function_c(x, y, z):
    return (1/np.sqrt(np.pi)*(np.exp(-np.sqrt((x+1)**2 + y**2 + z**2)) + np.exp(-np.sqrt((x-1)**2 + y**2 + z**2))))**2

def print_results_a(function, x):
    # Prints the results of the trapezoid and simpson integrations for 1a

    print("{0: <20}".format("trapezoid"), "{0: <20}".format("simpson"), "{0: <20}".format("num_points"))
    for i in range(len(x)):
        print("{0: <20}".format(trapezoid(function(x[i]), x[i])),
              "{0: <20}".format(simpson(function(x[i]), x[i])), "{0: <20}".format(len(x[i])))

def print_results_b():
    # Prints the results of the trapezoid and simpson integrations for 1b

    x = [np.linspace(0, 2, 10), np.linspace(0, 2, 20), np.linspace(0, 2, 40), np.linspace(0, 2, 100),
           np.linspace(0, 2, 200), np.linspace(0, 2, 400)]
    y = [np.linspace(-2, 2, 10, ), np.linspace(-2, 2, 20, ), np.linspace(-2, 2, 40, ), np.linspace(-2, 2, 100, ),
           np.linspace(-2, 2, 200, ), np.linspace(-2, 2, 400, )]

    print("{0: <20}".format("trapezoid"), "{0: <20}".format("simpson"), "{0: <20}".format("num_points"))
    for i in range(len(x)):
        print("{0: <20}".format(trapezoid(trapezoid(function_b(x[i][:, None], y[i][None, :]), y[i]), x[i])),
              "{0: <20}".format(simpson(simpson(function_b(x[i][:, None], y[i][None, :]), y[i]), x[i])), "{0: <20}".format(len(x[i])))

def print_results_c():
    # Prints the results of the trapezoid and simpson integrations for 1c
    # Generate the points using linspace

    x = [np.linspace(-10, 10, 10), np.linspace(-10, 10, 20), np.linspace(-10, 10, 40), np.linspace(-10, 10, 100),
           np.linspace(-10, 10, 200), np.linspace(-10, 10, 400)]
    y = x
    z = y
    print("{0: <20}".format("trapezoid"), "{0: <20}".format("simpson"), "{0: <20}".format("num_points"))
    for i in range(len(x)):
        print("{0: <20}".format(trapezoid(trapezoid(trapezoid(function_c(x[i][:, None, None], y[i][None, :, None], z[i][None, None, :]), x[i]), y[i]), z[i])),
              "{0: <20}".format(simpson(simpson(simpson(function_c(x[i][:, None, None], y[i][None, :, None], z[i][None, None, :]), x[i]), y[i]), z[i])), "{0: <20}".format(len(x[i])))


def main():
    # Creates the integration limits into lists for problem 1a and runs the other functions.

    x_a1 = [np.linspace(0, 10, 10), np.linspace(0, 10, 20), np.linspace(0, 10, 40), np.linspace(0, 10, 100), np.linspace(0, 10, 200), np.linspace(0, 10, 400)]
    x_a2 = [np.linspace(0.0001, 1, 10), np.linspace(0.0001, 1, 20), np.linspace(0.0001, 1, 40),
            np.linspace(0.0001, 1, 100), np.linspace(0.0001, 1, 200), np.linspace(0.0001, 1, 400)]
    x_a3 = [np.linspace(0, 5, 10), np.linspace(0, 5, 20), np.linspace(0, 5, 40), np.linspace(0, 5, 100),
            np.linspace(0, 5, 200), np.linspace(0, 5, 400)]

    print("1a.1")
    print_results_a(function_a1, x_a1)
    print("quad = ", si.quad(function_a1, 0, np.inf))
    print()
    print("1a.2")
    print_results_a(function_a2, x_a2)
    print("quad = ", si.quad(function_a2, 0, 1))
    print()
    print("1a.3")
    print_results_a(function_a3, x_a3)
    print("quad = ", si.quad(function_a3, 0, 5))
    print()

    print("1b")
    print_results_b()
    print("nquad= ", si.nquad(function_b, [[0, 2], [-2, 2]]))
    print()

    print("1c")
    print_results_c()
    print("nquad= ", si.nquad(function_c, [[-10, 10], [-10, 10], [-10, 10]]))

if __name__=="__main__":

    main()