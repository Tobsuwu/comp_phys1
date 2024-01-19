""" This program calculates a given 2D-integral using scipys Simpson rule and nquad. """

# Importing the needed packages.
import numpy as np
from scipy.integrate import simps, nquad

def f(x, y):
    # Defining the 2D integral which we will integrate

    return (x+y)*np.exp(-1/4*np.sqrt(x**2+y**2))

def calculate_simpson(x0, x1, y0, y1, n):
    # Calculates the 2D integral of the function f. The number of grid points depends on parameter n.
    # In case of even number of points, we use the parameter even="avg" to calculate the average between calculations
    # 1) and 2).
    # 1) calculates the first n-2 intervals with simpsons rule and uses trapezoidal rule for the last one.
    # 2) calculates the last n-2 intervals with simpsons rule and uses trapezoidal rule for the first one.

    x = np.linspace(x0, x1, n)
    y = np.linspace(y0, y1, n)
    integral = simps(simps(f(x[:, None], y[None, :]), x, even="avg"), y, even="avg")
    return integral

def main():
    # Creates a list containing the number of points we wish to test our simpson integral and compares it to the nquad
    # function

    n = [5, 10, 15, 20, 25, 30, 100, 1000]
    print("{0: <20}".format("Simpson integral"), "{0: <20}".format("Number of points"))
    for i in range(len(n)):
        print("{0: <20}".format(calculate_simpson(0, 2, -2, 2, n[i])), "{0: <20}".format((n[i])))
    print("nquad: ", nquad(f, [[0, 2], [-2, 2]]))


main()