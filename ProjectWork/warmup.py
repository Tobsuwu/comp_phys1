""" This program calculates a given 2D-integral using scipys Simpson rule and plots the real space domain
in a contourf plot. """

# Importing the needed packages.
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

def fun(x, y):
    """
    This function returns the 2D integrand we wish to integrate with simps.

    :param x: value of x-coordinate
    :param y: value of y-coordinate
    :return: integrand
    """

    return (x+3*y)*np.exp(-1/2*np.sqrt(x**2+y**2))

def calculate_2d_simpson(n):
    """
    This function calculates the 2D integral with scipys simpson rule. It also plots the contourf plot over the
    domain omega defined by the calculated real space vectors.

    :param n: number of points used in integration. n*n is also the grid size.
    :return: None
    """

    # "Lattice" vectors
    a = np.array([1.3, 0.2])
    b = np.array([0.45, 1])
    A = np.array([a, b]).transpose()
    jacobian = np.linalg.det(A)  # Linear transform. Jacobian is the determinant of matrix A.

    # Create the alpha space vectors.
    alpha = np.linspace(0, 1, n)

    # Initialize the real space vectors and the integrand.
    x = np.zeros_like(alpha)
    y = np.zeros_like(alpha)
    integrand = np.zeros([n, n])

    for i in range(n):
        for j in range(n):
            x[i], y[j] = A.dot(np.array([alpha[i], alpha[j]]))  # real space r = A dot alpha
            integrand[i, j] = fun(x[i], y[j])

    integral = simps(simps(integrand, alpha, even="avg"), alpha, even="avg")*jacobian
    print("Value of the integral:", integral)

    # Plotting the contourf plot
    fig, ax = plt.subplots(1, 1)
    X, Y = np.meshgrid(x, y)  # Over the real space. e.g. domain Omega
    cb = ax.contourf(X.T, Y.T, integrand)
    fig.colorbar(cb)
    ax.set_title('2D-integral over the domain Omega')
    ax.set_xlabel('x, real space coordinates')
    ax.set_ylabel('y, real space coordinates')
    plt.show()

    return None


def main():
    """
    Calls the integration function and determinbes how many grid points are used.

    :return:
    """

    n = 101
    calculate_2d_simpson(n)

main()