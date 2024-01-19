""" Poisson equation and relaxation methods. """

import numpy as np
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D


def jacobi_method(phi, rho, N, h):
    """
    This function implements the relaxation algorithm using the Jacobi method.

    :param phi: the approximated solution for the PDE we wish to iterate
    :param rho: charge density, zero in this exercise
    :param N: number of grid points
    :param h: grid spacing
    :return: phi_new: The new value of our PDE solution
    """

    phi_new = np.copy(phi)
    for i in range(1, N-1):
        for j in range(1, N-1):
            phi_new[i, j] = 0.25*(phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1] - rho[i, j]*h**2)

    return phi_new


def gauss_seidel_method(phi, rho, N, h):
    """
    This function implements the relaxation algorithm using the Gauss-Seidel method.

    :param phi: the approximated solution for the PDE we wish to iterate
    :param rho: charge density, zero in this exercise
    :param N: number of grid points
    :param h: grid spacing
    :return: phi_new The new value of our PDE solution
    """

    phi_new = np.copy(phi)
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            phi_new[i, j] = 0.25*(phi[i+1, j] + phi_new[i-1, j] + phi[i, j+1] + phi_new[i, j-1] - rho[i, j]*h**2)

    return phi_new


def SOR_method(phi, rho, N, h, omega):
    """

    This function implements the relaxation algorithm using the SOR method.

    :param phi: the approximated solution for the PDE we wish to iterate
    :param rho: charge density, zero in this exercise
    :param N: number of grid points
    :param h: grid spacing
    :return: phi_new The new value of our PDE solution
    """

    phi_new = np.copy(phi)
    for i in range(1, N-1):
        for j in range(1, N-1):
            phi_new[i, j] = (1-omega)*phi[i, j] + omega/4*(phi[i+1, j] + phi_new[i-1, j] + phi[i, j+1] + phi_new[i, j-1]
                                                           -rho[i, j]*h**2)

    return phi_new


def main():
    """
    Main function that initializes the grids and boundary conditions for the solution given in the problem.
    Executes the while loops that call the relaxation algorithms and plots the results.

    :return:
    """
    L = 1  # Boundary condition. Max value of x and y.
    N = 21  # Number of points in the grid.
    h = L/(N-1)  # Grid spacing.
    x = y = np.linspace(0, L, N)  # Initializing the uniform grids.
    X, Y = np.meshgrid(x, y)  # Meshgrid for plotting.
    rho = np.zeros((N,N))  # Charge density, zero in this problem
    phi = np.zeros((N,N))  # Initialization of our PDE solution

    phi[0, :] = 1  # Boundary conditions from the problem setup. Initial value of the solution.

    # Initializing the solutions computed by different methods
    phi_jacobi = np.copy(phi)
    phi_GS = np.copy(phi)
    phi_SOR = np.copy(phi)

    error = 1000  # Initializing the error and determining the tolerance level.
    tol = 1.0e-6
    num_iters_jacobi = 0  # Calculating the number of iterations it takes to reach the desired tolerance.
    while error > tol:
        phi_old = np.copy(phi_jacobi)
        phi_jacobi = jacobi_method(phi_jacobi, rho, N, h)
        error = abs(np.mean(phi_jacobi)-np.mean(phi_old))
        num_iters_jacobi += 1

    error = 1000
    num_iters_GS = 0
    while error > tol:
        phi_old = np.copy(phi_GS)
        phi_GS = gauss_seidel_method(phi_GS, rho, N, h)
        error = abs(np.mean(phi_GS)-np.mean(phi_old))
        num_iters_GS += 1

    error = 1000
    num_iters_SOR = 0
    omega = 1.8
    while error > tol:
        phi_old = np.copy(phi_SOR)
        phi_SOR = SOR_method(phi_SOR, rho, N, h, omega)
        error = abs(np.mean(phi_SOR)-np.mean(phi_old))
        num_iters_SOR += 1

    #  Plotting the results in a figure alongside the number of required iterations.
    fig = figure()

    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot_wireframe(X, Y, phi, rstride=1, cstride=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Initial solution')

    ax2 = fig.add_subplot(222, projection='3d')
    ax2.plot_wireframe(X, Y, phi_jacobi, rstride=1, cstride=1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title('Solution solved by Jacobi method')
    ax2.text(1, 0.1, 1.2, 'Number of iterations: ' + str(num_iters_jacobi))

    ax3 = fig.add_subplot(223, projection='3d')
    ax3.plot_wireframe(X, Y, phi_GS, rstride=1, cstride=1)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_title('Solution solved by Gauss-Sedel method')
    ax3.text(1, 0.1, 1.2, 'Number of iterations: ' + str(num_iters_GS))

    ax4 = fig.add_subplot(224, projection='3d')
    ax4.plot_wireframe(X, Y, phi_SOR, rstride=1, cstride=1)
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_title('Solution solved by SOR method')
    ax4.text(1, 0.1, 1.2, 'Number of iterations: ' + str(num_iters_SOR))
    fig.suptitle("Laplace equation in 2D solved by different relaxation algorithms")

    show()

main()