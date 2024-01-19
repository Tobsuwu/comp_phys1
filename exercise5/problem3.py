""" Matrix form of 1D Poisson equation """
""" See gitlab for the poisson1D_definition.pdf """
""" The results don't quite match the given figure with n = 11 but as the n grows the results get better. """

# Import the needed packages
import numpy as np
import matplotlib.pyplot as plt


def poisson1d(n):
    """
    This function calculates the solution of the poisson equation numerically using the given Dirichlet boundary
    conditions: phi(0) = phi(1) = 0. It initializes the matrices needed to solve the linear equation Az = b.
    The derivation for this is carried out in the poisson1D_definition_pdf available in the course repo.
    This function also calls the test function.

    :param n: The number of grid points we wish to use. The numerical approximation becomes more accurate when n grows.
    :return: None
    """

    h = 1/(n-1)  # Spacing
    x = np.linspace(0, 1, n)  # Creating the grid using the boundary conditions and the number of points

    # Creating the coefficient matrix for the second derivative approximation
    off_diag = np.ones((n - 1,))
    K = 2*np.identity(n) - np.diag(off_diag,-1) - np.diag(off_diag, 1)
    K[0, 0] = 1  # Implementing the boundary conditions to the matrix.
    K[0, 1] = 0
    K[-1, -1] = 1
    K[-1, -2] = 0

    F = np.zeros((n, 1))
    for j in range(n):
        F[j] = np.pi**2*np.sin(np.pi*x[j])  # Assigning the grid values to the test function

    num_phi = np.linalg.solve(1/(h**2)*K, F)  # Solving the linear matrix equation

    test_poisson1d(num_phi, x)


def test_poisson1d(num_phi, grid):
    """
    This function plots the numerical solution of the Poisson equation in respect to the uniform grid and compares it
    to the known analytical solution of the equation over the same grid. It also calculates the mean absolute error
    (MAE) for the grid.

    :param num_phi: The numerical solution for the Poisson equation calculated by the function poisson1D.
    :param grid: the grid
    :return: None
    """

    ana_phi = np.sin(np.pi*grid)  # The analytical solution
    num_points = str(len(grid))

    plt.figure()
    plt.plot(grid, ana_phi, '--', grid, num_phi)
    plt.legend(["Analytical solution", "Finite difference: N = " + num_points])
    plt.xlabel("x")
    plt.ylabel("phi(x)")
    plt.suptitle("1D Poisson equation")
    plt.show()

    # Calculating the error
    error = np.zeros_like(grid)
    for i in range(len(grid)):
        error[i] = abs(ana_phi[i] - num_phi[i])
    MAE = np.mean(error)

    print("Mean absolute error of the grid with N = " + num_points + ": ", MAE)


def main():
    """
    Main function that controls the number of points given to the calculation.
    :return:
    """
    n = [5, 11, 21, 40]
    for i in range(len(n)):
        poisson1d(n[i])


main()