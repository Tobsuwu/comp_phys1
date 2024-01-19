""" Problem 2: Gauss law and electric flux """
""" This program calculates the electric flux over a side/plane of a cubic box generated by a point charge 
    inside the box. If the charge is considered as a plane inside the box, this program treats the grid points of 
    the plane as point charges and calculates the total flux as a sum of fluxes generated by the point charges over the 
    grid. """
""" Derivations for the equations used in this program are available in the repo. """

# Importing the needed packages.
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt


def integrand(dx, dy, q_pos, side_plane):
    """
    This function Formulates the integrand which the calculate_flux function integrates.

    :param dx: first differential of the area element.
    :param dy: second differential of the area element.
    :param q_pos: position of the point charge inside the box.
    :param side_plane: the plane/side of the box which flux we wish to calculate.
    :return: E: the integrand which we integrate.
    """

    # Charge and permittivity. Constants in this problem
    q = 1
    epsilon = 1

    # Initializing the plane point in order to calculate the directional vector r and the unit normal vector regards
    # to the plane element.

    plane_point = np.array([0, 0, 0])
    n_unit = np.array([0, 0, 0])

    # Formulate the proper unit normal vector of the area element and plane point depending on the given parameter
    # side/plane. Unit normal vector regards to the plane element is always [1, 0, 0], [0, 1 ,0], [0, 0, 1] depending
    # which side/plane flux we wish to calculate.

    if side_plane[0] != 0:
        plane_point = np.array([side_plane[0], dx, dy])
        n_unit = np.array([1, 0, 0])
    elif side_plane[1] != 0:
        plane_point = np.array([dx, side_plane[1], dy])
        n_unit = np.array([0, 1, 0])
    elif side_plane[2] != 0:
        plane_point = np.array([dx, dy, side_plane[2]])
        n_unit = np.array([0, 0, 1])

    r = -q_pos + plane_point  # Directional vector determined by two points of the vector.
    r_unit = r/np.linalg.norm(r)  # Unit directional vector

    # Calculating the point charge in regards to area element normal vector. Pdf in the repo for definition.
    E = 1/(4*np.pi) * q/epsilon * 1/(np.linalg.norm(r)**2) * np.dot(r_unit, n_unit)
    return E


def calculate_point_charge_flux(L, q_pos, side_plane, N):
    """
    This function initializes the area element grid and calls the integrand function on the grid points in order to
    calculate the electric field over it e.g. flux.

    :param q_pos: the position of the point charge inside the LxLxL box as numpy array
    :param side_plane: the plane/side of the box which flux we wish to calculate as numpy array
    :param N: number of points used in the simpson integration and grid size N*N in area element
    :return: flux: Calculated flux value
    """

    # Creating the area element with respect to dx and dy. dA = dx*dy
    dx = np.linspace(-L/2, L/2, N)
    dy = np.linspace(-L/2, L/2, N)

    # Initializing the integrand.
    E = np.zeros([N, N])

    # Calculating the values of integrand.
    for i in range(N):
        for j in range(N):
            E[i, j] = integrand(dx[i], dy[j], q_pos, side_plane)

    # Calculating the flux over area elements.
    flux = simps(simps(E, dx, even='avg'), dy, even='avg')
    return flux

def print_results(q_pos, side_plane, flux):
    """
    This function prints the flux through a certain plane and charge position.

    :param q_pos: charge position inside the box
    :param side_plane: side/plane through which we calculate the flux
    :param flux: the flux
    :return: None
    """
    print("Flux through plane ", side_plane, " when the point charge is at", q_pos)
    print(abs(flux), "(absolute value)")
    return None

def main():
    """
    This function calls the other functions, plots figures and prints some details.
    :return:
    """

    N = 11  # Number of points in simpson integral. Grid N x N in the area element.
    L = 1  # Measurement of the box. L*L*L, constant in this problem.

    q_pos = np.array([0, 0, 0])  # Position of the point charge is in origin
    # Indicates the plane/side in the cube which flux we wish to calculate. e.g. z = L/2
    side_plane = np.array([0, 0, L/2])
    flux = calculate_point_charge_flux(L, q_pos, side_plane, N)
    print_results(q_pos, side_plane, flux)
    print("Which equates to 1/6 of enclosed Q, in line with Gauss's Law.")
    print()

    q_pos = np.array([-L/4, 0, 0])  # Position of point charge (-L/4, 0, 0)
    side_plane = np.array([L/2, 0, 0])  # Flux through plane at (L/2, x, z)
    flux = calculate_point_charge_flux(L, q_pos, side_plane, N)
    print_results(q_pos, side_plane, flux)
    print()

    q_pos = np.array([-L / 4, 0, 0])  # Position of point charge (-L/4, 0, 0)
    side_plane = np.array([0, L / 2, 0])  # Flux through plane (x, L/2, z)
    flux = calculate_point_charge_flux(L, q_pos, side_plane, N)
    print_results(q_pos, side_plane, flux)

    # Calculating the flux through plane (L/2, y, z), when the charge is distributed in the xy-plane. e.g sum of
    # point charge fluxes through the grid in xy-plane. The code seems to take a lot of time when the number of grid
    # points increases. Works fine until N in np.linspace is >41. The function calculate_point_charge_flux works
    # smoothly with bigger grid in the other 3 parts of the problem.

    x = np.linspace(-L/4, L/4, 11)
    y = np.linspace(-L/4, L/4, 11)

    side_plane = np.array([L/2, 0, 0])  # Flux through plane (L/2, y, z)
    total_flux = np.zeros([len(x), len(y)])
    for i in range(len(x)):
        for j in range(len(x)):
            q_pos = np.array([x[i], y[j], 0])  # Charge position is the xy-plane
            grid_point_flux = calculate_point_charge_flux(L, q_pos, side_plane, len(x))
            total_flux[i, j] = grid_point_flux

    # Plotting the contourf plot

    fig, ax = plt.subplots(1, 1)
    X, Y = np.meshgrid(x, y)
    cb = ax.contourf(X.T, Y.T, total_flux)
    fig.colorbar(cb)
    ax.set_title('Electric flux through the side/plane locating at x = L/2, with charge in xy-plane')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

main()