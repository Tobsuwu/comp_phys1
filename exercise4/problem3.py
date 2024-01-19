""" This program calculates the electron density along the line. """

import numpy as np
import read_xsf_example as rxe


def main():

    r_0 = np.array([[0.1], [0.1], [2.8528]])  # Unit: Ångström
    r_1 = np.array([[4.45], [4.45], [2.8528]])  # The parametrized equation of 3D line is given as: x = (-r_0[x]_0+r_1[x])*t etc.
    t = np.linspace(0, 1, 400)
    r = np.zeros((len(t), 3))
    for i in range(len(t)):
        r[i][0] = 4.35*t[i] + 0.1
        r[i][1] = 4.35*t[i] + 0.1
        r[i][2] = 2.8528
    rho, lattice, grid, shift = rxe.read_example_xsf_density('dft_chargedensity1.xsf')
    # Using the function from problem 2 we get the numerical density rho. We would then use our line to get
    # points from the rho matrix. z-axis stays constant while we take the values corresponding to the line.

main()