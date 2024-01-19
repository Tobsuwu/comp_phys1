""" This program loads "experimental data" and interpolates it using spline interpolations and visualizes it."""

# Import the needed packages

import numpy as np
import matplotlib.pyplot as plt
import spline_class


def experimental_data():
    # Loads the data from the given files and returns the x- and y-values as [1, 10] and the experimental data Fexp
    # as a 10*10 grid.

    xexp = np.loadtxt('x_grid.txt')
    yexp = np.loadtxt('y_grid.txt')
    Fexp = np.loadtxt('exp_data.txt')
    Fexp = Fexp.reshape([len(xexp), len(yexp)])
    return xexp, yexp, Fexp

def visualize(xexp, yexp, Fexp, f, x):
    # Function that visualizes the data using subplots.

    X, Y = np.meshgrid(xexp, yexp)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Problem 2')
    ax1.contourf(X.T, Y.T, Fexp)
    ax1.scatter(X.T, Y.T, color='red')
    ax1.plot(x, 2*x**2, linewidth=5, color='black')
    ax1.set_title('Experimental data, grid points and desired path')
    ax2.plot(x, f)
    ax2.scatter(x, np.loadtxt('ref_interpolated.txt'), color='red')
    ax2.set_title('Result of the interpolation')
    ax2.legend(["Interpolated data", "reference points"])
    plt.show()

def main():
    # Calculates the interpolation using the spline class module. Experimental data is given as parameters to the class
    # spline and interpolated along the path y=2*x^2

    xexp, yexp, Fexp = experimental_data()
    spl2d = spline_class.spline(x=xexp, y=yexp, f=Fexp, dims=2)
    x = np.linspace(0, 1, 100) # y has to be in range [0, 2]
    f = np.zeros_like(x)
    for i in range(len(x)):
        f[i] = spl2d.eval2d(x[i], 2*x[i]**2)
    visualize(xexp, yexp, Fexp, f, x)

main()