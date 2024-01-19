import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse


def make_matrix(x):
    # Creation of the "Mystery matrix" on the given uniform grid. The mystery matrix calculates the derivative of the
    # function deduced by the plotted graphs.
    # Derivation: Nani?

    N = len(x)
    h = x[1]-x[0]
    off_diag = np.ones((N-1,))  # Creates an array with only values 1.
    A = np.zeros((N,N)) - np.diag(off_diag,-1) + np.diag(off_diag,1)
    A[0,0:3] = [-3.0,4.0,-1.0] 
    A[-1,-3:] = [1.0,-4.0,3.0]
    return A/h/2

def main():
    # Defines the functions which we imply the mystery matrix

    N = 50
    grid = np.linspace(0,np.pi,N)
    x1 = np.sin(grid)
    x2 = np.cos(grid)
    x3 = np.exp(grid)
    x4 = np.sin(grid)**2+np.cos(grid)**2


    A = make_matrix(grid)
    
    # calculate here b=Ax as described in the problem setup
    b1 = np.dot(A, x1)
    b2 = np.dot(A, x2)
    b3 = np.dot(A, x3)
    b4 = np.dot(A, x4)

    # Plotting the 4 plots in the same figure
    fig = plt.figure()
    f1 = fig.add_subplot(221)
    f2 = fig.add_subplot(222)
    f3 = fig.add_subplot(223)
    f4 = fig.add_subplot(224)

    f1.plot(grid, b1, grid, x1, '--')
    f1.legend(["b=Ax", "x = sin(x)"])
    f2.plot(grid, b2, grid, x2, '--')
    f2.legend(["b=Ax", "x = cos(x)"])
    f3.plot(grid, b3, grid, x3, '--')
    f3.legend(["b=Ax", "x = e^x"])
    f4.plot(grid, b4, grid, x4, '--')
    f4.legend(["b=Ax", "x = sin(x)^2 + cos(x)^2"])
    fig.suptitle("Different types of matrix vector products")
    plt.show()

if __name__=="__main__":
    main()



