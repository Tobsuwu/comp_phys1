"""
Computational physics 1

1. Add code to function 'largest_eig'
- use the power method to obtain the largest eigenvalue and the 
  corresponding eigenvector of the provided matrix

2. Compare the results with scipy's eigs
- this is provided, but you should use that to validating your 
  power method implementation

Notice: 
  dot(A,x), A.dot(x), A @ x could be helpful for performing 
  matrix operations

"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.integrate import simps

def eigenvalue(A, v):
    # Calculates the eigenvalue using the dot product

    Av = A.dot(v)
    return v.dot(Av)

def largest_eig(A,tol=1e-9):
    """
    - Add simple power method code in this function
    - Add appropriate commenting especially in this comment section
    """

    v = np.ones(A.shape[1])/np.linalg.norm(A.shape[1])  # Creating a normalized unit vector for initial guess

    error = lambda_old = 100  # Initializing the error term for the first iteration of the loop

    while error > tol:  # Repeating until desired accuracy is received
        Av = A.dot(v)
        v_new = Av/np.linalg.norm(Av)  # The normalized iteration vector (lec3 slides equation (32))
        lambda_new = eigenvalue(A, v_new)  # Calculates the new eigenvalue
        error = abs(lambda_new - lambda_old)  # Calculates the error and exits the while loop if error < tol
        lambda_old = lambda_new  # Initializing values for next iteration
        v = v_new

    eig_value = lambda_new
    eig_vector = v_new

    return eig_value, eig_vector

def main():
    # Initializes the matrix which eigenvector and eigenvalues we wish to calculate and does plots on them.
    # The acquired plot from power method eigenvectors doesn't quite match the scipy version.


    grid = np.linspace(-5,5,100)
    grid_size = grid.shape[0]
    dx = grid[1]-grid[0]
    dx2 = dx*dx
    
    H0 = sp.diags(
        [
            -0.5 / dx2 * np.ones(grid_size - 1),
            1.0 / dx2 * np.ones(grid_size) - 1.0/(abs(grid)+1.5),
            -0.5 / dx2 * np.ones(grid_size - 1)
        ],
        [-1, 0, 1])
    
    eigs, evecs = sla.eigsh(H0, k=1, which='LA')
    
    l,vec=largest_eig(H0)
    
    print('largest_eig estimate: ', l)
    print('scipy eigsh estimate: ', eigs)
    
    psi0=evecs[:,0]
    norm_const=simps(abs(psi0)**2,x=grid)
    psi0=psi0/norm_const
    print(norm_const)

    psi0_ = vec*1.0
    norm_const=simps(abs(psi0_)**2,x=grid)
    psi0_=psi0_/norm_const
    print(norm_const)

    plt.plot(grid,abs(psi0)**2,label='scipy eig. vector squared')
    plt.plot(grid,abs(psi0_)**2,'r--',label='largest_eig vector squared')
    plt.text(-1, 2, "largest_eig estimate eigenvalue: " + str(l))
    plt.text(-1, 1.9, "scipy eigsh estimate eigenvalue : " + str(eigs))
    plt.legend(loc=0)
    plt.show()


if __name__=="__main__":
    main()
