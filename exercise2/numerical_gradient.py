""" This program calculates the N dimensional numerical gradient of the given function. """
import numpy as np


def f(x, N):
    F = np.zeros((N, 1))
    if N < 3:
        for i in range(N):
            F[i] = np.sin(x[i]) + np.cos([N])
        return F
   #for i in range(N-1):
        #f[i] = np.sin(x[i]) + np.cos([N-1]) + np.sum

    #return np.sin(x[0]) + np.cos([N-1]) + np.sum(x[k]**2, where=(k>=2 and k<=N-1))

def main():
    N = 1
    x = np.zeros((N, 1))
    function = f(x, N)
    gradient = np.gradient(np.array([function[0]], dtype=object))
    print(gradient)

if __name__=="__main__":

    main()