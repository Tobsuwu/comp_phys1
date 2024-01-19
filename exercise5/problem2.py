""" Trajectory of a charged particle influenced by both electric and magnetic field. """
""" Lorentz force: F = q*(E + v x B) """
""" Newton 2: F_tot = ma -> q*(E + v x B) = ma -> q*(E + dx/dt x B) = m*dx^2/dt^2 """
""" -> dv/dt = (q*(E + v(t) x B))/m """

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

def func(x, t, E, B):
    # Calculates the Lorentz force F/m. E = qE/m and B = qB/m. x is the vector containing the velocities and
    # corresponding positions. e.g. x = [vx, vy, vz, rx, ry, rz]
    # xx contains the accelerations calculated from the motion equation and the velocities copied from the x vector
    # e.g. xx = [ax, ay, az, vx, vy, vz]

    xx = np.zeros_like(x)
    xx[0:3] = E + np.cross(x[0:3], B)
    xx[3:] = x[0:3]
    return xx

def main():
    t_0 = 0  # Initial time
    x_0 = np.array([0.1, 0.1, 0.1, 0, 0, 0])  # Initial position. [vx, vy, vz, rx, ry, rz]
    E = np.array([0.05, 0, 0])  # qE/m = 0.05
    B = np.array([0, 4.0, 0])  # qB/m = 4
    t = np.linspace(t_0, 5, 100)  # Creating the time during we wish to calculate the trajectory

    # Using Scipys odeint to solve the ODE. It calculates the derivative of the x in respect to time, returning us the
    # velocity vectors and position vectors with respect to time t.
    # e.g ode_sol[i] = [vx[i], vy[i], vz[i], rx[i], ry[i], rz[i]

    ode_sol = odeint(func, x_0, t, args=(E, B))

    # Acquiring the space coordinates in respect to time from the ODE solution
    x = ode_sol[:, 3]
    y = ode_sol[:, 4]
    z = ode_sol[:, 5]

    # Plotting the trajectory of the particle using the space coordinates in respect to time
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label='r(t)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Trajectory of the particle under Lorentz force')
    ax.legend()
    plt.show()

    print('Velocity vector at time t = 5: ', ode_sol[-1, 0:3])

main()