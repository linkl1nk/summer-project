'In the following part, three diferent Euler method will be used to decribe the motion of a pendulun.'
'Four plottings will be genrated to show the defference between the three method'
'three methods are: 1) explicit Euler 2) implicit Euler 3) symplectic Euler'
'plottings show: 1) trajectory 2) total energy of the system 3)phase portrait 4) theta -t'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# governing equation is theta_2nd_derevative + a*sin(theta) = 0'

def explicit_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1,n):
        theta[i] = theta[i-1] + h * omega[i-1]
        omega[i] = omega [i-1] - h * a * np.sin(theta[i-1])
    return theta, omega

def symplectic_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1,n):
        theta[i] = theta[i-1] + h * omega[i-1]
        omega[i] = omega [i-1] - h * a * np.sin(theta[i])
    return theta, omega

def implicit_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    
    for i in range(1,n):
        def func(x):
            return x + h * a * np.sin (theta[i-1] + h * x) - omega[i-1]
        omega[i] = fsolve(func,omega[i-1])[0]
        theta[i] = theta[i-1] + h * omega[i]
    return theta, omega

def pendulum(L, h, T, theta0, g= 9.8):
    n = int(T/h)
    t = np.linspace(0, T, n)
    a = g/L

    # Run all three methods
    theta_exp, omega_exp = explicit_euler(L, h, n, theta0, a, g)
    theta_imp, omega_imp = implicit_euler(L, h, n, theta0, a, g)
    theta_sym, omega_sym = symplectic_euler(L, h, n, theta0, a, g)

    # Compute energy: E = (1/2) * (L * omega)^2 + g * L * (1 - cos(theta))
    energy_exp = 0.5 * (L * omega_exp)**2 + g * L * (1 - np.cos(theta_exp))
    energy_imp = 0.5 * (L * omega_imp)**2 + g * L * (1 - np.cos(theta_imp))
    energy_sym = 0.5 * (L * omega_sym)**2 + g * L * (1 - np.cos(theta_sym))

    # Plotting results
    plt.figure(figsize=(12, 10))

    # Plot 1: Real trajectory of the pendulum (theta -> x, y)
    plt.subplot(2, 2, 1)
    plt.plot(L * np.sin(theta_exp), -L * np.cos(theta_exp), label='Explicit Euler')
    plt.plot(L * np.sin(theta_imp), -L * np.cos(theta_imp), label='Implicit Euler')
    plt.plot(L * np.sin(theta_sym), -L * np.cos(theta_sym), label='Symplectic Euler')
    plt.title('Pendulum Trajectory')
    plt.xlabel('x-position (m)')
    plt.ylabel('y-position (m)')
    plt.legend()
    plt.axis('equal')

    # Plot 2: Total energy over time
    plt.subplot(2, 2, 2)
    plt.plot(t, energy_exp, label='Explicit Euler')
    plt.plot(t, energy_imp, label='Implicit Euler')
    plt.plot(t, energy_sym, label='Symplectic Euler')
    plt.title('Total Energy of the System')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (J)')
    plt.legend()

    # Plot 3: Phase portrait (omega vs theta)
    plt.subplot(2, 2, 3)
    plt.plot(theta_exp, omega_exp, label='Explicit Euler')
    plt.plot(theta_imp, omega_imp, label='Implicit Euler')
    plt.plot(theta_sym, omega_sym, label='Symplectic Euler')
    plt.title('Phase Portrait')
    plt.xlabel('Theta (rad)')
    plt.ylabel('Omega (rad/s)')
    plt.legend()

    # Plot 4: Theta over time
    plt.subplot(2, 2, 4)
    plt.plot(t, theta_exp, label='Explicit Euler')
    plt.plot(t, theta_imp, label='Implicit Euler')
    plt.plot(t, theta_sym, label='Symplectic Euler')
    plt.title('Theta over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Theta (rad)')
    plt.legend()

    plt.tight_layout()
    plt.show()

# Example of how to run the simulation
pendulum(L=1.0, h=0.01, T=10.0, theta0=np.pi/4,g=9.8)
    

    
