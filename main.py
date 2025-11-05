import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sys
import os

# --- Add parent directory to sys.path ---
# This allows the script to find and import the 'src' module
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.join(script_dir, '..')
sys.path.append(parent_dir)

# --- Import required solvers from the 'src' module ---
# Assumes 'src/solvers.py' now contains the 6-ARGUMENT symplectic_euler
from src.solvers import explicit_euler, implicit_euler, symplectic_euler

#constant
G = 9.8

def simulate_pendulum(L=1.0, h=0.01, T=10.0, theta0=np.pi/4):
    a = G / L  # coefficient in the differential equation
    N = int(T / h)  # number of time steps
    t = np.linspace(0, T, N)

    # Initial angular velocity is zero
    omega0 = 0.0

    # Run all three methods
    theta_exp, omega_exp = explicit_euler(theta0, omega0, h, N, a)
    theta_imp, omega_imp = implicit_euler(theta0, omega0, h, N, a)
    
    # --- THIS IS THE FIX ---
    # Call the 6-argument version from src/solvers.py
    # (L, h, n, theta0, a, g)
    theta_sym, omega_sym = symplectic_euler(L, h, N, theta0, a, G)
    # --- END OF FIX ---

    # Compute energy: E = (1/2) * (L * omega)^2 + g * L * (1 - cos(theta))
    energy_exp = 0.5 * (L * omega_exp)**2 + G * L * (1 - np.cos(theta_exp))
    energy_imp = 0.5 * (L * omega_imp)**2 + G * L * (1 - np.cos(theta_imp))
    energy_sym = 0.5 * (L * omega_sym)**2 + G * L * (1 - np.cos(theta_sym))

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
if __name__ == "__main__":
    # We use h=0.001 (as decided) for a more accurate 'main' plot
    simulate_pendulum(L=1.0, h=0.001, T=10.0, theta0=np.pi/4)