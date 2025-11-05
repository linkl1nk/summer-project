# analysis/plot_01_phase_necessity.py

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# --- Add parent directory to sys.path ---
# This allows the script to find and import the 'src' module
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.join(script_dir, '..')
sys.path.append(parent_dir)

# --- Import required solvers from the 'src' module ---
# We assume src/solvers.py has:
#   symplectic_euler(theta0, omega0, h, N, a) -> 5 args
#   rk4(L, h, n, theta0, a, g) -> 6 args
from src.solvers import rk4, symplectic_euler

# --- Local Helper Functions ---

def energy_contour(L, g):
    """
    Generates a grid of energy values for the phase portrait background.
    """
    theta_grid = np.linspace(-2*np.pi, 2*np.pi, 300)
    omega_grid = np.linspace(-4, 4, 300) # Increased omega range
    Theta, Omega = np.meshgrid(theta_grid, omega_grid)
    E = 0.5 * (L * Omega)**2 + g * L * (1 - np.cos(Theta))
    return Theta, Omega, E

# --- Main Analysis Function ---

def run_analysis_phase_necessity(L, h, T, theta0, g=9.8):
    """
    Runs the 'Phase Necessity' analysis.
    
    GOAL: Prove that trajectories can look identical
    while their phase portraits reveal underlying differences.
    """
    
    n = int(T / h)
    t = np.linspace(0, T, n)
    a = g / L

    # --- 1. Run Solvers ---
    # Call using the argument signatures from src/solvers.py
    theta_sym, omega_sym = symplectic_euler(L, h, n, theta0, a, g)
    theta_rk4, omega_rk4 = rk4(L, h, n, theta0, a, g)

    # --- 2. Calculate Background Contours ---
    Theta_bg, Omega_bg, E_bg = energy_contour(L, g)

    # --- 3. Plotting ---
    print(f"Running 'Phase Necessity' analysis (h={h}, T={T})...")
    plt.figure(figsize=(18, 5))
    plt.suptitle(f'Argument 1: Phase Portrait Necessity (h={h}, T={T}s)', fontsize=16)

    # --- Plot 1: Trajectory (Demonstrates Goal 1: Overlap) ---
    plt.subplot(1, 3, 1)
    
    # GOAL 1: Plot RK4 in a light color, then over-plot Symplectic
    plt.plot(L * np.sin(theta_rk4), -L * np.cos(theta_rk4), 'b-', alpha=0.3, label='RK4 (Underneath)')
    plt.plot(L * np.sin(theta_sym), -L * np.cos(theta_sym), 'r--', label='Symplectic (On Top)')
    
    plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5) # Pivot
    plt.scatter(0, 0, c='k', marker='o', s=50) # Pivot
    plt.title('Trajectory (They look identical!)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)

    # --- Plot 2: Phase Portrait (RK4) ---
    plt.subplot(1, 3, 2)
    plt.contour(Theta_bg, Omega_bg, E_bg, levels=30, cmap='gray', alpha=0.3)
    plt.plot(theta_rk4, omega_rk4, 'b-', label='RK4')
    plt.title('Phase Portrait - RK4')
    plt.xlabel(r'$\theta$ (rad)')
    plt.ylabel(r'$\omega$ (rad/s)')
    plt.legend()
    plt.grid(True)
    # plt.axis('equal') # DO NOT USE: This hides the subtle differences

    # --- Plot 3: Phase Portrait (Symplectic) ---
    plt.subplot(1, 3, 3)
    plt.contour(Theta_bg, Omega_bg, E_bg, levels=30, cmap='gray', alpha=0.3)
    plt.plot(theta_sym, omega_sym, 'r--', label='Symplectic Euler')
    plt.title('Phase Portrait - Symplectic (Different Shape!)')
    plt.xlabel(r'$\theta$ (rad)')
    plt.ylabel(r'$\omega$ (rad/s)')
    plt.legend()
    plt.grid(True)
    # plt.axis('equal') # DO NOT USE: This hides the subtle differences

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# --- Run the Script ---
if __name__ == "__main__":
    # "Large h, Small T" parameters to prove Goal 2a
    pendulum_L = 1.0
    pendulum_h = 0.01      # "Coarse" step size
    pendulum_T = 20        # "Short" time
    pendulum_theta0 = np.pi/3 # 60 degrees
    
    run_analysis_phase_necessity(L=pendulum_L, h=pendulum_h, T=pendulum_T, theta0=pendulum_theta0)