import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def implicit_euler(theta0, omega0, h, N, a):
    theta = np.zeros(N)
    omega = np.zeros(N)
    theta[0], omega[0] = theta0, omega0

    for i in range(1, N):
        # Solve for omega_{n+1} using fsolve
        def func(w_new):
            return w_new - omega[i-1] + h * a * np.sin(theta[i-1] + h * w_new)
        
        omega_new = fsolve(func, omega[i-1])[0]
        theta_new = theta[i-1] + h * omega_new

        omega[i] = omega_new
        theta[i] = theta_new

    return theta, omega

def pendulum(L, h, g=9.8, theta0=0.1745):  
    """
    Simulate a pendulum motion using RK4 and Implicit Euler.
    """
    # Simulation duration (long enough to show divergence behavior)
    small_angle_period = 2 * np.pi * np.sqrt(L / g)
    t_max = 100 * small_angle_period  
    num_steps = int(t_max / h) + 1
    t = np.linspace(0, t_max, num_steps)

    # RK4 solution (reference)
    h_rk4 = h / 10
    num_steps_rk4 = int(t_max / h_rk4) + 1
    t_rk4 = np.linspace(0, t_max, num_steps_rk4)
    theta_rk4 = np.zeros(num_steps_rk4)
    omega_rk4 = np.zeros(num_steps_rk4)
    theta_rk4[0] = theta0
    omega_rk4[0] = 0.0

    for i in range(1, num_steps_rk4):
        th = theta_rk4[i-1]
        om = omega_rk4[i-1]
        dt = t_rk4[i] - t_rk4[i-1]
        
        k1_th = om
        k1_om = -(g / L) * np.sin(th)
        
        k2_th = om + 0.5 * dt * k1_om
        k2_om = -(g / L) * np.sin(th + 0.5 * dt * k1_th)
        
        k3_th = om + 0.5 * dt * k2_om
        k3_om = -(g / L) * np.sin(th + 0.5 * dt * k2_th)
        
        k4_th = om + dt * k3_om
        k4_om = -(g / L) * np.sin(th + dt * k3_th)
        
        theta_rk4[i] = th + (dt / 6) * (k1_th + 2*k2_th + 2*k3_th + k4_th)
        omega_rk4[i] = om + (dt / 6) * (k1_om + 2*k2_om + 2*k3_om + k4_om)

    # Interpolate RK4 onto Euler grid
    theta_exact = np.interp(t, t_rk4, theta_rk4)
    omega_exact = np.interp(t, t_rk4, omega_rk4)

    # Implicit Euler solution
    theta_ie, omega_ie = implicit_euler(theta0, 0.0, h, num_steps, g/L)

    # Visualization
    plt.figure(figsize=(10, 8))

    # Trajectory plot
    plt.subplot(2, 1, 1)
    circle_theta = np.linspace(-np.pi, np.pi, 100)
    plt.plot(L * np.sin(circle_theta), -L * np.cos(circle_theta), 'k--', alpha=0.3, label='Ideal path')
    plt.plot(L * np.sin(theta_exact), -L * np.cos(theta_exact), 'b-', label='RK4')
    plt.plot(L * np.sin(theta_ie), -L * np.cos(theta_ie), 'r-.', label='Implicit Euler')
    plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5)
    plt.scatter(0, 0, c='k', marker='o', s=50)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Pendulum Trajectory (L={L}m, θ₀={theta0:.3f} rad, h={h}s)')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)

    # Energy plot
    plt.subplot(2, 1, 2)
    energy_exact = 0.5 * (L * omega_exact)**2 + g * L * (1 - np.cos(theta_exact))
    energy_ie = 0.5 * (L * omega_ie)**2 + g * L * (1 - np.cos(theta_ie))
    
    plt.plot(t, energy_exact, 'b-', label='RK4')
    plt.plot(t, energy_ie, 'r-.', label='Implicit Euler')
    plt.xlabel('Time (s)')
    plt.ylabel('Total energy (J/kg)')
    plt.title('Energy conservation comparison')
    plt.legend()
    plt.grid(True)

    plt.ylim(0, 1.5 * energy_exact[0])

    plt.tight_layout()
    plt.show()

# Example usage:
pendulum(L=1.0, h=0.0001, theta0=np.pi/4)
