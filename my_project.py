import numpy as np
import matplotlib.pyplot as plt

def pendulum(L, h, g=9.8, theta0=0.1745):  # theta0 default ~10° in radians
    """
    Simulate a pendulum motion using three methods (all in radians):
    1. Exact numerical solution (RK4)
    2. Small angle approximation (linearized solution)
    3. Euler method (explicit Euler integration)
    
    """
    
    # Simulation duration (4 periods of small-angle approximation)
    small_angle_period = 2 * np.pi * np.sqrt(L / g)
    t_max = 4 * small_angle_period  
    num_steps = int(t_max / h) + 1  # Ensure integer number of steps
    t = np.linspace(0, t_max, num_steps)  # Time array with linspace

    # ========================================================================
    # 1. Small Angle Approximation (analytical solution)
    def small_angle_solution(t, theta0, L, g):
        omega_lin = np.sqrt(g / L)
    omega_lin = np.sqrt(g / L)  # Natural frequency (rad/s)
    theta_sa = theta0 * np.cos(omega_lin * t)  # Angular displacement (rad)
    omega_dot_sa = -theta0 * omega_lin * np.sin(omega_lin * t)  # Angular velocity (rad/s)

    # ========================================================================
    # 2. Euler Method
    theta_eu = np.zeros(num_steps)
    omega_eu = np.zeros(num_steps)
    theta_eu[0] = theta0
    omega_eu[0] = 0.0  # Released from rest

    for i in range(1, num_steps):
        alpha_eu = -(g / L) * np.sin(theta_eu[i-1])  # Angular acceleration (rad/s²)
        omega_eu[i] = omega_eu[i-1] + alpha_eu * h  # Explicit Euler step
        theta_eu[i] = theta_eu[i-1] + omega_eu[i-1] * h

    # ========================================================================
    # 3. RK4 Method (reference solution)
    # Use finer step size for accuracy (fixed ratio to h)
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
        dt = t_rk4[i] - t_rk4[i-1]  # Actual time step (constant for linspace)
        
        # RK4 coefficients
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

    # Interpolate RK4 to match Euler time points
    theta_exact = np.interp(t, t_rk4, theta_rk4)
    omega_exact = np.interp(t, t_rk4, omega_rk4)

    # ========================================================================
    # Visualization
    plt.figure(figsize=(14, 10))

    # Trajectory plot (convert to Cartesian coordinates)
    plt.subplot(2, 1, 1)
    # Ideal circle reference
    circle_theta = np.linspace(-np.pi, np.pi, 100)
    plt.plot(L * np.sin(circle_theta), -L * np.cos(circle_theta), 'k--', alpha=0.3, label='Ideal path')
    # Solutions
    plt.plot(L * np.sin(theta_exact), -L * np.cos(theta_exact), 'b-', label='Exact (RK4)')
    plt.plot(L * np.sin(theta_sa), -L * np.cos(theta_sa), 'g--', label='Small angle')
    plt.plot(L * np.sin(theta_eu), -L * np.cos(theta_eu), 'r-.', label='Euler method')
    plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5)  # Pivot line
    plt.scatter(0, 0, c='k', marker='o', s=50)  # Pivot point
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Pendulum Trajectory (L={L}m, θ₀={theta0:.3f} rad, h={h}s)')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)

    # Energy plot
    plt.subplot(2, 1, 2)
    energy_exact = 0.5 * (L * omega_exact)**2 + g * L * (1 - np.cos(theta_exact))
    energy_sa = 0.5 * (L * omega_dot_sa)**2 + g * L * (1 - np.cos(theta_sa))
    energy_eu = 0.5 * (L * omega_eu)**2 + g * L * (1 - np.cos(theta_eu))
    
    plt.plot(t, energy_exact, 'b-', label='Exact (RK4)')
    plt.plot(t, energy_sa, 'g--', label='Small angle')
    plt.plot(t, energy_eu, 'r-.', label='Euler method')
    plt.xlabel('Time (s)')
    plt.ylabel('Total energy (J/kg)')
    plt.title('Energy conservation comparison')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# Example usage (direct radians input):
pendulum(L=1.0, h=0.0001, theta0=np.pi/50) 