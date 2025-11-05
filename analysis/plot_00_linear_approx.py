import numpy as np
import matplotlib.pyplot as plt

def simulate_pendulum(L, h, theta0, g=9.8):
    small_angle_period = 2 * np.pi * np.sqrt(L / g)
    t_max = 4 * small_angle_period  
    num_steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, num_steps)

    # Small angle approximation
    omega_lin = np.sqrt(g / L)
    theta_sa = theta0 * np.cos(omega_lin * t)
    omega_sa = -theta0 * omega_lin * np.sin(omega_lin * t)

    # RK4 method
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

    theta_exact = np.interp(t, t_rk4, theta_rk4)
    omega_exact = np.interp(t, t_rk4, omega_rk4)

    energy_exact = 0.5 * (L * omega_exact)**2 + g * L * (1 - np.cos(theta_exact))
    energy_sa = 0.5 * (L * omega_sa)**2 + g * L * (1 - np.cos(theta_sa))

    return t, theta_exact, theta_sa, energy_exact, energy_sa

# =============================
# 统一计算两组数据
L = 1.0
h = 0.0001
g = 9.8

# 小角度：π/10
t_small, theta_exact_small, theta_sa_small, energy_exact_small, energy_sa_small = simulate_pendulum(L, h, np.pi/10, g)

# 大角度：π/6
t_large, theta_exact_large, theta_sa_large, energy_exact_large, energy_sa_large = simulate_pendulum(L, h, np.pi/6, g)

# 统一能量图纵坐标
max_energy = max(np.max(energy_exact_small), np.max(energy_sa_small),
                 np.max(energy_exact_large), np.max(energy_sa_large))
ylim = (0, max_energy * 1.1)

# =============================
# 小角度作图
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
circle_theta = np.linspace(-np.pi, np.pi, 200)
plt.plot(L * np.sin(circle_theta), -L * np.cos(circle_theta), 'k--', alpha=0.3, label='Ideal path')
plt.plot(L * np.sin(theta_exact_small), -L * np.cos(theta_exact_small), color='#66CCFF', label='RK4')
plt.plot(L * np.sin(theta_sa_small), -L * np.cos(theta_sa_small), 'g--', label='Small angle')
plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5)
plt.scatter(0, 0, c='k', marker='o', s=50)
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Small Angle (θ₀ = 18°)')
plt.axis('equal')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(t_small, energy_exact_small, color='#66CCFF', label='RK4')
plt.plot(t_small, energy_sa_small, 'g--', label='Small angle')
plt.ylim(ylim)
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.title('Energy (Small Angle)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# =============================
# 大角度作图
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(L * np.sin(circle_theta), -L * np.cos(circle_theta), 'k--', alpha=0.3, label='Ideal path')
plt.plot(L * np.sin(theta_exact_large), -L * np.cos(theta_exact_large), color='#66CCFF', label='RK4')
plt.plot(L * np.sin(theta_sa_large), -L * np.cos(theta_sa_large), 'g--', label='Small angle')
plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5)
plt.scatter(0, 0, c='k', marker='o', s=50)
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Large Angle (θ₀ = 30°)')
plt.axis('equal')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(t_large, energy_exact_large, color='#66CCFF', label='RK4')
plt.plot(t_large, energy_sa_large, 'g--', label='Small angle')
plt.ylim(ylim)
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.title('Energy (Large Angle)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
