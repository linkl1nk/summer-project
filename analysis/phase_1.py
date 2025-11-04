import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Explicit Euler
def explicit_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1, n):
        theta[i] = theta[i-1] + h * omega[i-1]
        omega[i] = omega[i-1] - h * a * np.sin(theta[i-1])
    return theta, omega

# Symplectic Euler
def symplectic_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1, n):
        theta[i] = theta[i-1] + h * omega[i-1]
        omega[i] = omega[i-1] - h * a * np.sin(theta[i])
    return theta, omega

# Implicit Euler
def implicit_euler(L, h, n, theta0, a, g):
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1, n):
        def func(x):
            return x + h * a * np.sin(theta[i-1] + h * x) - omega[i-1]
        omega[i] = fsolve(func, omega[i-1])[0]
        theta[i] = theta[i-1] + h * omega[i]
    return theta, omega

# 计算背景能量等高线
def energy_contour(L, g):
    theta_grid = np.linspace(-2*np.pi, 2*np.pi, 300)
    omega_grid = np.linspace(-3, 3, 300)
    Theta, Omega = np.meshgrid(theta_grid, omega_grid)
    E = 0.5 * (L * Omega)**2 + g * L * (1 - np.cos(Theta))
    return Theta, Omega, E

# 主程序
def pendulum(L, h, T, theta0, g=9.8):
    n = int(T / h)
    t = np.linspace(0, T, n)
    a = g / L

    # 三种方法求解
    theta_exp, omega_exp = explicit_euler(L, h, n, theta0, a, g)
    theta_imp, omega_imp = implicit_euler(L, h, n, theta0, a, g)
    theta_sym, omega_sym = symplectic_euler(L, h, n, theta0, a, g)

    # 额外求一个高精度RK4用于做"ideal path"
    def rk4(L, h, n, theta0, a, g):
        theta = np.zeros(n)
        omega = np.zeros(n)
        theta[0] = theta0
        omega[0] = 0
        for i in range(1, n):
            th, om, dt = theta[i-1], omega[i-1], h
            k1_th = om
            k1_om = -a * np.sin(th)
            k2_th = om + 0.5 * dt * k1_om
            k2_om = -a * np.sin(th + 0.5 * dt * k1_th)
            k3_th = om + 0.5 * dt * k2_om
            k3_om = -a * np.sin(th + 0.5 * dt * k2_th)
            k4_th = om + dt * k3_om
            k4_om = -a * np.sin(th + dt * k3_th)
            theta[i] = th + (dt / 6) * (k1_th + 2*k2_th + 2*k3_th + k4_th)
            omega[i] = om + (dt / 6) * (k1_om + 2*k2_om + 2*k3_om + k4_om)
        return theta, omega

    theta_rk4, omega_rk4 = rk4(L, h, n, theta0, a, g)

    # 计算背景能量等高线
    Theta_bg, Omega_bg, E_bg = energy_contour(L, g)

    # 画图
    plt.figure(figsize=(18, 5))

    # 图 1: Trajectory plot
    plt.subplot(1, 3, 1)
    circle_theta = np.linspace(-np.pi, np.pi, 100)
    plt.plot(L * np.sin(circle_theta), -L * np.cos(circle_theta), 'k--', alpha=0.3, label='Ideal path')
    plt.plot(L * np.sin(theta_rk4), -L * np.cos(theta_rk4), 'b-', label='RK4')
    plt.plot(L * np.sin(theta_sym), -L * np.cos(theta_sym), 'r-.', label='Symplectic Euler')
    plt.plot([0, 0], [-L*1.1, L*0.1], 'k-', alpha=0.5)
    plt.scatter(0, 0, c='k', marker='o', s=50)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title(f'Pendulum Trajectory (L={L}m, θ₀={theta0:.3f} rad, h={h}s)')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)

    # 图 2: Phase portrait (RK4)
    plt.subplot(1, 3, 2)
    plt.contour(Theta_bg, Omega_bg, E_bg, levels=30, cmap='gray', alpha=0.3)
    plt.plot(theta_rk4, omega_rk4, 'b-', label='RK4')
    plt.xlabel(r'$\theta$ (rad)')
    plt.ylabel(r'$\omega$ (rad/s)')
    plt.title('Phase Portrait - RK4')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')

    # 图 3: Phase portrait (Symplectic)
    plt.subplot(1, 3, 3)
    plt.contour(Theta_bg, Omega_bg, E_bg, levels=30, cmap='gray', alpha=0.3)
    plt.plot(theta_sym, omega_sym, 'r-.', label='Symplectic Euler')
    plt.xlabel(r'$\theta$ (rad)')
    plt.ylabel(r'$\omega$ (rad/s)')
    plt.title('Phase Portrait - Symplectic Euler')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')

    plt.tight_layout()
    plt.show()

# 运行示例
pendulum(L=1.0, h=0.001, T=20, theta0=np.pi/3, g=9.8)
