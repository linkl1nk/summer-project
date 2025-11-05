import numpy as np
from scipy.optimize import fsolve


def explicit_euler(theta0, omega0, h, N, a):
    theta = np.zeros(N)
    omega = np.zeros(N)
    theta[0], omega[0] = theta0, omega0
    
    for i in range(1, N):
        omega[i] = omega[i-1] - h * a * np.sin(theta[i-1])
        theta[i] = theta[i-1] + h * omega[i-1]
    
    return theta, omega

def implicit_euler(theta0, omega0, h, N, a):
    theta = np.zeros(N)
    omega = np.zeros(N)
    theta[0], omega[0] = theta0, omega0

    for i in range(1, N):
        # Define a function whose root gives the new omega
        def func(w_new):
            return w_new - omega[i-1] + h * a * np.sin(theta[i-1] + h * w_new)
        
        omega_new = fsolve(func, omega[i-1])[0]
        theta_new = theta[i-1] + h * omega_new

        omega[i] = omega_new
        theta[i] = theta_new

    return theta, omega

def symplectic_euler(L, h, n, theta0, a, g):    
    theta = np.zeros(n)
    omega = np.zeros(n)
    theta[0] = theta0
    omega[0] = 0
    for i in range(1,n):
        theta[i] = theta[i-1] + h * omega[i-1]
        omega[i] = omega[i-1] - h * a * np.sin(theta[i]) # Uses NEW theta
    return theta, omega

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

def small_angle_solution(t, theta0, L, g):
    """
    Calculates the analytical (exact) solution for the
    LINEARIZED small-angle approximation.
    """
    omega_lin = np.sqrt(g / L)
    theta_sa = theta0 * np.cos(omega_lin * t)
    omega_sa = -theta0 * omega_lin * np.sin(omega_lin * t)
    return theta_sa, omega_sa