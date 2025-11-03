import numpy as np
import matplotlib.pyplot as plt
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

def symplectic_euler(theta0, omega0, h, N, a):
    theta = np.zeros(N)
    omega = np.zeros(N)
    theta[0], omega[0] = theta0, omega0

    for i in range(1, N):
        omega[i] = omega[i-1] - h * a * np.sin(theta[i-1])
        theta[i] = theta[i-1] + h * omega[i]

    return theta, omega