import numpy as np
import matplotlib.pyplot as plt

def plot_pendulum(L, h, g=9.8, theta0=np.pi/6):
    # in the function, I will use three different methods to slolve the pendulum problem:
    # 1. Exact nonlinear solution using solve_ivp



    # 2. Small-angle linear approximation
    a = np.sqrt(g / L)
    theta_linear = theta0 * np.cos(a * np.arange(0, 10, h))
    omega_linear = -theta0 * a * np.sin(a * np.arange(0, 10, h))


    # 3. Euler method
    

    
