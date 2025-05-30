''' in this file we implement the symplectic Euler method for solving ODEs in pendulum dynamics. '''
import numpy as np
import matplotlib.pyplot as plt

def symplectic_pendulum(L,theta0,h,g = 9.81):
    period = 2 * np.pi * np.sqrt(L /g)
    t_max = 4 * period  
    steps = int(t_max / h) + 1


    theta = np.zeros(steps)

