''' in this file we implement the symplectic Euler method for solving ODEs in pendulum dynamics.

 the  key of this method is [theta_n+1, omega_n+1] = [theta_n, omega_n] + h*[oemga_n, -a*sin(theta_n+1)]

 '''
import numpy as np
import matplotlib.pyplot as plt

def symplectic_pendulum(L,theta0,h,g = 9.81):
    period = 2 * np.pi * np.sqrt(L /g)
    t_max = 4 * period  
    steps = int(t_max / h) + 1
    t = np.linspace(0, t_max, steps)
    print(t[2]-t[1])
    

    theta = np.zeros(steps)
    omega = np.zeros(steps)
    theta[0] = theta0
    omega[0] = 0

    for i in range(1, steps):
        alpha = np.sqrt(g / L)
        theta[i]= theta[i-1] + h * omega[i-1]
        omega[i]= omega[i-1] - alpha *h* np.sin(theta[i])

    # =============================== Visualization ===============================
    plt.figure(figsize=(10, 6))

    # === 1. Pendulum trajectory in Cartesian coordinates ===
    plt.subplot(2, 1, 1)

    x = L * np.sin(theta)
    y = -L * np.cos(theta)

    plt.plot(x, y, color='royalblue', linewidth=1.5, label='Symplectic Euler')
    plt.plot([0, 0], [-L * 1.1, L * 0.1], 'k--', linewidth=1)  # Pendulum axis
    plt.scatter(0, 0, c='black', s=50, label='Pivot')  # Pivot point

    plt.xlabel('x (m)', fontsize=12)
    plt.ylabel('y (m)', fontsize=12)
    plt.title(f'Pendulum Trajectory\nL = {L} m, θ₀ = {theta0:.2f} rad, h = {h}s', fontsize=14)
    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()

    # === 2. Total energy over time ===
    plt.subplot(2, 1, 2)

    # Energy per unit mass: E = ½(Lω)^2 + gL(1 - cos(θ))
    energy = 0.5 * (L * omega)**2 + g * L * (1 - np.cos(theta))

    plt.plot(t, energy, color='darkred', linewidth=1.5, label='Total Energy')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Energy (J/kg)', fontsize=12)
    plt.title('Energy Conservation Over Time', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()

    # === Display the full figure ===
    # plt.show()


# example
symplectic_pendulum(1,np.pi/3,0.0001,g = 9.81)



' this part will show how to use Newton method to solve cos(x) + x/2 = 0 '

def f(x):
    return np.cos(x) + x/2

def f_1(x):
    return -np.sin(x)+ 1/2

x = 0
while abs(f(x)) > 10 ** -5:
    x = x - f(x)/f_1(x)

print(x, f(x))