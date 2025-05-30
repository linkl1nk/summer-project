import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# basic value setting
# ----------------------------

g = 9.81       # gravitational acceleration
L = 1.0        # length of string
a = g * L      # consatant
dt = 0.01      # step size
T = 10         # total time
N = int(T / dt)   # number of steps

# initial condition
theta0 = np.pi/120  # initial angle: pi/120
omega0 = 0          # initial angular speed

# ----------------------------
# initiallize array
# y: exact value ， ya: approximation

# we have t0,t1,……,tN so we have N points on graph. For each point they have [ω, θ] 
# so we have the shape of N rows and 2 columns
# ----------------------------
y = np.zeros((N, 2))
ya = np.zeros((N, 2))
t = np.linspace(0, T, N)

y[0] = [omega0, theta0]
ya[0] = [omega0, theta0]


# ----------------------------
# transfer 2nd ODE to sequence of 1st ODEs
# exact: d(omega)/dt = -a * sin(theta)
#        d(theta)/dt = omega
# approximation: d(omega)/dt = -a * theta
#                d(theta)/dt = omega

# use explicit Euler 
# exact : omega n+1 = omega n - dt*(g/L)*sin(theta n)
#         theta n+1 = theta n + dt*omega n
# approximation: omega n+1 = omega n - dt*(g/L)*(theta n)
#                theta n+1 = theta n + dt*omega n
# ----------------------------

for i in range(N - 1):
# exact value:
    omega, theta = y[i] # the omega and theta for ith point
    # omega n+1 = omega n - dt*(g/L)*sin(theta n)
    y[i+1, 0] = omega + dt * (-g/L) * np.sin(theta)
    # theta n+1 = theta n + dt*omega n
    y[i+1, 1] = theta + dt * omega

# approximation
    omega_a, theta_a = ya[i]
    # omega n+1 = omega n - dt*(g/L)*(theta n)
    ya[i+1, 0] = omega_a + dt * (-g/L) * theta_a
    # theta n+1 = theta n + dt*omega n
    ya[i+1, 1] = theta_a + dt * omega_a

# ----------------------------
# check energy conservation
# ----------------------------
Ek = 0.5 * y[:, 0]**2              # first column of y ---- angular speed
Ep = a * (1 - np.cos(y[:, 1]))     # second column of y ---- theta
Et = Ek + Ep

Eka = 0.5 * ya[:, 0]**2
Epa = a * (1 - np.cos(ya[:, 1]))
Eta = Eka + Epa

# ----------------------------
# make graph
# ----------------------------

plt.figure(figsize=(12, 5))

# trajectory
plt.subplot(1, 2, 1)
plt.plot(np.sin(y[:,1]), -np.cos(y[:,1]), label='Non-linear')
plt.plot(np.sin(ya[:,1]), -np.cos(ya[:,1]), label='Approx (linear)', linestyle='--')
plt.title('Pendulum Trajectory')
plt.xlabel('x = sin(θ)')
plt.ylabel('y = -cos(θ)')
plt.legend()
plt.axis('equal')

# energy
plt.subplot(1, 2, 2)
plt.plot(t, Et, label='Non-linear Energy')
plt.plot(t, Eta, label='Approx Energy', linestyle='--')
plt.title('Total Energy Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.legend()

plt.tight_layout()
plt.show()
