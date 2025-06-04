' this part will show how to use Newton method to solve cos(x) + x/2 = 0 '
import numpy as np

def f(x):
    return np.cos(x) + x/2

def f_1(x):
    return -np.sin(x)+ 1/2

x = 0
while abs(f(x)) > 10 ** -3:
    x = x - f(x)/f_1(x)

print(x, f(x))


# 2D newton

import numpy as np

def newton_system(F, J, x0, tol=1e-8, max_iter=50):
    """
    Solves a system of nonlinear equations F(x) = 0 using Newton's method in 2D.

    Parameters:
    - F: function that takes a vector x and returns a vector F(x)
    - J: function that takes a vector x and returns the Jacobian matrix J_F(x)
    - x0: initial guess vector (numpy array of shape (2,))
    - tol: convergence tolerance
    - max_iter: maximum number of iterations

    Returns:
    - x: the estimated root vector of the system F(x) = 0
    """
    x = x0.astype(float)  # Ensure using float for division

    for i in range(max_iter):
        Fx = F(x)                     # Evaluate F at current x
        Jx = J(x)                     # Evaluate Jacobian at current x

        # Solve J_F(x_n) * delta = -F(x_n)
        try:
            delta = np.linalg.solve(Jx, -Fx)
        except np.linalg.LinAlgError:
            raise ValueError("Jacobian is singular at iteration {}".format(i))

        x = x + delta                 # Update x

        # Check for convergence
        if np.linalg.norm(delta, ord=2) < tol:
            print(f"Converged in {i+1} iterations.")
            return x

    raise ValueError("Did not converge within the maximum number of iterations")

# Example usage
if __name__ == "__main__":
    # Define a sample system of equations F(x) = 0 where x = [x, y]
    def F(x):
        return np.array([
            x[0]**2 + x[1]**2 - 4,     # Circle: x^2 + y^2 = 4
            x[0] - x[1] - 1            # Line: x - y = 1
        ])

    # Jacobian matrix of F
    def J(x):
        return np.array([
            [2*x[0], 2*x[1]],
            [1,     -1    ]
        ])

    # Initial guess
    x0 = np.array([1.0, 1.0])

    # Solve
    solution = newton_system(F, J, x0)
    print("Solution:", solution)
