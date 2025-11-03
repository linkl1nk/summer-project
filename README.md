# Numerical Simulation of a Pendulum

This project was completed as part of the UCL Summer Project (2025). It explores and compares various numerical methods for simulating the motion of a simple pendulum, which is a classic non-linear dynamical system.

---

## 1. Project Goal

The primary goal is to solve the pendulum's ordinary differential equation (ODE) and compare the results of different numerical integrators. The key metrics for comparison are:
* **Accuracy:** How close the simulation stays to a known analytical solution (e.g., the small-angle approximation or a high-precision RK4 solution).
* **Stability & Energy Conservation:** How the numerical solution behaves over long time periods, particularly whether it artificially gains or loses energy.

---

## 2. Core Methods Implemented

This project implements and analyzes the following integrators:

1.  **Small-Angle Approximation (Analytical Solution)**
    * **File:** `small_angel.py` (You may want to rename this to `small_angle.py`)
    * **Description:** The exact analytical solution for the *linearized* ODE ($\sin(\theta) \approx \theta$). This serves as a baseline for comparison at small initial angles.

2.  **Explicit Euler Method**
    * **File:** `...` (Which file implements this? `pendulum.py` perhaps?)
    * **Description:** A basic, first-order numerical method. It is simple to implement but known for its instability and failure to conserve energy.

3.  **Symplectic Euler Method (or Semi-Implicit Euler)**
    * **File:** `symplectic_euler.py`
    * **Description:** A first-order method that is popular in physics. It is designed to (almost) conserve the system's energy over long simulations, making it far more stable than the explicit Euler method.

4.  **(Other Methods - Help me identify these)**
    * **File:** `newton.py`
        * **Description:** `...` (What does this file do? Does it use Newton's method to solve something? Or is it related to the `Newton.py` file in the original project?)
    * **File:** `phase_2.py`
        * **Description:** `...` (Does this plot the "phase space" or "phase diagram" of the pendulum? i.e., plotting angular velocity `omega` vs. angle `theta`?)

---

## 3. How to Set Up and Run

This project requires Python and standard scientific libraries.

**1. Clone the repository:**
```bash
git clone [https://github.com/linkl1nk/summer-project.git](https://github.com/linkl1nk/summer-project.git)
cd summer-project