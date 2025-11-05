# Pendulum ODE: A Reflection on Numerical Solvers

This project began as a UCL Summer Project to compare numerical methods for solving the pendulum ODE. It has been extended into a reflection on the "correct" choice of algorithm, balancing short-term accuracy, long-term stability, and computational complexity.

This project's narrative argues that for physical simulations, "structure-preserving" algorithms (like Symplectic Euler) are often superior to "high-accuracy" algorithms (like RK4), especially over long time horizons.

---

## Project Narrative & Key Analyses

This repository is organized around a central narrative. The scripts in the `analysis/` folder are numbered to follow this "story."

### Act 0: The Limit of Linearization
* **Script:** `analysis/plot_00_linear_approx.py`
* **Goal:** Justify why numerical methods are needed.
* **Analysis:** Compares the analytical "small angle" solution (sin(θ) ≈ θ) against the high-accuracy RK4 solver.
* **Finding:** The approximation holds for small angles (e.g., 18°) but fails significantly at larger angles (e.g., 30°), proving we must solve the full non-linear ODE.

### Act 1: The 3 Eulers - A Search for Stability
* **Script:** `analysis/plot_01_euler_comparison.py`
* **Goal:** Compare the three base Euler methods.
* **Analysis:** Uses Energy Plots and Phase Portraits to compare Explicit, Implicit, and Symplectic Euler.
* **Finding:** Explicit Euler is unstable (gains energy). Implicit Euler is "damped" (loses energy). Symplectic Euler is the only one that conserves energy (oscillates around the true value), making it the "winner" for physical simulation.
* **Evidence:** `plot_02_implicit_stability.py` (proves Implicit loses energy long-term) and `plot_03_symplectic_stability.py` (proves Symplectic conserves energy long-term).

### Act 2: The Reflection - RK4 vs. Symplectic
* **Script:** `analysis/plot_04_rk4_reflection.py`
* **Goal:** A deeper reflection comparing the "best" Euler (Symplectic) against the "gold standard" (RK4).
* **Analysis:** Runs a long-term simulation (`T=200`) with a large step-size (`h=0.05`).
* **Finding:** This reveals the core reflection of the project:
    1.  The **Trajectory Plot** becomes a useless, overlapping mess.
    2.  The **Phase Portrait** is the *only* tool that shows the truth: `RK4` (non-symplectic) slowly drifts and spirals away (energy drift). `Symplectic Euler` (structure-preserving) remains perfectly stable and bounded.
    3.  **Conclusion:** `RK4` wins on short-term accuracy, but `Symplectic Euler` wins on long-term stability.

### Explorations
* **Script:** `scratchpad/newton_method_exploration.py`
* **Purpose:** A "learning note" exploring how Newton's Method works from scratch. This was done to understand the `scipy.fsolve` solver used by the Implicit Euler method.

---

## Project Structure

* `main.py`: A simple "dashboard" script that runs all 3 Euler methods and plots the 4-panel comparison.
* `src/solvers.py`: **The Single Source of Truth.** Contains all core solver functions (Euler, RK4, small_angle).
* `analysis/`: Contains all the specialized, "publication-ready" plots for the narrative.
* `scratchpad/`: Contains early drafts and learning explorations (like `newton.py`).
* `requirements.txt`: A list of all Python dependencies.
* `.gitignore`: Configured to ignore `.venv/` and `__pycache__/`.

---

## How to Run

1.  Clone the repository and `cd` into it.
2.  Create and activate a virtual environment:
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```
3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
4.  Run the main dashboard:
    ```bash
    python main.py
    ```
5.  Run a specific analysis from the narrative:
    ```bash
    python analysis/plot_01_euler_comparison.py
    ```

---

## Future Work

This project's structure is now clean, but the analysis is not yet complete.

* [ ] **Refactor Imports:** The scripts in `analysis/` still need to be fully refactored to `import` their solvers from `src/solvers.py`.
* [ ] **Parameter Tuning:** The `analysis/` plots (especially `plot_04_rk4_reflection.py`) need their `h` and `T` parameters carefully "tuned" to ensure their trends are "not blurry" and persuasively demonstrate the scientific argument.
* [ ] **Big O Analysis:** Add a formal Big O complexity analysis to this README, comparing the computational cost of each method (e.g., `Implicit` vs. `Symplectic`).