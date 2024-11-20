"""
iterator.py is the primary logic loop of the software, which iterates upon the design variables
in order to converge on the most optimal design
given a certain initial set of design variables, it calls the calculator and evaluates the
out-coming safety factors and mass to adjust the design variables and run the calculator again
employing a certain optimisation strategy, iterator.py tries to converge
on the optimal design configuration
"""

# External imports
import numpy as np
from random import random

# Internal imports
from calculator import *
from checks import *


"""Pseudocode idea made with help from GPT-4

max_iterations = _
material_used = _
start_point = _
current_solution = start_point
for iteration in range(max_iterations):
    # Generate a candidate solution
    candidate = generate_neighbor(current_solution)

    # Verify the solution
    if not verify_solution(candidate):
        continue

    # Compute objective function
    objective = compute_objective(candidate)

    # Simulated annealing step
    if accept_solution(objective, current_objective, temperature):
        current_solution = candidate
        current_objective = objective

    # Update temperature
    temperature = cool_down(temperature)
"""


def random_neighbour(candidate, step_size=0.001):
    """Given a certain solution, randomly generate a new one"""
    return candidate + np.random.randint(-1, 2, size=candidate.shape) * step_size


def simulated_annealing(start_point, max_iterations):
    """The primary iteration loop, which TBD"""

    # Some constants for now
    material = "Al 7075 T6"
    target_safety_margins = np.array([0, 1, 1, 1, 1])  # TODO fill this out

    # F_x, F_y, F_z, M_x, M_y, M_z
    applied_load_vector = np.array([0, 0, 426, 0, 0, 0])

    # Finding the initial condition
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = start_point
    candidate = start_point
    best_mass = mass(material_properties[material][3], t_1, t_2, w, X, D_1, D_2)

    print(candidate)
    print(best_mass)

    for step in range(max_iterations):
        new_candidate = random_neighbour(candidate)

        try:
            safety_margins, candidate_mass = calculator(new_candidate, material, applied_load_vector)
        except ValueError:
            continue

        if np.any(safety_margins < target_safety_margins):
            continue

        # Valid solution found, perform annealing
        T = max(.0001, ((max_iterations - step) / max_iterations) ** 3 - .005)

        if candidate_mass < best_mass:
            candidate = new_candidate
            best_mass = candidate_mass
            continue
        else:
            prob = np.exp(- abs(best_mass - candidate_mass) / T)
            if random() < prob:
                candidate = new_candidate
                best_mass = candidate_mass

    return candidate, best_mass, safety_margins


def test():
    first_design = np.array([.01, .01, .15, .01, .01, .01, .10, .02, .02])
    max_iterations = 100000
    final_design, final_mass, final_margins = simulated_annealing(first_design, max_iterations)
    print("Optimization outcome: ")
    print(final_design)
    print(final_mass)
    print(final_margins)


if __name__ == "__main__":
    test()

"""
iterator.py is the primary logic loop of the software, which iterates upon the design variables
in order to converge on the most optimal design
given a certain initial set of design variables, it calls the calculator and evaluates the
out-coming safety factors and mass to adjust the design variables and run the calculator again
employing a certain optimisation strategy, iterator.py tries to converge
on the optimal design configuration
"""

# External imports
import numpy as np

# Internal imports
from calculator import *
from checks import *

"""Pseudocode idea made with help from GPT-4"""
# Example of main loop
max_iterations = _
material_used = _
start_point = _
current_solution = start_point
for iteration in range(max_iterations):
    # Generate a candidate solution
    candidate = generate_neighbor(current_solution)

    # Verify the solution
    if not verify_solution(candidate):
        continue

    # Compute objective function
    objective = compute_objective(candidate)

    # Simulated annealing step
    if accept_solution(objective, current_objective, temperature):
        current_solution = candidate
        current_objective = objective

    # Update temperature
    temperature = cool_down(temperature)





