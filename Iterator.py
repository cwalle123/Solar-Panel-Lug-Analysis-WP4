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





