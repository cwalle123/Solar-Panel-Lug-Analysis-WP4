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
import logging


# Internal imports
from calculator import *
from checks import *


def random_neighbour(candidate, safety_factors, step_size=0.001):
    """Given a certain solution, randomly generate a new one"""
    # Main random step
    new_candidate = candidate + np.random.randint(-1, 2, size=candidate.shape) * step_size

    # Accelerated stepping for very high safety margins
    margins_parameters_relation = np.array([[0, 1, 1, 0, 0, 0, 1, 0, 0],
                                            [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                            [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                            [0, 1, 0, 0, 1, 0, 1, 0, 0],
                                            [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                            [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    if max(safety_factors) > 1:
        max_index = np.argmax(safety_factors)
        adjustment = margins_parameters_relation[max_index] * np.random.randint(0, 2, size=candidate.shape) * step_size
        return new_candidate - adjustment

    return new_candidate


def simulated_annealing_logged(start_point, max_iterations):
    """The primary iteration loop with logging. Logs the outcome of every 100th iteration"""
    logging.basicConfig(
        filename="simulated_annealing.log",  # Name of the log file
        filemode="w",  # Overwrite the file on each run
        level=logging.INFO,  # Logging level
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Some constants for now
    material = "Al 7075 T6"
    fastener = "M6x1"
    target_safety_margins = np.array([0, 0, 0, 0, -10, 0.15])  # TODO correct once appropriate safety factors are known

    # Finding the initial condition
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = start_point
    candidate = start_point
    safety_margins, best_mass = calculator(candidate, material, fastener, applied_load_vector)

    logging.info(f"Initial candidate: {candidate}, Mass: {best_mass}, Margins:{safety_margins}")

    # Primary annealing loop
    for step in range(max_iterations):
        # Generate a new solution candidate
        new_candidate = random_neighbour(candidate, safety_margins)

        # Evaluate validity, then the safety margins and mass (physically incompatible design give a ValueError)
        try:
            safety_margins, candidate_mass = calculator(new_candidate, material, fastener, applied_load_vector)
        except ValueError:
            if step % 100 == 0:
                logging.info(f"Step {step}: Configuration doesn't pass physical compatibility")
            continue

        if step % 100 == 0:
            logging.info(f"Step {step}: {candidate}, Mass: {candidate_mass}, Margins: {safety_margins}")

        if np.any(safety_margins < target_safety_margins):
            # Margins fail to meet criteria
            continue

        # Valid solution found, perform annealing
        T = max(.0001, ((max_iterations - step) / max_iterations) ** 5 - .005)

        if candidate_mass < best_mass:
            candidate = new_candidate
            best_mass = candidate_mass
        else:
            prob = np.exp(- abs(best_mass - candidate_mass) / T)
            if random() < prob:
                candidate = new_candidate
                best_mass = candidate_mass

    logging.info(f"Final candidate: {candidate}, Mass: {best_mass}, Margins: {safety_margins}")
    logging.info(f"Final candidate [mm]: {candidate*1000}, Mass [g]: {best_mass*1000}, Margins: {safety_margins}")
    return candidate, best_mass, safety_margins


def single_run():
    first_design = np.array([.01, .01, .15, .01, .01, .03, .10, .02, .02])
    first_design = np.array([0.001, 0.001, 0.198, 0.016, 0.006, 0.045, 0.033, 0.01, 0.02])

    max_iterations = 100000
    final_design, final_mass, final_margins = simulated_annealing_logged(first_design, max_iterations)
    print("Optimization outcome: ")
    print(final_design)
    print(final_mass)
    print(final_margins)


if __name__ == "__main__":
    single_run()

