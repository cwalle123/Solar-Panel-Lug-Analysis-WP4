"""
calculator.py verifies the validity of a single design, i.e. a single set of design parameters.
it runs instances of every check function and compiles the resulting safety margins to either
print them out or pass them on to the iterator.py
it also computes the mass of the design
"""

# External imports
import numpy as np

# Internal imports
from checks import *


# F_x, F_y, F_z, M_x, M_y, M_z
applied_load_vector = np.array([0, 0, 426, 0, 0, 0])

def verify_validity(candidate_vector):
    """Runs a multitude of different checks for verifying the validity of the design candidate
    Includes physical compatibility checks like non-negative lengths, no clipping of thread holes etc. """
    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector

    if not positivity_check(candidate_vector):
        return False

    if not thread_edge_distance_check(e_1, e_2, D_2):
        return False

    if not w_compliance_check(w, D_2):
        return False

    if not X_compliance_check(t_1, X, D_2, h):
        return False

    return True


def compute_safety_margins(candidate_vector, material):
    """Checks the physically relevant failure criteria and produces relevant safety factors. If any are negative, the
    computation immediately halts"""

    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector
    # Collect material data
    elastic_modulus, sigma_ult, tau_ult, rho, thermal_expansion = material_properties[material]

    safety_margins = np.array([1, 2, 3, 4, 5])

    return safety_margins


def calculator(candidate_vector, material):

    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector
    # Collect material data
    elastic_modulus, sigma_ult, tau_ult, rho, thermal_expansion = material_properties[material]

    # Compatibility verification
    if not verify_validity(candidate_vector):
        # Configuration failed, check case be case to print out which criteria failed:
        if not positivity_check(candidate_vector):
            print("Positivity check failed")

        if not thread_edge_distance_check(e_1, e_2, D_2):
            print("Thread-edge distance check failed")

        if not w_compliance_check(w, D_2):
            print("w compliance check failed")

        if not X_compliance_check(t_1, X, D_2, h):
            print("X compliance check failed")

        raise ValueError("Configuration fails physical compatibility")

    # Safety margin computation
    safety_margins = compute_safety_margins(candidate_vector, material)

    # Compute the objective function
    mass_of_design = mass(rho, t_1, t_2, w, X, D_1, D_2)

    # Show results
    # print(safety_margins)
    # print(mass_of_design)

    return safety_margins, mass_of_design


def main():
    """Runs one iteration of the calculator with user defined inputs"""
    print("Please input your design parameters, all in meters")
    t_1 = input("t_1: ")
    t_2 = input("t_2: ")
    X = input("X: ")
    D_1 = input("D_1: ")
    D_2 = input("D_2: ")
    h = input("h: ")
    w = input("w: ")
    e_1 = input("e_1: ")
    e_2 = input("e_2: ")

    candidate_vector = np.array([t_1, t_2, X, D_1, D_2, h, w, e_1, e_2]).astype(float)
    print("Using material Al 7075 T6")

    dummy = input("Press Enter to run analysis")

    safety_margins, mass_of_design = calculator(candidate_vector, "Al 7075 T6")

    print("Here is the safety margins: ")
    print(safety_margins)

    print(f"And the mass of the design is {mass_of_design} kg")


if __name__ == "__main__":
    main()
