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


def verify_validity(candidate_vector):
    """Runs a multitude of different checks for verifying the validity of the design candidate
    Includes physical compatibility checks like non-negative lengths, no clipping of thread holes etc. """
    if not positivity_check(candidate_vector):
        return False
    ...
    return True

def compute_safety_margins(candidate_vector, material):
    """Checks the physically relevant failure criteria and produces relevant safety factors. If any are negative, the
    computation immediately halts"""

    # Collect material data
    elastic_modulus, sigma_ult, tau_ult, rho, thermal_expansion = material_properties[material]

    ...


def calculator(candidate_vector):

    # temporary variable assignment, all in SI, use the reader for reference
    t_1, t_2, t_3, D_1, D_2, h, w, e_1, e_2, D_fo, D_fi = candidate_vector

    # Compatibility verification
    if not verify_validity(candidate_vector):
        # Configuration failed, check case be case to print out which criteria failed:
        ...
        # raise error to stop the program
        raise ValueError("Configuration fails physical compatibility")

    # Safety margin computation
    safety_margins = compute_safety_margins(candidate_vector)


    # Compute the objective function
    mass_of_design = mass(candidate_vector)

    return safety_margins, mass_of_design

