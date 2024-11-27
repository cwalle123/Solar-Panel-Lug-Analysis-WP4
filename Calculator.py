
"""
calculator.py verifies the validity of a single design, i.e. a single set of design parameters.
it runs instances of every check function and compiles the resulting safety margins to either
print them out or pass them on to the iterator.py
it also computes the mass of the design
"""

# External imports
import numpy as np
import pandas as pd

# Internal imports
from checks import *


# F_x, F_y, F_z, M_x, M_y, M_z
applied_load_vector = np.array([134, 134, 409, 2.39, 2.39, 2.39])


def verify_validity(candidate_vector, fastener):
    """Runs a multitude of different checks for verifying the validity of the design candidate
    Includes physical compatibility checks like non-negative lengths, no clipping of thread holes etc. """
    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector

    if not positivity_check(candidate_vector):
        return False

    if not thread_size_check(D_1, D_2, fastener):
        return False

    if not thread_edge_distance_check(e_1, e_2, D_2):
        return False

    if not w_compliance_check(w, D_1, D_2, e_1):
        return False

    if not X_compliance_check(t_1, X, e_2, h):
        return False

    return True


def compute_safety_margins(candidate_vector, material, fastener, load_vector):
    """Checks the physically relevant failure criteria and produces relevant safety factors. If any are negative, the
    computation immediately halts"""

    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector
    # Collect material data
    elastic_modulus, sigma_ult, tau_ult, rho, material_thermal_expansion = material_properties[material]
    # Collect temperature data
    T_min, T_max, T_ref = temperatures
    # Collect fastener data
    d_head, d_shank, l_shank, d_sm, pitch, A_fast, A_stiff, fast_strength, fast_thermal_expansion = fastener_properties[fastener]

    bearing_check_SM, SC_wall_SM = bearing_check(candidate_vector, load_vector, material)

    # TODO FIGURE OUT THE FORCE RATIO
    thermal_check_SM = thermal_load_fastener(fast_thermal_expansion, material_thermal_expansion, elastic_modulus, A_stiff, .5, T_ref, T_max, T_min, A_fast, fast_strength)

    lug_pullthrough_SM, wall_pullthrough_SM = pull_through_check(load_vector, candidate_vector, material, fastener)

    safety_margins = np.array([bearing_check_SM, SC_wall_SM, thermal_check_SM, lug_pullthrough_SM, wall_pullthrough_SM])

    return safety_margins

def calculator(candidate_vector, material, fastener, load_vector):

    # Unpack the candidate vector
    t_1, t_2, X, D_1, D_2, h, w, e_1, e_2 = candidate_vector
    # Collect material data
    elastic_modulus, sigma_ult, tau_ult, rho, thermal_expansion = material_properties[material]
    # Unpack the load vector
    F_x, F_y, F_z, M_x, M_y, M_z = load_vector

    # Compatibility verification
    if not verify_validity(candidate_vector, fastener):
        # Configuration failed, check case be case to print out which criteria failed:
        # if not positivity_check(candidate_vector):
        #     print("Positivity check failed")
        #
        # if not thread_edge_distance_check(e_1, e_2, D_2):
        #     print("Thread-edge distance check failed")
        #
        # if not w_compliance_check(w, D_1, D_2, e_1):
        #     print("w compliance check failed")
        #
        # if not X_compliance_check(t_1, X, e_2, h):
        #     print("X compliance check failed")

        raise ValueError("Configuration fails physical compatibility")

    # Safety margin computation
    safety_margins = compute_safety_margins(candidate_vector, material, fastener, load_vector)

    # Compute the objective function
    mass_of_design = mass(rho, t_1, t_2, w, X, D_1, D_2)

    # Show results
    # print(safety_margins)
    # print(mass_of_design)

    return safety_margins, mass_of_design

def main():
    # Read input CSV using pandas
    input_df = pd.read_csv('input.csv', encoding='utf-8-sig')

    # Extract relevant columns for processing
    candidate_vectors = input_df.iloc[:, :9].astype(float).to_numpy()
    applied_load_vectors = input_df.iloc[:, 9:15].astype(float).to_numpy()
    materials = input_df.iloc[:, -1]

    # Process each row and store the results
    results = []
    for candidate_vector, applied_load_vector, material in zip(candidate_vectors, applied_load_vectors, materials):
        safety_margins, mass_of_design = calculator(candidate_vector, material, applied_load_vector)
        results.append([safety_margins, mass_of_design])

    # Create a new DataFrame for the results
    output_df = pd.DataFrame(results, columns=['safety_margins', 'mass_of_design'])

    # Write the output to a CSV
    output_df.to_csv('output.csv', index=False)


if __name__ == "__main__":
    main()
