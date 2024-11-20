"""
checks.py contains verification of all different failure modes
each failure mode is given its own unique function that will be used
by calculator.py to verify an overall design
"""

# External imports
import numpy as np
import math


# elastic modulus [Pa], yield strength [Pa], shear strength [Pa], density [kg/m^3], thermal expansion coefficient [1/K]
material_properties = {
    'Al 7075 T6': (71.7e9, 430e6, .5*430e6, 2810, 2.36e-5)
}


def vectorized(func):
    """decorator used to easily vectorise any function"""
    return np.vectorize(func)


def convention_for_safety_margin_calculation(inputs):
    """this one just checks one failure criterion and outputs the name of the criterion as well as the safety margin"""
    ...
    # first is the name of the failure mode check, then the safety margin
    return str, float


def convention_for_verification_check(inputs):
    """this one is a true/false output on whether some condition is valid. for eg. we should have a check for whether
    all the lengths are actually greater than zero or stuff like that"""
    ...
    # first is name of the verification check, then the True/False on whether the check is passed
    return str, bool


def positivity_check(design_vector):
    """Verifies no length is given as negative, as that would be nonsensical"""
    if np.all(design_vector > 0):
        return True
    else:
        return False


def X_compliance_check(t_1, X, e_2, h):
    """Verifies the horizontal direction passes thread clearance"""
    if h + 2*t_1 + 4*e_2 <= X:
        return True
    else:
        return False


def w_compliance_check(w, D_1, D_2, e_1):
    """Verifies the vertical direction passes thread clearance"""
    if 2*e_1 + 2*D_2 <= w and 3*D_1 <= w:
        return True
    else:
        return False


def thread_edge_distance_check(e_1, e_2, D_2):
    if e_1 >= 1.5*D_2 and e_2 >= 1.5*D_2:
        return True
    else:
        return False


def bearing_check(candidate_vector, force_vector, material):
    # Extract material property
    sigma_bearing = material_properties[material][1]

    # Extract candidate variables
    t_2 = candidate_vector[1]
    x = candidate_vector[2]
    D_2 = candidate_vector[4]
    w = candidate_vector[6]
    e_1 = candidate_vector[7]
    e_2 = candidate_vector[8]

    # Extract applied forces
    Fx = force_vector[0]
    Fz = force_vector[2]
    My = force_vector[4]

    # Constants
    nf = 4  # Number of fasteners

    # Compute geometric parameters
    dist = np.sqrt((x / 2 - e_2) ** 2 + (w / 2 - e_1) ** 2)
    A_hole = 0.25 * np.pi * D_2 ** 2

    # Compute forces on fasteners
    F_in_plane_x = Fx / nf
    F_in_plane_z = Fz / nf
    F_in_plane_M_y = (My * A_hole * dist) / (nf * A_hole * dist ** 2)

    # Fastener locations
    fasteners = np.array([
        [-(x / 2 - e_2), -(w / 2 - e_1)],
        [-(x / 2 - e_1), (w / 2 - e_1)],
        [(x / 2 - e_2), -(w / 2 - e_2)],
        [(x / 2 - e_2), (w / 2 - e_2)]
    ])

    # Compute angles for fasteners
    angles = np.arctan2(fasteners[:, 1], fasteners[:, 0])

    # Compute force components
    force_components = np.vstack([
        F_in_plane_M_y * np.cos(angles),
        F_in_plane_M_y * np.sin(angles)
    ]).T

    # Combine forces
    P = force_components + np.array([F_in_plane_x, F_in_plane_z])

    # Compute magnitudes of resultant forces
    P_magnitudes = np.sqrt(P[0] ** 2 + P[1] ** 2)

    # Maximum resultant force
    P_max = np.max(P_magnitudes)

    # Compute bearing stress and margin of safety
    sigma_bearable = P_max / (t_2 * D_2)
    MS = (sigma_bearing / sigma_bearable) - 1

    return MS

# Objective function
@vectorized
def mass(Density_Lug, Thickness_Lug_1, Thickness_Lug_2, w_Lug, TotalLength_Lug, Diameter_Lug_1, Diameter_Lug_2): #all inputs should be given in SI units
    """This function calculates the total mass of the lug by calculating and summing the volumes of
    the hinges and the backplate, and then multiplies the total volume by the density of the material"""
    HingeArea_Lug = w_Lug**(2) + (w_Lug**(2) * math.pi) / (8) - (Diameter_Lug_1**(2) * math.pi) / (4)  # m^2
    TotalHingeVolume_Lug = 2 * Thickness_Lug_1 * HingeArea_Lug  # m^3
    BackplateArea_Lug = TotalLength_Lug * w_Lug - Diameter_Lug_2**(2) * math.pi  # m^2
    BackplateVolume_Lug = Thickness_Lug_2 * BackplateArea_Lug  # m^3
    TotalVolume = TotalHingeVolume_Lug + BackplateVolume_Lug  # m^3
    return Density_Lug * TotalVolume  # kg

