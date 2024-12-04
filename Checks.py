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
    'Al 7075 T6': (71.7e9, 430e6, .5 * 430e6, 2810, 2.36e-5),
    'Ti-6Al-4V': (113.8e9, 880e6, .5 * 880e6, 4430, 8.6e-6),
    '18-8 Stainless Steel': (193e9, 205e6, .5 * 205e6, 7900, 16.5e-6)
}

# head diameter [m], shank diameter [m], shank length [m], d_sm [m], pitch [m], fastener area [m^2], stiffness area [m^2], tensile strength [Pa], thermal expansion coefficient [1/K]
fastener_properties = {
    'M6x1': (0.01, 0.006, 0.01, 0.0048, 0.001, 2.83e-5, 1.79e-5, 8.96e8, material_properties['Ti-6Al-4V'][-1]),
    'M0.5x0.125': (0.00094, 0.0005, 0.0015, 0.000347, 0.000125, 1.96e-7, 9.44e-8, 4.83e8,
                   material_properties['18-8 Stainless Steel'][-1]),
    'M1x0.25': (
    0.002, 0.001, 0.0015, 0.000693, 0.00025, 7.85e-7, 3.77e-7, 4.83e8, material_properties['18-8 Stainless Steel'][-1])
}

# minimum operating, maximum operating, reference (manufacturing) [K]
temperatures = (250, 300, 288)

# TODO SC wall thickness
t_3 = .01  # [m]


# Objective function
def mass(Density_Lug, Thickness_Lug_1, Thickness_Lug_2, w_Lug, TotalLength_Lug, Diameter_Lug_1,
         Diameter_Lug_2):  #all inputs should be given in SI units
    """This function calculates the total mass of the lug by calculating and summing the volumes of
    the hinges and the backplate, and then multiplies the total volume by the density of the material"""
    HingeArea_Lug = w_Lug ** (2) + (w_Lug ** (2) * math.pi) / (8) - (Diameter_Lug_1 ** (2) * math.pi) / (4)  # m^2
    TotalHingeVolume_Lug = 2 * Thickness_Lug_1 * HingeArea_Lug  # m^3
    BackplateArea_Lug = TotalLength_Lug * w_Lug - Diameter_Lug_2 ** (2) * math.pi  # m^2
    BackplateVolume_Lug = Thickness_Lug_2 * BackplateArea_Lug  # m^3
    TotalVolume = TotalHingeVolume_Lug + BackplateVolume_Lug  # m^3
    return Density_Lug * TotalVolume  # kg


def positivity_check(design_vector):
    """Verifies no length is given as negative, as that would be nonsensical"""
    if np.all(design_vector > 0):
        return True
    else:
        return False


def beam_fit_check(h):
    beam_diameter = 0.03
    if h >= beam_diameter:
        return True
    else:
        return False


def thread_size_check(D_1, D_2, fastener):
    """Verifies the lug is big enough to fit the threads"""
    d_shank = fastener_properties[fastener][1]
    if D_1 >= d_shank and D_2 >= d_shank:
        return True
    else:
        return False


def X_compliance_check(t_1, X, e_2, h):
    """Verifies the horizontal direction passes thread clearance"""
    if h + 2 * t_1 + 4 * e_2 <= X:
        return True
    else:
        return False


def w_compliance_check(w, D_1, D_2, e_1):
    """Verifies the vertical direction passes thread clearance"""
    if 2 * e_1 + 2 * D_2 <= w and 2 * D_1 <= w:
        return True
    else:
        return False


def thread_edge_distance_check(e_1, e_2, D_2):
    """Verifies the thread holes have sufficient clearance from lug edges"""
    if e_1 >= 1.5 * D_2 and e_2 >= 1.5 * D_2:
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
    F_in_plane_M_y = (My) / (nf * dist)

    # Fastener locations
    fasteners = np.array([
        [-(x / 2 - e_2), -(w / 2 - e_1)],
        [-(x / 2 - e_2), (w / 2 - e_1)],
        [(x / 2 - e_2), -(w / 2 - e_1)],
        [(x / 2 - e_2), (w / 2 - e_1)]
    ])

    # Compute angles for fasteners
    angles = np.arctan2(fasteners[:, 1], fasteners[:, 0])
    if My >= 0:
        angles_forces = angles + np.pi / 2
    else:
        angles_forces = angles - np.pi / 2

    # Compute force components
    force_components = np.vstack([
        F_in_plane_M_y * np.cos(angles_forces),
        F_in_plane_M_y * np.sin(angles_forces)]).T

    # Combine forces
    P = force_components + np.array([F_in_plane_x, F_in_plane_z])

    # Compute magnitudes of resultant forces
    P_magnitudes = np.sqrt(P[0] ** 2 + P[1] ** 2)

    # Maximum resultant force
    P_max = np.max(P_magnitudes)

    # Compute bearing stress and margin of safety
    sigma_bearable = P_max / (t_2 * D_2)
    MS = (sigma_bearing / sigma_bearable) - 1

    # Compute the bearing of the wall
    sigma_yield_wall = 345e6
    sigma_wall = P_max / (t_3 * D_2)
    MS_wall = (sigma_yield_wall / sigma_wall) - 1

    return MS, MS_wall


def thermal_load_fastener(thermalExp_coef_fastener, thermalExp_coef_clamped, E, A_stiffness, force_ratio,
                          temp_ref, temp_working_max, temp_working_min, area_fastener, stress_max_fastener):
    """calculate the force induced by thermal expansion
        -Material of both fastener and clamped is needed
        -"""

    min_temp = temp_working_min - temp_ref

    max_temp = temp_working_max - temp_ref

    min_thermal_load = (thermalExp_coef_fastener - thermalExp_coef_clamped) * min_temp * E * A_stiffness * (
            1 - force_ratio)

    max_thermal_load = (thermalExp_coef_fastener - thermalExp_coef_clamped) * max_temp * E * A_stiffness * (
            1 - force_ratio)

    safety_factor = stress_max_fastener / (max([min_thermal_load, max_thermal_load]) / area_fastener) - 1

    # we only care about maximum so we take max

    return safety_factor


def pull_through_check(applied_load_vector, candidate_vector, material, fastener):
    # Material properties
    yield_stress = material_properties[material][1]  # Shear yield strength [Pa]
    D_fo = fastener_properties[fastener][0]  # Outer fastener head diameter [m]
    D_fi = fastener_properties[fastener][1]  # Inner fastener hole diameter [m]

    # Extract applied loads
    F_x = applied_load_vector[0]
    F_y = applied_load_vector[1]  # Force in y-direction [N]
    F_z = applied_load_vector[2]
    M_x = applied_load_vector[3]
    M_y = applied_load_vector[4]
    M_z = applied_load_vector[5]  # Moment about z-axis [Nm]
    nf = 4  # Number of fasteners

    # Extract candidate parameters
    t_2 = candidate_vector[1]  # Plate thickness [m]
    x = candidate_vector[2]  # Fastener pitch in x [m]
    D_2 = candidate_vector[4]  # Fastener diameter [m]
    w = candidate_vector[6]  # Plate width [m]
    e_1 = candidate_vector[7]  # Edge distance in y [m]
    e_2 = candidate_vector[8]  # Edge distance in x [m]

    fasteners = np.array([
        [-(x / 2 - e_2), -(w / 2 - e_1)],
        [-(x / 2 - e_2), (w / 2 - e_1)],
        [(x / 2 - e_2), -(w / 2 - e_1)],
        [(x / 2 - e_2), (w / 2 - e_1)]
    ])

    # Compute geometric properties
    dist = np.sqrt((x / 2 - e_2) ** 2 + (w / 2 - e_1) ** 2)  # Distance to fasteners [m]
    A_hole = 0.25 * np.pi * D_2 ** 2  # Area of a fastener hole [m²]
    A_shear = np.pi * D_2 * (t_2 + t_3)  # Combined shear area for plate and wall [m²]
    A_normal = np.pi * (D_fo ** 2 - D_fi ** 2) / 4  # Annular shear area under the fastener head [m²]

    # Compute forces on fasteners
    M_z_total = - F_x * w + M_z
    M_x_total = F_z * w + M_x
    F_in_plane_y = F_y / nf - M_x_total / ((w - 2 * e_1) * 2) + M_z_total / (
                (x - 2 * e_2) * 2)  # In-plane force due to F_y [N]

    # Compute stresses
    sigma_normal = F_in_plane_y / A_normal  # Normal stress due to F_y [Pa]
    tau_t_2 = F_in_plane_y / (np.pi * D_2 * t_2)  # Shear stress on the plate [Pa]
    tau_t_3 = F_in_plane_y / (np.pi * D_2 * t_3)  # Shear stress on the wall [Pa]

    # Von Mises stresses
    sigma_VM_t_2 = np.sqrt(sigma_normal ** 2 + 3 * tau_t_2 ** 2)  # Von Mises stress on the plate [Pa]
    sigma_VM_t_3 = np.sqrt(sigma_normal ** 2 + 3 * tau_t_3 ** 2)  # Von Mises stress on the wall [Pa]

    # Margins of safety
    MS_t_2 = (yield_stress / sigma_VM_t_2) - 1  # Margin of safety for the plate
    MS_t_3 = (yield_stress / sigma_VM_t_3) - 1  # Margin of safety for the wall

    return MS_t_2, MS_t_3


def flange_check(applied_load_vector, candidate_vector, material, fastener):
    # Material properties
    yield_stress = material_properties[material][1]  # Shear yield strength [Pa]

    # Extract applied loads
    F_x = applied_load_vector[0]  # Force in x-direction [N]
    F_y = applied_load_vector[1]  # Force in y-direction [N]
    F_z = applied_load_vector[2]  # Force in z-direction [N]

    # Extract candidate parameters
    D_1 = candidate_vector[3]  # Hole diameter [m]
    t_1 = candidate_vector[0]  # Flange thickness [m]
    w = candidate_vector[6]  # Flange height [m]
    e = w / 2

    # Find the relevant areas
    def compute_Aav(w, D, e, t):
        A1 = (w/2 - D/2 * np.cos(np.pi/4)) * t
        A4 = A1
        A2 = (w/2 - D/2) * t
        A3 = (e - D/2) * t
        Aav = 6 / (3/A1 + 1/A2 + 1/A3 + 1/A4)
        return Aav

    A_t = (w - D_1) * t_1  # Axial area
    A_br = D_1 * t_1  # Bearing area
    A_av = compute_Aav(w, D_1, e, t_1)

    #  Pre-defined fitting functions obtained by sampling points from the pictures in the reader
    def compute_Ktu(area_ratio):
        if area_ratio < 0:
            print("K_tu computation outside of valid range")
            return np.NaN
        elif area_ratio > 1.4:
            area_ratio = 1.4
        return 1.642 * area_ratio - 1.307 * area_ratio**2 + 0.342 * area_ratio**3

    def compute_Kbru(eD):
        if eD < 0.5:
            print("K_bru computation outside of valid range")
            return np.NaN
        elif eD > 2.25:
            return 1.7
        return -1.443 + 3.647 * eD - 1.648 * eD**2 + 0.277 * eD**3

    def compute_Kt(wD):
        if wD < 1 or wD > 5:
            print("K_t computation outside of valid range")
            return np.NaN
        return 1.244 - 0.351 * wD + 0.127 * wD**2 - 0.015 * wD**3

    K_tu = compute_Ktu(A_av/A_br)
    K_bru = compute_Kbru(e/D_1)
    K_t = compute_Kt(w/D_1)

    # Find the relevant margins
    P_tu = K_tu * A_br * yield_stress
    P_bru = K_bru * yield_stress * A_br
    P_u = K_t * yield_stress * A_t

    R_a = F_y / min(P_u, P_bru)
    R_tr = F_z / P_tu

    # Compute the safety margin
    deviation = np.power(R_a, 1.6) + np.power(R_tr, 1.6)
    SM = 1/np.power(deviation, 0.625) - 1

    return SM



