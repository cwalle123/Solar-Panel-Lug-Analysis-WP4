import numpy as np

def bearing_check(candidate_vector, force_vector, material):
    # Extract material property
    sigma_bearing = material['sigma_bearing']  # [Pa]

    # Extract candidate variables
    t_2 = candidate_vector[1]  # Plate thickness
    x = candidate_vector[2]  # Fastener pitch in x
    D_2 = candidate_vector[4]  # Fastener diameter
    w = candidate_vector[6]  # Width of the plate
    e_1 = candidate_vector[7]  # Edge distance in y
    e_2 = candidate_vector[8]  # Edge distance in x

    # Extract applied forces
    Fx = force_vector[0]  # Force in x-direction
    Fz = force_vector[2]  # Force in z-direction
    My = force_vector[4]  # Moment about y-axis

    # Constants
    nf = 4  # Number of fasteners

    # Compute geometric parameters
    dist = np.sqrt((x / 2 - e_2) ** 2 + (w / 2 - e_1) ** 2)  # Distance to fasteners
    A_hole = 0.25 * np.pi * D_2 ** 2  # Area of a fastener hole

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
        F_in_plane_M_y * np.cos(angles),  # x-components
        F_in_plane_M_y * np.sin(angles)  # y-components
    ]).T

    # Combine forces
    P = force_components + np.array([F_in_plane_x, F_in_plane_z])

    # Compute magnitudes of resultant forces
    P_magnitudes = np.sqrt(P[:, 0] ** 2 + P[:, 1] ** 2)

    # Maximum resultant force
    P_max = np.max(P_magnitudes)

    # Compute bearing stress and margin of safety
    sigma_bearable = P_max / (t_2 * D_2)
    MS = (sigma_bearing / sigma_bearable) - 1

    return MS









