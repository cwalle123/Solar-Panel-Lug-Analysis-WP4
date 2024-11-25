import numpy as np


def pull_through_check(applied_load_vector, candidate_vector, material, fastener):

    # Material properties
    yield_stress = material_properties[material][1]  # Shear yield strength [Pa]
    # TODO Update once fasteners are implemented
    D_fo = fastener_properties[fastener][0]  # Outer fastener head diameter [m]
    D_fi = fastener_properties[fastener][0]  # Inner fastener hole diameter [m]

    # Extract applied loads
    F_x = candidate_vector[0]
    F_y = applied_load_vector[1]  # Force in y-direction [N]
    F_z =applied_load_vector[2]
    M_x = candidate_vector[3]
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

    # TODO DELETE LATER WHEN FASTENER PROPERTIES WORK
    D_fi = D_2
    D_fo = 1.5*D_fi

    # Fixed wall thickness
    t_3 = 0.0005  # Wall thickness [m]

    # Compute geometric properties
    dist = np.sqrt((x / 2 - e_2) ** 2 + (w / 2 - e_1) ** 2)  # Distance to fasteners [m]
    A_hole = 0.25 * np.pi * D_2 ** 2  # Area of a fastener hole [m²]
    A_shear = np.pi * D_2 * (t_2 + t_3)  # Combined shear area for plate and wall [m²]
    A_normal = np.pi * (D_fo ** 2 - D_fi ** 2) / 4  # Annular shear area under the fastener head [m²]

    # Compute forces on fasteners
    M_z_total = F_x * (w / 2 - e_1) + M_z
    M_x_total = F_z * (w / 2 - e_1) + M_x
    F_in_plane_y = F_y / nf - M_x_total/(x-2*e_1) + M_z_total/(w - 2*e_1) # In-plane force due to F_y [N]

    # Compute stresses
    sigma_normal = F_in_plane_y / A_normal  # Normal stress due to F_y [Pa]
    tau_t_2 = F_in_plane_y / (A_shear * t_2)  # Shear stress on the plate [Pa]
    tau_t_3 = F_in_plane_y / (A_normal * t_3)  # Shear stress on the wall [Pa]

    # Von Mises stresses
    sigma_VM_t_2 = np.sqrt(sigma_normal ** 2 + 3 * tau_t_2 ** 2)  # Von Mises stress on the plate [Pa]
    sigma_VM_t_3 = np.sqrt(sigma_normal ** 2 + 3 * tau_t_3 ** 2)  # Von Mises stress on the wall [Pa]

    # Margins of safety
    MS_t_2 = (yield_stress / sigma_VM_t_2) - 1  # Margin of safety for the plate
    MS_t_3 = (yield_stress / sigma_VM_t_3) - 1  # Margin of safety for the wall

    return MS_t_2, MS_t_3
