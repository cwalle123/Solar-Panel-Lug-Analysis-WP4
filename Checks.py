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
    all of the lengths are actually greater than zero or stuff like that"""
    ...
    # first is name of the verification check, then the True/False on whether the check is passed
    return str, bool


@vectorized
def positivity_check(design_vector):
    """Ensures no length is given as negative, as that would be nonsensical"""
    if np.all(design_vector > 0):
        return True
    else:
        return False


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



