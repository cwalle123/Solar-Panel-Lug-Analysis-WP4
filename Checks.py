"""
checks.py contains verification of all different failure modes
each failure mode is given its own unique function that will be used
by calculator.py to verify an overall design
"""

# External imports
import numpy as np
import math

# elastic modulus [Pa], yield strength [Pa], shear strength [Pa], density [kg/m^3], thermal expansion coefficient [1/K]
materials = {
    'Al 7075 T6': (71.7e9, 430e6, .5*430e6, 2810, 2.36e-5)
}

def vectorized(func):
    """decorator used to easily vectorise any function"""
    return np.vectorize(func)

#Objective function
@vectorized
def mass(Density_Lug, Thickness_Lug_1, Thickness_Lug_2, w_Lug, h_Lug, Diameter_Lug_1, Diameter_Lug_2): #all inputs should be given in SI units
    """This function calculates the total mass of the lug by calculating and summing the volumes of
    the hinges and the backplate, and then multiplies the total volume by the density of the material"""
    HingeArea_Lug = w_Lug**(2) + (w_Lug**(2) * math.pi) / (8) - (Diameter_Lug_1**(2) * math.pi) / (4)  # m^2
    TotalHingeVolume_Lug = 2 * Thickness_Lug_1 * HingeArea_Lug  # m^3
    BackplateArea_Lug = (3 * h_Lug + 2 * Thickness_Lug_1) * w_Lug - Diameter_Lug_2**(2) * math.pi  # m^2
    BackplateVolume_Lug = Thickness_Lug_2 * BackplateArea_Lug  # m^3
    TotalVolume = TotalHingeVolume_Lug + BackplateVolume_Lug  # m^3
    return Density_Lug * TotalVolume  # kg
