"""
checks.py contains verification of all different failure modes
each failure mode is given its own unique function that will be used
by calculator.py to verify an overall design
"""

# External imports
import numpy as np


# elastic modulus [Pa], yield strength [Pa], shear strength [Pa], density [kg/m^3], thermal expansion coefficient [1/K]
materials = {
    'Al 7075 T6': (71.7e9, 430e6, .5*430e6, 2810, 2.36e-5)
}


def vectorized(func):
    """decorator used to easily vectorise any function"""
    return np.vectorize(func)

#Objective function
def mass(Density_Lug,Thickness_Lug_1,Thickness_Lug_2,w_Lug,h_Lug,Diameter_Lug_1,Diameter_Lug_2): #Inputs for the function.
    HingeArea_Lug = w_Lug**(2) + (w_Lug**(2) * math.pi) / (8) - (Diameter_Lug_1**(2) * math.pi) / (4) #This calculates the area of one hinge with the assumption that it is made up of a half circle and a square attached to each other on their flat edges with a smaller hole in the middle of where the two edges meet. Both edges are the same lenght "w_Lug" and the diameter of the smaller hole is "Diameter_Lug_1".
    TotalHingeVolume_Lug = 2 * Thickness_Lug_1 * HingeArea_Lug #This calculates the volume of a hinge and multiplies it by 2 as there are two hinges.
    BackplateArea_Lug = (3 * h_Lug + 2 * Thickness_Lug_1) * w_Lug - Diameter_Lug_2**(2) * math.pi #This calculates the area of the backplate with the assumption that the total width is "(3 * h_Lug + 2 * Thickness_Lug_1)".
    BackplateVolume_Lug = Thickness_Lug_2 * BackplateArea_Lug #This calculates the volume of the backplate
    TotalVolume = TotalHingeVolume_Lug + BackplateVolume_Lug #This calculates the total volume by adding the volume of the 2 hinges and the backplate
    return Density_Lug * TotalVolume #This returns the total mass of the lug (= density * volume)
