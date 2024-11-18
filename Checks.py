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

