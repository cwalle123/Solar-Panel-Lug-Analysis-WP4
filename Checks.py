"""
checks.py contains verification of all different failure modes
each failure mode is given its own unique function that will be used
by calculator.py to verify an overall design
"""

# External imports
import numpy as np


def vectorized(func):
    """decorator used to easily vectorise any function"""
    return np.vectorize(func)

