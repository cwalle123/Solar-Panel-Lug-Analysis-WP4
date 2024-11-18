"""
calculator.py verifies the validity of a single design, i.e. a single set of design parameters.
it runs instances of every check function and compiles the resulting safety margins to either
print them out or pass them on to the iterator.py
it also computes the mass of the design
"""

# External imports
import numpy as np

# Internal imports
from checks import *

