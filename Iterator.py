"""
iterator.py is the primary logic loop of the software, which iterates upon the design variables
in order to converge on the most optimal design
given a certain initial set of design variables, it calls the calculator and evaluates the
out-coming safety factors and mass to adjust the design variables and run the calculator again
employing a certain optimisation strategy, iterator.py tries to converge
on the optimal design configuration
"""

# External imports
import numpy as np

# Internal imports
from calculator import *
from checks import *


