"""
Created on Wed Nov 20 14:32:43 2024

@author: lauke
"""

#to function properly, the instruments need be between -23.15◦C to +26.85◦C


#imports
import numpy as np
from math import pi
'''inputs'''

'''dummy values'''
#Stiffness area of fastener
stiffnessArea_fastener = 1
#force ratio to be determined at fastening section
force_ratio = 1
thermalExp_coef_fastener = 1
thermalExp_coef_clamped = 1
E = 1


def thermal_load_fastener(thermalExp_coef_fastener, thermalExp_coef_clamped, temp, E_fastener, E_backplate, A_stiffness, temp_ref, temp_working_max, temp_working_min, area_fastener, stress_max_fastener, d_outer, d_inner, distance_thickness, d_thread, shank_length, d_shank):
    """calculate the force induced by thermal expansion
        -Material of both fastener and clamped is needed
        -"""
        
    delta_a = 4 * distance_thickness / (E_backplate * pi (d_outer**2 + d_inner**2))
    
    delta_b = 1/E_fastener * (0.4 * d_thread / (d_outer ** 2 * pi / 4) + 0.4 * d_thread / (d_inner ** 2 * pi / 4) + shank_length / (d_shank **2 * pi / 4) ) + 0.4 * d_thread / (E_fastener * d_outer ** 2 * pi / 4)
    
    force_ratio = delta_b/(delta_a + delta_b)
    
    min_temp = temp_working_min - temp_ref
    
    max_temp = temp_working_max - temp_ref
    
    min_thermal_load = (thermalExp_coef_fastener-thermalExp_coef_clamped) * min_temp * E * A_stiffness * (1-force_ratio)
    
    max_thermal_load = (thermalExp_coef_fastener-thermalExp_coef_clamped) * max_temp * E * A_stiffness * (1-force_ratio)
    
    safety_factor = stress_max_fastener/(max([min_thermal_load,max_thermal_load])/area_fastener)-1
    
    #we only care about maximum so we take max
    
    return np.vectorize(safety_factor)
    
