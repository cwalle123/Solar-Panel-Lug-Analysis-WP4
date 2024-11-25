# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:32:43 2024

@author: lauke
"""

#to function properly, the instruments need be between -23.15◦C to +26.85◦C


#imports
import numpy as np
'''inputs'''


#Stiffness area of fastener
stiffnessArea_fastener = 1
force_ratio = 1
thermalExp_coef_fastener = 1
thermalExp_coef_clamped = 1
E = 1


def thermal_load_fastener(thermalExp_coef_fastener, thermalExp_coef_clamped, temp, E, A_stiffness, force_ratio, temp_ref, temp_working_max, temp_working_min, area_fastener, stress_max_fastener):
    """calculate the force induced by thermal expansion
        -Material of both fastener and clamped is needed
        -"""
    
    min_temp = temp_working_min - temp_ref
    
    max_temp = temp_working_max - temp_ref
    
    min_thermal_load = (thermalExp_coef_fastener-thermalExp_coef_clamped) * min_temp * E * A_stiffness * (1-force_ratio)
    
    max_thermal_load = (thermalExp_coef_fastener-thermalExp_coef_clamped) * max_temp * E * A_stiffness * (1-force_ratio)
    
    safety_factor = stress_max_fastener/(max([min_thermal_load,max_thermal_load])/area_fastener)-1
    
    #we only care about maximum so we take max
    
    return np.vectorize(safety_factor)
    
