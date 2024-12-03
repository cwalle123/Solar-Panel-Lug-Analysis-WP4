# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:20:29 2024

@author: mkaha
"""

# flanges design
import numpy as np
# first order itration code in utf-8- bits
#check the flanges thickness under axial-load in the (y-z) plane
# to do check for oblique loads and add MS
# to do interpolate the width w, the height h , the e(location of the hole)
# and the D hole diametre 

fx = 1000 #transvers comp
fy = 1000 # axial comp
fz = 1000 # transvers comp
Mx = 1000
My = 1000
Mz = 1000
w = 30e-2 # base width assumed
D1 = 10e-2 # hole diameter assumed
W1D1 = w/D1
h = 40e-2  # flange profile height assumed 
e = 15e-3  # location of the hole also was assumed 
sigma_yield = 430e6
k1 = 0.5
k2 = 0.1
k3 = 0.25
  
def t_min_flage_axial_eq_6(w,D1,sigma_yield,fy):
    k1 = 0.5
    #f_axila = np.sqrt(fy**2+fz**2)
    #sigma_yield = 430e6
    #w =30 #np.arange(1,1000,1 )
    #D_1 = 15 #np.arange(1,1000,1 )
    #tot_param = f_axila/sigma_yield
    #find the smallest value for w ,D_1,t that is equal to tot_param 
    t = fy/(sigma_yield*(w-D1)*k1) 
    return t

t_min1 = t_min_flage_axial_eq_6(w,D1,sigma_yield , fy)*1000 # in mm
#checking the tension-shear


def t_min_flange_axial_shear(w,D1,sigma_yield,fy):
    k2 = 0.1
    t = fy/(sigma_yield*0.01*k2) 
    return t

t_min2 = t_min_flange_axial_shear(w, D1, sigma_yield , fy)*1000 # in mm

def Aav_over_Abearin_calculator(D1,h,e,w):
    c = 6/((4*D1/(h-e))+((2*D1)/(2*np.pi*0.5*w)))
    return c 


def t_min_flage_transverse(w,D_1,sigma_yield,fy):
    k3 = 0.25
    t = fy/(sigma_yield*0.01*k3) 
    return t

t_min3 = t_min_flage_transverse(w, D1, sigma_yield , fy)*1000 # in mm

#bearing check
def bearing_check(sigma_yield,D1,t):
    F = sigma_yield * D1 * t
    return F
f = bearing_check(sigma_yield , D1, t_min2)/1000

Pu = sigma_yield * k1 * D1 * t_min2 / 1000 # eq 6
pbru = sigma_yield  *k2 * D1 * t_min2 /1000# eq 7
Ptu = sigma_yield * k3 * D1 * t_min2/1000 # eq 10

#oblique loads
Ra = (fy / min(Pu, pbru))**1.6    # axial comp / min(Pu eq6 ,Pbru eq 7 )    
Rt = (np.sqrt(fx**2+fz**2) / Ptu)**1.6    
oblique = Ra + Rt 

#M.S. saftey margins

MS = (1/(oblique)**0.625) - 1
    
         
