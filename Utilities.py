# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 19:59:25 2024

@author: herbi
"""
import numpy as np

epsilon = 1e-5 #UV cutoff
R = 1 #AdS radius
R_EFF = R-epsilon #Radius at which the boundary points live
S_NUMBER = int(1000) #Resolution of curve.

def dot_product(x1, x2):
    #Embedded
    return -x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] - x1[3]*x2[3]


def boundary_to_global(theta):
    return polar_to_global([R_EFF,theta])

def polar_to_global(polars):
    r = polars[0]
    phi = polars[1]
    
    coshp = (1+r**2)/(1-r**2)
    sinhp = 2*r/(1-r**2)
    
    X0 = R*coshp
    X1 = R*sinhp*np.cos(phi)
    X2 = R*sinhp*np.sin(phi)
    X3 = 0 
    
    return np.array([X0,X1,X2,X3])

def find_distance_global(XA,XB):
    zet = -dot_product(XA, XB)/R**2
    return R*np.arccosh(zet)

def find_distance_polars(thet1,thet2):
    XA = boundary_to_global(thet1)
    XB = boundary_to_global(thet2)
    return find_distance_global(XA, XB)

def cartesian_to_global(x,y):
    r = np.sqrt(x**2+y**2)
    
    coshp = (1+r**2)/(1-r**2)
    sinhp = 2*r/(1-r**2)
    
    X0 = R*coshp
    X1 = R*sinhp * x/r
    X2 = R*sinhp * y/r
    X3 = 0

    return np.array([X0,X1,X2,X3])

def global_to_cartesian(X):
    X0 = X[0]/R
    r = np.sqrt((X0-1)/(X0+1))
    sinhp = (2*r)/(1-r**2)
    
    
    
    cos_phi = X[1]/(R*sinhp)   
    sin_phi = X[2]/(R*sinhp)
    
    return r*cos_phi,r*sin_phi


def get_geodesic_global_to_cartesian(XA,XB):
    SB = find_distance_global(XA, XB)
    
    s_arr = np.linspace(0, SB, S_NUMBER)
    
    xs = []
    
    q = (1/np.sinh(SB/R)) * (XB-np.cosh(SB/R)*XA)
    
    for s in s_arr:
        x = XA*np.cosh(s/R) + q*np.sinh(s/R)
        xs.append(x)
        
    
    x_arr = []
    y_arr = []
    for point in xs:
        x,y = global_to_cartesian(point)
        x_arr.append(x)
        y_arr.append(y)
        
    return x_arr,y_arr

def get_geodesic_global_to_global(XA,XB):
    SB = find_distance_global(XA, XB)
    
    s_arr = np.linspace(0, SB, S_NUMBER)
    
    xs = []
    
    q = (1/np.sinh(SB/R)) * (XB-np.cosh(SB/R)*XA)
    
    for s in s_arr:
        x = XA*np.cosh(s/R) + q*np.sinh(s/R)
        xs.append(x)
           
    return np.array(xs)



def get_geodesic_polar_to_cartesian(polarA,polarB):
    XA = polar_to_global(polarA)
    XB = polar_to_global(polarB)

    return get_geodesic_global_to_cartesian(XA, XB)

def get_geodesic_cartesian_to_cartesian(cartesianA,cartesianB):
    XA = cartesian_to_global(cartesianA[0], cartesianA[1])
    XB = cartesian_to_global(cartesianB[0],cartesianB[1])
    return get_geodesic_global_to_cartesian(XA, XB)

def combinatoric_matching(mapping_array, template_array):
    flattened_template = [item for sublist in template_array for item in sublist]
    value_to_var = {v: flattened_template[i] for i, v in enumerate(
        sorted(set(v for pair in mapping_array for v in pair)))}
    final_matrix = [[value_to_var[pair[0]],
                     value_to_var[pair[1]]] for pair in mapping_array]
    
    return final_matrix