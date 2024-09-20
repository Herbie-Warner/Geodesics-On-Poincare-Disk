"""
Plotting and calculating lengths of geodesics on the Poincare disk given some
subsystem definitions.

Herbie Warner 10/07/2024
"""

import numpy as np
import matplotlib.pyplot as plt

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

def plot_entanglement_wedge_from_global(XA0,XA1,XB0,XB1,ax,name,colour):

    geoB = get_geodesic_global_to_global(XB0, XB1)
    

    
    minDist0 = np.inf
    minDist1 = np.inf
    for point in geoB:
        dist = find_distance_global(XA0, point)
        if dist < minDist0:
            minDist0 = dist
            pointMin0 = point
           
        dist = find_distance_global(XA1, point)
        if dist < minDist1:
            minDist1 = dist
            pointMin1 = point
  
        

    if minDist1+minDist0 < find_distance_global(XA0, XA1):
        x0,y0 = get_geodesic_global_to_cartesian(XA0, pointMin0) 
        ax.plot(x0,y0,color=colour,label=name)
        x1,y1 = get_geodesic_global_to_cartesian(XA1, pointMin1)
        ax.plot(x1,y1,color=colour)

        
    else:
        x,y =  get_geodesic_global_to_cartesian(XA0, XA1)
        ax.plot(x,y,color=colour,label=name)
       
    
def plot_entanglement_wedge_from_polars(pointA0,pointA1,pointB0,pointB1,ax,name,colour):
    XA0 = polar_to_global(pointA0)
    XA1 = polar_to_global(pointA1)
    XB0 = polar_to_global(pointB0)
    XB1 = polar_to_global(pointB1)
    
    plot_entanglement_wedge_from_global(XA0, XA1, XB0, XB1,ax,name,colour)

        
    
        
def plot_entanglement_wedge_from_cartesian(pointA0,pointA1,pointB0,pointB1,ax,name,colour):
    XA0 = cartesian_to_global(pointA0[0], pointA0[1])
    XA1 = cartesian_to_global(pointA1[0], pointA1[1])
    XB0 = cartesian_to_global(pointB0[0], pointB0[1])
    XB1 = cartesian_to_global(pointB1[0], pointB1[1])
    plot_entanglement_wedge_from_global(XA0, XA1, XB0, XB1,ax,name,colour)


def get_EW_length_from_polars(pointA0,pointA1,pointB0,pointB1):
    XA0 = polar_to_global(pointA0)
    XA1 = polar_to_global(pointA1)
    XB0 = polar_to_global(pointB0)
    XB1 = polar_to_global(pointB1)
    geoB = get_geodesic_global_to_global(XB0, XB1)
    

    minDist0 = np.inf
    minDist1 = np.inf
    for point in geoB:
        dist = find_distance_global(XA0, point)
        if dist < minDist0:
            minDist0 = dist
           
        dist = find_distance_global(XA1, point)
        if dist < minDist1:
            minDist1 = dist
  
        
    return min(minDist1+minDist0 ,find_distance_global(XA0, XA1))

    

def get_EW_length_from_polars_E_AD(pointA0,pointA1,pointB0,pointB1,pointC0,pointC1):
    XA0 = polar_to_global(pointA0)
    XA1 = polar_to_global(pointA1)
    XB0 = polar_to_global(pointB0)
    XB1 = polar_to_global(pointB1)
    XC0 = polar_to_global(pointC0)
    XC1 = polar_to_global(pointC1)
    
    
    geoB = get_geodesic_global_to_global(XB0, XB1)
    
    minDistB0 = np.inf
    minDistB1 = np.inf
    
    for point in geoB:
        dist = find_distance_global(XA0, point)
        if dist < minDistB0:
            minDistB0 = dist
           
        dist = find_distance_global(XA1, point)
        if dist < minDistB1:
            minDistB1 = dist
            
    tot_B_dist = minDistB1 + minDistB0
    
    geoC = get_geodesic_global_to_global(XC0, XC1)
    
    minDistC0 = np.inf
    minDistC1 = np.inf
    
    for point in geoC:
        dist = find_distance_global(XA0, point)
        if dist < minDistC0:
            minDistC0 = dist
           
        dist = find_distance_global(XA1, point)
        if dist < minDistC1:
            minDistC1 = dist
            
    tot_C_dist = minDistC1 + minDistC0
    
    
    
    minDistBC = np.inf

    
    for pointB in geoB:
        for pointC in geoC:
            dist = find_distance_global(pointB, pointC)
            if dist < minDistBC:
                minDistBC = dist

                
              
                
    total_o_length = min(minDistC1,minDistB1) + min(minDistC0,minDistB0) + minDistBC
    
    
    SA = find_distance_global(XA0, XA1)
    
    minimum = min(SA,total_o_length,tot_B_dist,tot_C_dist)
    return minimum

def find_entanglement_wedge_dist_E_AD_polars(pointA0,pointA1,pointB0,pointB1,pointC0,pointC1):
    XA0 = polar_to_global(pointA0)
    XA1 = polar_to_global(pointA1)
    XB0 = polar_to_global(pointB0)
    XB1 = polar_to_global(pointB1)
    XC0 = polar_to_global(pointC0)
    XC1 = polar_to_global(pointC1)
    
    
    geoB = get_geodesic_global_to_global(XB0, XB1)
    
    minDistB0 = np.inf
    minDistB1 = np.inf
    pointMinB0 = None
    pointMinB1 = None
    
    for point in geoB:
        dist = find_distance_global(XA0, point)
        if dist < minDistB0:
            minDistB0 = dist
            pointMinB0 = point
           
        dist = find_distance_global(XA1, point)
        if dist < minDistB1:
            minDistB1 = dist
            pointMinB1 = point
            
    tot_B_dist = minDistB1 + minDistB0
    
    geoC = get_geodesic_global_to_global(XC0, XC1)
    
    minDistC0 = np.inf
    minDistC1 = np.inf
    pointMinC0 = None
    pointMinC1 = None
    
    for point in geoC:
        dist = find_distance_global(XA0, point)
        if dist < minDistC0:
            minDistC0 = dist
            pointMinC0 = point
           
        dist = find_distance_global(XA1, point)
        if dist < minDistC1:
            minDistC1 = dist
            pointMinC1 = point
            
    tot_C_dist = minDistC1 + minDistC0
    
    
    
    minDistBC = np.inf
    pointMinB2 = None
    pointMinC2 = None
    
    for pointB in geoB:
        for pointC in geoC:
            dist = find_distance_global(pointB, pointC)
            if dist < minDistBC:
                minDistBC = dist
                pointMinB2 = pointB
                pointMinC2 = pointC
                
                
                
    total_o_length = min(minDistC1,minDistB1) + min(minDistC0,minDistB0) + minDistBC
    
    
    SA = find_distance_global(XA0, XA1)
    
    minimum = min(SA,total_o_length,tot_B_dist,tot_C_dist)
    
   
    return minimum
    
    
        

def plot_entanglement_wedge_from_polars_E_AD(pointA0,pointA1,pointB0,pointB1,pointC0,pointC1,ax,name,colour):
    XA0 = polar_to_global(pointA0)
    XA1 = polar_to_global(pointA1)
    XB0 = polar_to_global(pointB0)
    XB1 = polar_to_global(pointB1)
    XC0 = polar_to_global(pointC0)
    XC1 = polar_to_global(pointC1)
    
    
    geoB = get_geodesic_global_to_global(XB0, XB1)
    
    minDistB0 = np.inf
    minDistB1 = np.inf
    pointMinB0 = None
    pointMinB1 = None
    
    for point in geoB:
        dist = find_distance_global(XA0, point)
        if dist < minDistB0:
            minDistB0 = dist
            pointMinB0 = point
           
        dist = find_distance_global(XA1, point)
        if dist < minDistB1:
            minDistB1 = dist
            pointMinB1 = point
            
    tot_B_dist = minDistB1 + minDistB0
    
    geoC = get_geodesic_global_to_global(XC0, XC1)
    
    minDistC0 = np.inf
    minDistC1 = np.inf
    pointMinC0 = None
    pointMinC1 = None
    
    for point in geoC:
        dist = find_distance_global(XA0, point)
        if dist < minDistC0:
            minDistC0 = dist
            pointMinC0 = point
           
        dist = find_distance_global(XA1, point)
        if dist < minDistC1:
            minDistC1 = dist
            pointMinC1 = point
            
    tot_C_dist = minDistC1 + minDistC0
    
    
    
    minDistBC = np.inf
    pointMinB2 = None
    pointMinC2 = None
    
    for pointB in geoB:
        for pointC in geoC:
            dist = find_distance_global(pointB, pointC)
            if dist < minDistBC:
                minDistBC = dist
                pointMinB2 = pointB
                pointMinC2 = pointC
                
                
                
    total_o_length = min(minDistC1,minDistB1) + min(minDistC0,minDistB0) + minDistBC
    
    
    SA = find_distance_global(XA0, XA1)
    
    minimum = min(SA,total_o_length,tot_B_dist,tot_C_dist)
    
    print(minDistB0+minDistB1)
    print(minDistC0+minDistC1)
    print(SA)
    print(total_o_length)
    
    if minimum == SA:
        x,y = get_geodesic_global_to_cartesian(XA0, XA1)
        ax.plot(x,y,color=colour,label=name)
    elif minimum == tot_B_dist:
        x0,y0 = get_geodesic_global_to_cartesian(XA0, pointMinB0) 
        ax.plot(x0,y0,color=colour,label=name)
        x1,y1 = get_geodesic_global_to_cartesian(XA1, pointMinB1)
        ax.plot(x1,y1,color=colour)
    elif minimum == tot_C_dist:
        x0,y0 = get_geodesic_global_to_cartesian(XA0, pointMinC0) 
        ax.plot(x0,y0,color=colour,label=name)
        x1,y1 = get_geodesic_global_to_cartesian(XA1, pointMinC1)
        ax.plot(x1,y1,color=colour)
        
    elif minimum == total_o_length:
        
        if minDistC1 < minDistB1:
            point1 = pointMinC1
        else:
            point1 = pointMinB1
            
        
        if minDistC0 < minDistB0:
            point0 = pointMinC0
        else:
            point0 = pointMinB0
        
        x0,y0 = get_geodesic_global_to_cartesian(XA0, point0) 
        ax.plot(x0,y0,color=colour,label=name)
        x1,y1 = get_geodesic_global_to_cartesian(XA1, point1)
        ax.plot(x1,y1,color=colour)
        x1,y1 = get_geodesic_global_to_cartesian(pointMinB2, pointMinC2)
        ax.plot(x1,y1,color=colour)
        
        
        


def plot_two_system():
    theta = np.linspace(0, 2 * np.pi, 200)
    #max_rho = R * np.sinh(1 / epsilon)
    
    pointA0 = [1-epsilon, 0.7*np.pi/3]
    pointA1 = [1-epsilon, 2.3*np.pi/3]
    
    pointB0 = [1-epsilon,  np.pi]
    pointB1 = [1-epsilon, 2 * np.pi]

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    
    
    x,y = get_geodesic_polar_to_cartesian(pointA0, pointA1)
    ax.plot(x,y,label="GEOA")
    
    x,y = get_geodesic_polar_to_cartesian(pointB0, pointB1)
    ax.plot(x,y,label="GEOB")
    
    
    plot_entanglement_wedge_from_polars(pointA0, pointA1, pointB0, pointB1, ax, "EW(A:B)","y")
    

    #x,y = get_entanglement_wedge_from_polars(pointA0, pointA1, pointB0, pointB1)
    #x,y = get_geodesic_global_to_cartesian(globalA0, global_min)
    #plt.plot(x,y)
    ax.set_title('Geodesics in the Poincaré Disk')
    ax.legend(loc='upper right') 
    plt.show()
    
def polar_to_cartesian(polar):
    return [polar[0]*np.cos(polar[1]),polar[0]*np.sin(polar[1])]

def cartesian_to_polar(cartesian):
    r = np.sqrt(cartesian[0]**2+cartesian[1]**2)
    theta = np.arctan(cartesian[1]/cartesian[0])
    return [r,theta]
    

def compute_inequality(params):
    pointA0 = [1-epsilon, params[0]]
    pointA1 = [1-epsilon, params[1]]
    
    pointB0 = [1-epsilon,  params[2]]
    pointB1 = [1-epsilon, params[3]]
    
    pointC0 = [1-epsilon, params[4]]
    pointC1 = [1-epsilon, params[5]]
    
    EW_A_CD = get_EW_length_from_polars(pointA0, pointA1, pointB0, pointB1)
    EW_A_BD = get_EW_length_from_polars(pointA0, pointA1, pointC0, pointC1)
    EW_A_D = get_EW_length_from_polars_E_AD(pointA0, pointA1, pointB0, pointB1, pointC0, pointC1)
    S_A = find_distance_global(polar_to_global(pointA0),polar_to_global(pointA1))
    
    ineq = EW_A_CD+EW_A_BD-EW_A_D-S_A
    return ineq
    

def plot_three_system(params):
    theta = np.linspace(0, 2 * np.pi, 200)
    #max_rho = R * np.sinh(1 / epsilon)
    
    pointA0 = [1-epsilon, params[0]]
    pointA1 = [1-epsilon, params[1]]
    
    pointB0 = [1-epsilon,  params[2]]
    pointB1 = [1-epsilon, params[3]]
    
    pointC0 = [1-epsilon, params[4]]
    pointC1 = [1-epsilon, params[5]]

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
    C_angles = np.linspace(pointC0[1], pointC1[1], 100)
    ax.plot(np.cos(C_angles),np.sin(C_angles), label="C", linewidth=5)
    
    
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    
    
    x,y = get_geodesic_polar_to_cartesian(pointA0, pointA1)
    ax.plot(x,y,label="GEOA")
    
    x,y = get_geodesic_polar_to_cartesian(pointB0, pointB1)
    ax.plot(x,y,label="GEOB")
    
    x,y = get_geodesic_polar_to_cartesian(pointC0, pointC1)
    ax.plot(x,y,label="GEOC")
    
    
    plot_entanglement_wedge_from_polars_E_AD(pointA0,pointA1,pointB0,pointB1,pointC0,pointC1,ax,"EW(A:CD)","g")
    
    plot_entanglement_wedge_from_polars(pointA0, pointA1, pointB0, pointB1, ax, "EW(A:BD)","y")
    plot_entanglement_wedge_from_polars(pointA0, pointA1, pointC0, pointC1, ax, "EW(A:CD)","r")
   
   
    
    EW_A_CD = get_EW_length_from_polars(pointA0, pointA1, pointB0, pointB1)
    EW_A_BD = get_EW_length_from_polars(pointA0, pointA1, pointC0, pointC1)
    EW_A_D = get_EW_length_from_polars_E_AD(pointA0, pointA1, pointB0, pointB1, pointC0, pointC1)
    S_A = find_distance_global(polar_to_global(pointA0),polar_to_global(pointA1))
    
    print("EW(A:CD): ", EW_A_CD)
    print("EW(A:BD): ", EW_A_BD)
    print("EW(A:D):  ", EW_A_D)
    print("S Of (A): ", S_A)
    print()
    print("In: ",EW_A_CD+EW_A_BD-EW_A_D-S_A)

    ax.set_title('Geodesics in the Poincaré Disk')
    ax.legend(loc='upper right') 
    plt.savefig("EW_sys.png",dpi=600)
    plt.show()
    
    
def plot_bipartite_wedge():
    """
    I force connected E(AB) for graphical purposes

    Returns
    -------
    None.

    """
    

    theta = np.linspace(0, 2 * np.pi, 200)

    pointA0 = [R_EFF,-np.pi/3]
    pointA1 = [R_EFF,np.pi/3]

    pointB0 = [R_EFF,2*np.pi/3]
    pointB1 = [R_EFF, 4* np.pi/3]

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
      
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    
    x0,y0 = get_geodesic_polar_to_cartesian(pointA1, pointB0)
    ax.plot(x0,y0)
    
    x1,y1 = get_geodesic_polar_to_cartesian(pointA0, pointB1)
    ax.plot(x1,y1)
    
    
    minPoint0 = None
    minDist0 = np.inf
    
    minPoint1 = None

    for i in range(len(x0)):
        point0 = cartesian_to_global(x0[i], y0[i])
        for j in range(len(x1)):
            point1 = cartesian_to_global(x1[i], y1[i])
            dist = find_distance_global(point0, point1)
            if dist < minDist0:
                minPoint0 = point0
                minPoint1 = point1
                minDist0 = dist
                
                
    ewx,ewy = get_geodesic_global_to_cartesian(minPoint0, minPoint1)
    ax.plot(ewx,ewy,color="g")
    

    plt.show()
    
def plot_tripartite_wedge():


    theta = np.linspace(0, 2 * np.pi, 200)

    pointA0 = [R_EFF,np.pi/3]
    pointA1 = [R_EFF,2*np.pi/3]

    pointB0 = [R_EFF,3*np.pi/3]
    pointB1 = [R_EFF, 4* np.pi/3]
    
    
    pointC0 = [R_EFF,5*np.pi/3]
    pointC1 = [R_EFF, 6* np.pi/3]

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
      
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    
    C_angles = np.linspace(pointC0[1], pointC1[1], 100)
    ax.plot(np.cos(C_angles),np.sin(C_angles), label="C", linewidth=5)
    
    x0,y0 = get_geodesic_polar_to_cartesian(pointA1, pointB0)
    ax.plot(x0,y0,color='k')
    
    x1,y1 = get_geodesic_polar_to_cartesian(pointA0, pointC1)
    ax.plot(x1,y1,color='k')
    
    x2,y2 = get_geodesic_polar_to_cartesian(pointB1, pointC0)
    ax.plot(x2,y2,color='k')
    
    
    
    minPoint0 = None
    minDist0 = np.inf
    
    minPoint1 = None

    for i in range(len(x0)):
        point0 = cartesian_to_global(x0[i], y0[i])
        for j in range(len(x1)):
            point1 = cartesian_to_global(x1[i], y1[i])
            dist = find_distance_global(point0, point1)
            if dist < minDist0:
                minPoint0 = point0
                minPoint1 = point1
                minDist0 = dist
                
                
    ewx,ewy = get_geodesic_global_to_cartesian(minPoint0, minPoint1)
    ax.plot(ewx,ewy,label="EW")
    
    
    #ax.set_title('Geodesics in the Poincaré Disk')
    #ax.legend(loc='upper right') 
    plt.show()


    
    
def calculate_lengths(params):
    """
    J(A|B)-J(A|C)+J(A|BC)
    """
    S_A = find_distance_polars(params[0], params[1])
    S_B = find_distance_polars(params[2], params[3])
    S_C = find_distance_polars(params[4], params[5])
    
    S_AB = S_A + S_B
    S_BC = find_distance_polars(params[2], params[5])
    S_AC = find_distance_polars(params[5], params[0]) + find_distance_polars(params[1], params[4])
    
    S_ABC = find_distance_polars(params[5], params[0]) + find_distance_polars(params[1], params[2])
    
    I3 = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    J_A_B = 0
    
    
def plot_PQ(pointA0,pointA1,pointB0,pointB1,P,Q):
    theta = np.linspace(0, 2 * np.pi, 200)

    
    

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
      
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    

    P_glob = polar_to_global(P)
    Q_glob = polar_to_global(Q)
    
    A0_glob = polar_to_global(pointA0)
    A1_glob = polar_to_global(pointA1)
    B0_glob = polar_to_global(pointB0)
    B1_glob = polar_to_global(pointB1)

    ax.scatter(polar_to_cartesian(P)[0],polar_to_cartesian(P)[1],label="P")
    ax.scatter(polar_to_cartesian(Q)[0],polar_to_cartesian(Q)[1],label="Q")
    
    
    x,y = get_geodesic_polar_to_cartesian(pointA0, P)
    ax.plot(x,y,color='k')
   
    x,y = get_geodesic_polar_to_cartesian(pointA1,Q)
    ax.plot(x,y,color='g')
    x,y = get_geodesic_polar_to_cartesian(pointB1, P)
    ax.plot(x,y,color='b')
    x,y = get_geodesic_polar_to_cartesian(pointB0, Q)
    ax.plot(x,y,color='y')
    x,y = get_geodesic_polar_to_cartesian(P, Q)
    ax.plot(x,y,color='r')
    
    x,y = get_geodesic_polar_to_cartesian(pointB0, pointB1)
    ax.plot(x,y,color='k',linestyle='--')
  
    
    
    #plot_entanglement_wedge_from_polars(pointA0, pointA1, pointB0, pointB1, ax, 'EW(A:C)', 'g')
  
    
    #ax.set_title('Geodesics in the Poincaré Disk')
    ax.legend(loc='upper right') 
    plt.show()
    
    
    
def length_function(coords,pointA0,pointA1,pointB0,pointB1,renyi):
    
    
    
    r1, theta1, r2, theta2 = coords
    
    P = [r1,theta1]
    Q = [r2,theta2]
 
 
    
        
    P_glob = polar_to_global(P)
    Q_glob = polar_to_global(Q)
    
    A0_glob = polar_to_global(pointA0)
    A1_glob = polar_to_global(pointA1)
    B0_glob = polar_to_global(pointB0)
    B1_glob = polar_to_global(pointB1)
    
    lengthA0P = find_distance_global(A0_glob, P_glob)
    lengthA1Q = find_distance_global(A1_glob, Q_glob)
    lengthB1P = find_distance_global(B1_glob, P_glob)
    lengthB0Q = find_distance_global(B0_glob, Q_glob)
    lengthPQ = find_distance_global(P_glob, Q_glob)
    lengthB = find_distance_global(B0_glob, B1_glob)
    
        
    summed = lengthA0P + lengthA1Q + ((renyi/(renyi-1)) *(lengthB1P+lengthPQ+lengthB0Q-lengthB))
    
    return summed


def minimisation():
    
    pointA0 = [R_EFF,2*np.pi/3]
    pointA1 = [R_EFF,4*np.pi/3]

    pointB1 = [R_EFF,np.pi/3]
    pointB0 = [R_EFF, - np.pi/3]
    
    renyi = 1.01
    
    
    bounds = [(0, R_EFF), (0, 2 * np.pi), (0, R_EFF), (0, 2 * np.pi)]

    res = differential_evolution(length_function, bounds, args=(pointA0,pointA1,pointB0,pointB1,renyi),disp=False, strategy='best1bin')
  
    optimized_coords = res.x
 
    
    P = [optimized_coords[0],optimized_coords[1]]
    Q = [optimized_coords[2],optimized_coords[3]]
    
    print(res.fun)
    print(get_EW_length_from_polars(pointA0, pointA1, pointB0, pointB1))
    
    
 
    
    #print(optimized_coords)
   
    plot_PQ(pointA0, pointA1, pointB0, pointB1, P, Q)
    
    #geoBX,geoBY = get_geodesic_polar_to_cartesian(pointB0, pointB1)   
    #P = cartesian_to_polar([geoBX[250],geoBY[250]])
    #Q = cartesian_to_polar([geoBX[200],geoBY[200]])
    #print(length_function([P[0],P[1],Q[0],Q[1]], pointA0, pointA1, pointB0, pointB1, renyi))
    
    
def convert_angle_to_latex_angle(polar):
    thet = polar[1]
    thet *= 180/np.pi
    #thet += 180
    return [polar[0]*1.5,thet]
    
def generate_points_for_graph_latex(params):
    theta = np.linspace(0, 2 * np.pi, 200)
    #max_rho = R * np.sinh(1 / epsilon)
    
    pointA0 = [1-epsilon, params[0]]
    pointA1 = [1-epsilon, params[1]]
    
    pointB0 = [1-epsilon,  params[2]]
    pointB1 = [1-epsilon, params[3]]
    
    pointC0 = [1-epsilon, params[4]]
    pointC1 = [1-epsilon, params[5]]

    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')
    

    A_angles = np.linspace(pointA0[1], pointA1[1], 100)
    ax.plot(np.cos(A_angles),np.sin(A_angles), label="A", linewidth=5)
    
    C_angles = np.linspace(pointC0[1], pointC1[1], 100)
    ax.plot(np.cos(C_angles),np.sin(C_angles), label="C", linewidth=5)
    
    
    B_angles = np.linspace(pointB0[1], pointB1[1], 100)
    ax.plot(np.cos(B_angles),np.sin(B_angles), label="B", linewidth=5)
    
    
    xA,yA = get_geodesic_polar_to_cartesian(pointA0, pointA1)
    #ax.plot(xA,yA,label="GEOA",color='k',linewidth=10)
    
    xB,yB = get_geodesic_polar_to_cartesian(pointB0, pointB1)
    #ax.plot(xB,yB,label="GEOB",color='k',linewidth=10)
    
    xC,yC = get_geodesic_polar_to_cartesian(pointC0, pointC1)
    #ax.plot(xC,yC,label="GEOC",color='k',linewidth=10)
    
    
    
    
    A0Glob = polar_to_global(pointA0)
    A1Glob = polar_to_global(pointA1)
        
    B0Glob = polar_to_global(pointB0)
    B1Glob = polar_to_global(pointB1)
        
    C0Glob = polar_to_global(pointC0)
    C1Glob = polar_to_global(pointC1)
    
    SA =  find_distance_global(A0Glob, A1Glob)
    SB =  find_distance_global(B0Glob, B1Glob)
    SC =  find_distance_global(C0Glob, C1Glob)
    
    distAB_connnected = find_distance_global(A0Glob, B1Glob)+find_distance_global(A1Glob, B0Glob)
    distAB_disconnected =SA+SB  
    distAC_connnected = find_distance_global(A0Glob, C1Glob)+find_distance_global(A1Glob, C0Glob)
    distAC_disconnected =SA+SC 
    distBC_connnected = find_distance_global(B0Glob, C1Glob)+find_distance_global(B1Glob, C0Glob)
    distBC_disconnected =SB+SC 
    
    S_AB = min(distAB_connnected,distAB_disconnected)
    S_AC = min(distAC_connnected,distAC_disconnected)
    S_BC = min(distBC_connnected,distBC_disconnected)
    
    
    
    
    case1 = SA+SB+SC
    case2 = SA+S_BC
    case3 = SB+S_AC
    case4 = SC+S_AB
    case5 =  find_distance_global(A1Glob, B0Glob) + find_distance_global(B1Glob, C0Glob)+find_distance_global(C1Glob, A0Glob)
    S_ABC = min(case1,case2,case3,case4,case5)
    
    """
    
   
    if distAB_connnected < distAB_disconnected:
        x,y = get_geodesic_global_to_cartesian(A0Glob, B1Glob)
        ax.plot(x,y,color="g",linewidth=8,label=r'$\mathcal{E}(AB)$')
        x,y = get_geodesic_global_to_cartesian(A1Glob, B0Glob)
        ax.plot(x,y,color="g",linewidth=8)
    else:
        ax.plot(xA,yA,color="g",linewidth=8)
        ax.plot(xB,yB,color="g",linewidth=8)
     
    distAC_connnected = find_distance_global(A0Glob, C1Glob)+find_distance_global(A1Glob, C0Glob)
    distAC_disconnected =SA+SC 
    if distAC_connnected < distAC_disconnected:
        x,y = get_geodesic_global_to_cartesian(A0Glob, C1Glob)
        ax.plot(x,y,color="r",linewidth=8,label=r'$\mathcal{E}(AC)$')
        x,y = get_geodesic_global_to_cartesian(A1Glob, C0Glob)
        ax.plot(x,y,color="r",linewidth=8)
    else:
        ax.plot(xA,yA,color="r",linewidth=4)
        ax.plot(xC,yC,color="r",linewidth=4)

        
    distBC_connnected = find_distance_global(B0Glob, C1Glob)+find_distance_global(B1Glob, C0Glob)
    distBC_disconnected =SB+SC 
    if distBC_connnected < distBC_disconnected:
        x,y = get_geodesic_global_to_cartesian(B0Glob, C1Glob)
        ax.plot(x,y,color="y",linewidth=8,label=r'$\mathcal{E}(BC)$')
        x,y = get_geodesic_global_to_cartesian(B1Glob, C0Glob)
        ax.plot(x,y,color="y",linewidth=8)
    else:
        ax.plot(xB,yB,color="y",linewidth=2)
        ax.plot(xC,yC,color="y",linewidth=2)
    """
        
    

    A0Cart = polar_to_cartesian(pointA0)
    A1Cart = polar_to_cartesian(pointA1)
    
    
    plot_entanglement_wedge_from_polars(pointA0, pointA1, pointB0, pointB1, ax, r'$E_W(A:CD)$', "g")
    plot_entanglement_wedge_from_polars(pointA0, pointA1, pointC0, pointC1, ax, r'$E_W(A:BD)$', "y")
    
    plot_entanglement_wedge_from_polars_E_AD(pointA0, pointA1, pointB0, pointB1,pointC0,pointC1,ax,r'$E_W(A:D)$',"b")
    
    
    
    EW_A_CD = get_EW_length_from_polars(pointA0, pointA1, pointB0, pointB1)
    EW_A_BD = get_EW_length_from_polars(pointA0, pointA1, pointC0, pointC1)
    EW_A_D = find_entanglement_wedge_dist_E_AD_polars(pointA0,pointA1,pointB0,pointB1,pointC0,pointC1)
    
    EQ = SB - S_AB + EW_A_CD + SC-S_AC+EW_A_BD-S_BC-EW_A_D + S_ABC
    print(EQ)
    
    
    """
    
    minDist0 = np.inf
    minDist1 = np.inf
    minPoint0 = None
    minPoint1 = None
  
    for i in range(len(x)):
        point = cartesian_to_global(x[i], y[i])
        dist0 = find_distance_global(point, A0Glob)
        dist1 = find_distance_global(point, A1Glob)
        if dist0 < minDist0:
            minDist0 = dist0
            minPoint0 = i
        if dist1 < minDist1:
            minDist1 = dist1
            minPoint1 = i
            
            
    #print(convert_angle_to_latex_angle(cartesian_to_polar([x[minPoint0],y[minPoint0]])))
   # print(convert_angle_to_latex_angle(cartesian_to_polar([x[minPoint1],y[minPoint1]])))
   
    xb,yb = get_geodesic_cartesian_to_cartesian(A0Cart, [x[minPoint0],y[minPoint0]])
    ax.plot(xb,yb)
    
    xb,yb = get_geodesic_cartesian_to_cartesian(A1Cart, [x[minPoint1],y[minPoint1]])
    ax.plot(xb,yb)
    
    x,y = get_geodesic_polar_to_cartesian(pointC0, pointC1)
    ax.plot(x,y,label="GEOC",color='k')
    
    minDist0 = np.inf
    minDist1 = np.inf
    minPoint0 = None
    minPoint1 = None
  
    for i in range(len(x)):
        point = cartesian_to_global(x[i], y[i])
        dist0 = find_distance_global(point, A0Glob)
        dist1 = find_distance_global(point, A1Glob)
        if dist0 < minDist0:
            minDist0 = dist0
            minPoint0 = i
        if dist1 < minDist1:
            minDist1 = dist1
            minPoint1 = i
            
            
    print(convert_angle_to_latex_angle(cartesian_to_polar([x[minPoint0],y[minPoint0]])))
    print(convert_angle_to_latex_angle(cartesian_to_polar([x[minPoint1],y[minPoint1]])))
   
    xb,yb = get_geodesic_cartesian_to_cartesian(A0Cart, [x[minPoint0],y[minPoint0]])
    ax.plot(xb,yb)
    
    xb,yb = get_geodesic_cartesian_to_cartesian(A1Cart, [x[minPoint1],y[minPoint1]])
    ax.plot(xb,yb)
    
    
    xb,yb = get_geodesic_polar_to_cartesian(pointB0, pointB1)    
    xc,yc = get_geodesic_polar_to_cartesian(pointC0, pointC1)

   
    minPoint0 = None
    minDist0 = np.inf
     
    minPoint1 = None
    
    for i in range(len(xb)):
         point0 = cartesian_to_global(xb[i], yb[i])
         for j in range(len(xc)):
             point1 = cartesian_to_global(xc[i], yc[i])
             dist = find_distance_global(point0, point1)
             if dist < minDist0:
                 minPoint0 = point0
                 minPoint1 = point1
                 minDist0 = dist
                 
    

    xb,yb = get_geodesic_global_to_cartesian(minPoint0, minPoint1)
    ax.plot(xb,yb)
    

   """

      
    #ax.set_title('Geodesics in the Poincaré Disk')
    ax.legend(loc='upper right') 
    plt.axis('off')

    plt.savefig("entanglementcrosssections.pdf", format="pdf")
    plt.show()


    

A0 = np.pi/6
A1 = 5*np.pi/6
B0 = np.pi
B1 = 3*np.pi/2 * 0.95
C0 = 3*np.pi/2 * 1.05
C1 = 2*np.pi


params = [A0,A1,B0,B1,C0,C1]

generate_points_for_graph_latex(params)