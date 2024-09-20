# -*- coding: utf-8 -*-
"""
Created on Fri Aug 5 20:00:10 2024

@author: herbi
NON-OVERLAPPING DEFINITIONS!!!

This is formalised version of correctextremal surfaces but need to find time
to implement properly including all the combinatorics.
"""

import numpy as np
import matplotlib.pyplot as plt
from Utilities import epsilon,R,R_EFF,S_NUMBER
from Utilities import combinatoric_matching, find_distance_polars
from Utilities import get_geodesic_polar_to_cartesian
from GraphCombinations import generate_valid_permutations




def plot_entanglement_wedge(systems,ax,colour,name):
    """
    Systems assumed to be array of [[a,b],[c,d]...] where each [a,b] define the
    angles that define one subsystem, and [c,d] another etc. Assumes non-
    overlapping definitions.
    """
    number = len(systems)
    perms = generate_valid_permutations(number)

    minDist = np.inf
    bestSys = None
    for perm in perms:
        dist = 0
        angles = combinatoric_matching(perm, systems)
        for system in angles:
            dist += find_distance_polars(system[0], system[1])
            
        if dist < minDist:
            minDist = dist
            bestSys = angles
            
       
    x,y = get_geodesic_polar_to_cartesian([R_EFF,bestSys[0][0]], [R_EFF,bestSys[0][1]])
    ax.plot(x,y,color=colour,label=name)
    
    for system in bestSys[1:]:
        x,y = get_geodesic_polar_to_cartesian([R_EFF,system[0]], [R_EFF,system[1]])
        ax.plot(x,y,color=colour)
    
            
def plot_systems(systems,ax,names=None):
    if names != None and len(names) == len(systems):
        for system,name in zip(systems,names):
            angles = np.linspace(system[0], system[1], 100)
            ax.plot(np.cos(angles),np.sin(angles),label=name, linewidth=5)
    else:
        for system in systems:
            angles = np.linspace(system[0], system[1], 100)
            ax.plot(np.cos(angles),np.sin(angles), linewidth=5)



def main():  
    systems = [[0,np.pi/4],[np.pi/4*1.1,np.pi],[3*np.pi/2,1.95*np.pi]]
        
    fig, ax = plt.subplots(figsize=(8, 8))  
    theta = np.linspace(0, 2 * np.pi, 200)
    
    x_boundary = np.cos(theta)
    y_boundary = np.sin(theta)
    
    ax.plot(x_boundary, y_boundary, color='k')  
    ax.set_aspect('equal')  
    plot_systems(systems,ax,["A","B","C"])
        
        
    
    plot_entanglement_wedge(systems,ax,'r',r'$\mathcal{E}(ABC)$')
    
    
    ax.legend()
    #plt.savefig("Entanglementwedge.pdf",format='pdf')
    plt.show()
    
main()