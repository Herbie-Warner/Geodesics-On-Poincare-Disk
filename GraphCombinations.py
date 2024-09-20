"""
Given N points defined on S^1 code generates all the non crossing connections
between the vertices. Used to consider all possible geometries of the
entanglement wedge cross sections. 

Herbie Warner 10/07/2024
"""
import matplotlib.pyplot as plt
import numpy as np


def is_crossing(pair1, pair2):
    a, b = pair1
    c, d = pair2
    return (a < c < b < d) or (c < a < d < b)

def generate_non_crossing_matchings(points):
    if len(points) == 2:
        return [[(points[0], points[1])]]
    
    matchings = []
    for i in range(1, len(points), 2):
        first_pair = (points[0], points[i])
        remaining_points = points[1:i] + points[i+1:]      
        for sub_matching in generate_non_crossing_matchings(remaining_points):
            valid = True
            for pair in sub_matching:
                if is_crossing(first_pair, pair):
                    valid = False
                    break
            if valid:
                matchings.append([first_pair] + sub_matching)
    
    return matchings

def generate_valid_permutations(n):
    if n < 1:
        return []
    points = [i + 1 for i in range(2 * n)]
    matchings = generate_non_crossing_matchings(points)  
    return matchings

def plot_permutation(permutation, n):
    angles = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    points = [(np.cos(angle), np.sin(angle)) for angle in angles]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
   
    circle = plt.Circle((0, 0), 1, color='black', fill=False)
    ax.add_artist(circle)

    for i, (x, y) in enumerate(points):
        ax.plot(x, y, 'o', color='blue')
        ax.text(x * 1.1, y * 1.1, str(i + 1), horizontalalignment='center', 
                verticalalignment='center')

    for (start, end) in permutation:
        x_values = [points[start - 1][0], points[end - 1][0]]
        y_values = [points[start - 1][1], points[end - 1][1]]
        ax.plot(x_values, y_values, 'r-')

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    plt.show()
 
    
def plot_all_permutations(N):
    #To plot all permutations
    perm = generate_valid_permutations(N)
    print("Number of permutations: {}".format(len(perm)))
    for val in perm:
        plot_permutation(val, N)