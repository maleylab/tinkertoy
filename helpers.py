import numpy as np

def unit_vector(vector):
    try:
        return vector/np.linalg.norm(vector)
    except:
        return 1.0
        
def angle_between(v1,v2):
    
    try:
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        angle = np.arccos(np.clip(np.dot(v1_u,v2_u), -1.0, 1.0))
        return np.rad2deg(angle)
    except:
        return 0.0
    
def parse_connectivity(connectivity):
    
    connect_points, func_points = connectivity.split(';')[0], connectivity.split(';')[1]

    return connect_points, func_points

def read_structure(xyzfile):

    atom_list = []
    coordinates_list = []
    
    
    with open(xyzfile, 'r') as f:
        lines = [i.strip() for i in f.readlines()]
    
    connectivity = lines[1]
    
    for line in lines[2:]:
        line = line.split()
        atom_list.append(line[0])
        coordinates = [line[1], line[2], line[3]]
        coordinates = [float(i) for i in coordinates]
        coordinates_list.append(coordinates)

    return atom_list, np.array(coordinates_list), connectivity
