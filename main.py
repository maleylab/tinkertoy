from molecule import *
from helpers import *
import argparse
import random


def prepare_core(core):

    core_atoms, core_coordinates, core_connectivity = read_structure(core)
    connection_points, functionalization_points = parse_connectivity(core_connectivity)
    core = Core()
    core.atoms = core_atoms
    core.coordinates = core_coordinates
    core.connectivity = connection_points
    core.functionalization = functionalization_points

    return core

    
def prepare_fragments(fragments, num_of_core_connections):

    fragments_for_building = []
    
    for _ in range(num_of_core_connections):
        fragment_choice = random.choice(fragments)
        frag_atoms, frag_coordinates, frag_connectivity = read_structure(fragment_choice)
        frag_cps, frag_fps = parse_connectivity(frag_connectivity)
        fragment = Fragment()
        fragment.instance = _+1
        fragment.atoms = frag_atoms
        fragment.coordinates = frag_coordinates
        fragment.connectivity = frag_cps
        fragment.read_functionalization(frag_fps)
        fragments_for_building.append(fragment)
        
    return fragments_for_building


def prepare_functionalgroups(functionalgroups):
    
    r_groups = []
    
    for count, r_group in enumerate(functionalgroups):
        rgroup_atoms, rgroup_coordinates, rgroup_connectivity = read_structure(r_group)
        rgroup_cps, rgroup_fps = parse_connectivity(rgroup_connectivity)
        functionalgroup = FunctionalGroup()
        functionalgroup.atoms = rgroup_atoms
        functionalgroup.coordinates = rgroup_coordinates
        functionalgroup.connectivity = rgroup_cps
        connection = functionalgroup.connectivity[0]
        translate, rotate = int(connection[0])-1, int(connection[1])-1
        functionalgroup.translate_coordinates(translate)
        functionalgroup.rotate_coordinates(rotate)
        r_groups.append(functionalgroup)
    
    return r_groups


def rotation(fragment, rotate):
       
    coordinates = fragment.coordinates
    vector = coordinates[rotate]
    x_unit = np.array([-1,0,0])
  
    # Normalize vector
    vector /= np.linalg.norm(vector)
    
    # Get axis
    uvw = np.cross(vector, x_unit)
    
    rcos = np.dot(vector, x_unit)
    rsin = np.linalg.norm(uvw)
    
    if not np.isclose(rsin, 0):
        uvw /= rsin
    u,v,w = uvw
    
    rotation_matrix = ( rcos * np.eye(3) + 
                       rsin * np.array([
                            [0, -w, v],
                            [w, 0, -u],
                            [-v, u, 0]   
                       ]) + 
                       (1.0 - rcos) * uvw[:,None] * uvw[None,:])
    
    rotated_coordinates = []
    for coordinate in coordinates:
        new_coordinate = np.dot(coordinate.T, rotation_matrix.T)
        rotated_coordinates.append(new_coordinate)
    
    fragment.coordinates = np.array(rotated_coordinates)

def rotation2(fragment, rotate):
       
    coordinates = fragment.coordinates
    vector = coordinates[rotate]
    z_unit = np.array([1,0,0])
  
    # Normalize vector
    vector /= np.linalg.norm(vector)
    
    # Get axis
    uvw = np.cross(vector, z_unit)
    
    rcos = np.dot(vector, z_unit)
    rsin = np.linalg.norm(uvw)
    
    if not np.isclose(rsin, 0):
        uvw /= rsin
    u,v,w = uvw
    
    rotation_matrix = ( rcos * np.eye(3) + 
                       rsin * np.array([
                            [0, -w, v],
                            [w, 0, -u],
                            [-v, u, 0]   
                       ]) + 
                       (1.0 - rcos) * uvw[:,None] * uvw[None,:])
    
    rotated_coordinates = []
    for coordinate in coordinates:
        new_coordinate = np.dot(coordinate.T, rotation_matrix.T)
        rotated_coordinates.append([-new_coordinate[0]+0.5, new_coordinate[2], new_coordinate[1]])
    
    fragment.coordinates = np.array(rotated_coordinates)

def rotation3(fragment, rotate):
       
    coordinates = fragment.coordinates
    vector = coordinates[rotate]
    z_unit = np.array([1,0,0])
  
    # Normalize vector
    vector /= np.linalg.norm(vector)
    
    # Get axis
    uvw = np.cross(vector, z_unit)
    
    rcos = np.dot(vector, z_unit)
    rsin = np.linalg.norm(uvw)
    
    if not np.isclose(rsin, 0):
        uvw /= rsin
    u,v,w = uvw
    
    rotation_matrix = ( rcos * np.eye(3) + 
                       rsin * np.array([
                            [0, -w, v],
                            [w, 0, -u],
                            [-v, u, 0]   
                       ]) + 
                       (1.0 - rcos) * uvw[:,None] * uvw[None,:])
    
    rotated_coordinates = []
    for coordinate in coordinates:
        new_coordinate = np.dot(coordinate.T, rotation_matrix.T)
        rotated_coordinates.append([new_coordinate[0], new_coordinate[1], new_coordinate[2]])
    
    fragment.coordinates = np.array(rotated_coordinates)

        
def main(core, fragments, functionalgroups):
    
    # Functional groups are being prepared correctly.
    
    r_groups = prepare_functionalgroups(functionalgroups)
    
    # Core is being prepare correctly.
    Core = prepare_core(core)
    num_core_connections = sum(1 for _ in Core.connectivity)
    
    # Fragments are being prepare correctly.
    fragments_for_building = prepare_fragments(fragments, num_core_connections)
    
    for frag_count, build_frag in enumerate(fragments_for_building):
        fcnl_points = build_frag.functionalization
        times_to_functionalize = 3
        #times_to_functionalize = random.choice(range(len(fcnl_points)+1))
        if times_to_functionalize == 0:
            pass
        else:
            for func_count,_ in enumerate(range(times_to_functionalize)):
                fcnl_point = random.choice(fcnl_points)
                translate, rotate = int(fcnl_point[0])-1, int(fcnl_point[1])-1
                index = fcnl_points.index(fcnl_point)
                build_frag.translate_coordinates(translate)
                rotation(build_frag, rotate)
                new_atoms = build_frag.atoms
                
                rgroup_graft = random.choice(r_groups)
                
                xe_stop = 0
                
                for count, atom in enumerate(rgroup_graft.atoms):
                    if atom == 'Xe':
                        xe_stop = count
                    else:
                        new_atoms.append(atom)

                new_coords = []
                
                for count, coordinate in enumerate(build_frag.coordinates):
                    if count == translate:
                        pass
                    else:
                        new_coords.append([coordinate[0], coordinate[1], coordinate[2]])
                
                for count, coordinate in enumerate(rgroup_graft.coordinates):
                    if count == xe_stop:
                        pass
                    else:
                        new_coords.append([coordinate[0], coordinate[1], coordinate[2]])
                
                new_atoms.pop(translate)
                build_frag.coordinates = np.array(new_coords)
                fcnl_points.pop(index)
                
                new_fps = []
                fragment_fps = build_frag.functionalization
                
                for fp in fcnl_points:
                    dummy_atom, bound_atom = int(fp[0]), int(fp[1])
                    if (bound_atom > translate):
                        bound_atom += -1
                    elif (dummy_atom > translate):
                        dummy_atom += -1
                    else:
                        pass
                    new_fps.append([dummy_atom, bound_atom])
                
                fcnl_points = new_fps
        
        fragment_cps = build_frag.connectivity
        dummy_atom, bound_atom = fragment_cps[0][0], fragment_cps[0][1]
        bound_atom = int(bound_atom) - times_to_functionalize
        build_frag.connectivity = f'{dummy_atom}:{bound_atom}'
        
    
    for count, frag in enumerate(fragments_for_building):
        connection = frag.connectivity
        rotate, translate = int(connection[0][0])-1, int(connection[0][1])-1
        frag.translate_coordinates(translate)
        rotation2(frag,rotate)

    core_connections = Core.connectivity
    for count, _ in enumerate(core_connections):
        
        core_point = core_connections[count]
        translate, rotate = int(core_point[0])-1, int(core_point[1])-1
        Core.translate_coordinates(translate)
        rotation3(Core,rotate)

        core_coords = Core.coordinates.tolist()
        core_atoms = Core.atoms
        # Core is ready here
        frag = fragments_for_building[count]
        frag_atoms = frag.atoms
        frag_coordinates = frag.coordinates
        for _ in frag_coordinates:
            core_coords.append([_[0], _[1], _[2]])
        Core.coordinates = np.array(core_coords)
        for _ in frag_atoms:
            core_atoms.append(_)
        Core.atoms = core_atoms

        
    Core.write_xyz('Prepared_Structure')
                
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, required=True)
    parser.add_argument('-f', nargs='+', required=True)
    parser.add_argument('-r', nargs='+', default = [], required=False)
    args = parser.parse_args()
    core = args.c
    fragments = args.f
    functionalgroups = args.r
    main(core, fragments, functionalgroups)