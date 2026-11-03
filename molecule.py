from helpers import unit_vector
from helpers import angle_between
import numpy as np
from scipy.spatial.transform import Rotation as R

class Molecule:
    def __init__(self):
        pass
    
    @property
    def atoms(self):
        return self.__atoms
        
    @property
    def coordinates(self):
        return self.__coordinates
    
    @atoms.setter
    def atoms(self, atoms:list):
        self.__atoms = atoms
    
    @coordinates.setter
    def coordinates(self, coordinates:list):
        self.__coordinates = coordinates
           
    def translate_coordinates(self, attachment_point:int):
        translation_matrix = -(self.__coordinates[attachment_point])
        translated_coordinates = translation_matrix + self.__coordinates
        self.__coordinates = translated_coordinates 
    
    def write_xyz(self, name:str):
   
        final_coordinates = []
        
        for i in (list(zip(self.__atoms, self.__coordinates))):
            if i[0] == 'Xe':
                pass
        
            elif i[0] == 'He':
                final_coordinates.append(f'H    {i[1][0]}    {i[1][1]}    {i[1][2]}')
            else:
                final_coordinates.append(f'{i[0]}    {i[1][0]}    {i[1][1]}    {i[1][2]}')
        
        with open(f'{name}.xyz', 'w') as f:
            f.write(f'{len(final_coordinates)}' +'\n')
            f.write(f'{name}' + '\n')
            for i in final_coordinates:
                f.write(i + "\n")


class Core(Molecule):
    
    @property
    def connectivity(self):
        return self.__connectivity

    @connectivity.setter
    def connectivity(self, cps:str):
        connection_list = []
        connections = cps.split(',')
        for connect in connections:
            connect = connect.split(':')
            dummy_atom, bound_atom = connect[0], connect[1]
            connection_list.append([dummy_atom, bound_atom])

        self.__connectivity = connection_list  
        
        
    @property
    def functionalization(self):
        return iter(self.__functionalization)   

    @functionalization.setter
    def functionalization(self, functionalization):
        functionalization_list = []
        functionalizations = functionalization.split(',')
        for fcnl_point in functionalizations:
            functionalization_list.append(fcnl_point)
        self.__functionalization = functionalization_list    
        
    def rotate_coordinates(self, bound_atom):
        vector = self.coordinates[bound_atom]
        
        rotation_angle = angle_between(vector, np.array([1,0,0]))
        r = R.from_rotvec(rotation_angle * np.array([0,0,1]), degrees=True)            

        rotated_coordinates = []
        
        for coordinate in self.coordinates:
            x,y,z = coordinate[0], coordinate[1], coordinate[2]
            vector = np.array([x,y,z])
            new_vector = np.dot(r.as_matrix(), vector)
            rotated_coordinates.append([new_vector[0], new_vector[1], new_vector[2]])
    
        self.coordinates = np.array(rotated_coordinates) 
        

class Fragment(Molecule):
    
    @property
    def instance(self):
        return self.__instance


    @property
    def functionalization(self):
        return self.__functionalization

    
    @property
    def connectivity(self):
        return self.__connectivity

    @connectivity.setter
    def connectivity(self, cps):
        connection_list = []
        connections = (cps.split(':'))
        dummy_atom, bound_atom = connections[0], connections[1]
        connection_list.append([dummy_atom, bound_atom])
        self.__connectivity = connection_list
        

    @instance.setter
    def instance(self, instance):
        self.__instance = instance
    
    
    @functionalization.setter
    def functionalization(self, functionalization):
        self.__functionalization = functionalization


    def read_functionalization(self, functionalization):
        functionalization_list = []
        functionalizations = functionalization.split(',')
        for fcnl_point in functionalizations:
            dummy_atom, bound_atom = fcnl_point.split(':')[1], fcnl_point.split(':')[0]
            functionalization_list.append([dummy_atom, bound_atom])
        self.__functionalization = functionalization_list

    def rotate_coordinates(self, bound_atom):
        vector = self.coordinates[bound_atom]
        
        rotation_angle = angle_between(vector, np.array([-1,0,0]))
        r = R.from_rotvec(rotation_angle * np.array([0,0,1]), degrees=True)            

        rotated_coordinates = []
        
        for coordinate in self.coordinates:
            x,y,z = coordinate[0], coordinate[1], coordinate[2]
            vector = np.array([x,y,z])
            new_vector = np.dot(r.as_matrix(), vector)
            rotated_coordinates.append([new_vector[0], new_vector[1], new_vector[2]])
    
        self.coordinates = np.array(rotated_coordinates) 


class FunctionalGroup(Molecule):

    @property
    def connectivity(self):
        return self.__connectivity

    @connectivity.setter
    def connectivity(self, cps:str):
        connection_list = []
        connections = cps.split(':')
        dummy_atom, bound_atom = connections[0], connections[1]
        connection_list.append([dummy_atom, bound_atom])

        self.__connectivity = connection_list    
        
    def rotate_coordinates(self, bound_atom):
        vector = self.coordinates[bound_atom]
        
        rotation_angle = angle_between(vector, np.array([-1,0,0]))
        r = R.from_rotvec(rotation_angle * np.array([0,0,1]), degrees=True)            

        rotated_coordinates = []
        
        for coordinate in self.coordinates:
            x,y,z = coordinate[0], coordinate[1], coordinate[2]
            vector = np.array([x,y,z])
            new_vector = np.dot(r.as_matrix(), vector)
            rotated_coordinates.append([-new_vector[0]-0.75, new_vector[1], new_vector[2]])
    
        self.coordinates = np.array(rotated_coordinates) 
        
        



    