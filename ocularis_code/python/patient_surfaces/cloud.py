# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:45:41 2022

@author: Daniel Björkman
"""



import pyvista as pv
from sklearn.neighbors import KDTree

import patient_surfaces
from utils.com import find_com

import math
import numpy as np



def define_plane(input_points_for_plane):
    
    assert len(input_points_for_plane) == 3
    
    point1, point2, point3 = input_points_for_plane
    
    # Calculate vectors from point2 to point1 and point3 to point1
    vector1 = np.array(point1) - np.array(point2)
    vector2 = np.array(point3) - np.array(point2)

    # Calculate the normal vector to the plane using cross product
    normal_vector = np.cross(vector1, vector2)
    
    # Check if the cross product is zero (indicating collinearity)
    if np.all(normal_vector == 0):
        # Choose a different pair of vectors for the calculation
        vector1 = np.array(point2) - np.array(point1)
        vector2 = np.array(point3) - np.array(point1)
        normal_vector = np.cross(vector1, vector2)    

    # Choose one of the input points as a point on the plane
    plane_point = np.array(point1)

    return normal_vector, plane_point
    
    

def shell_ray_first_intersect(shell, ray):
    points, ind = shell.ray_trace(ray.S, ray.T)

    if len(points) == 0:
        return

    shortest_distance = 1000000
    first_intersect = []
    for p in points:
        dist = math.dist(p, ray.S)
        if dist < shortest_distance:
            shortest_distance = dist
            first_intersect = p

    cloud_surface_interaction_point = np.asarray(first_intersect)
    
    return cloud_surface_interaction_point

def points_ray_intersect(input_points, ray):
    
    
    tree = KDTree(input_points)

    surface_com = find_com(input_points)
    
    start_point = ray.S
    end_point = ray.T
    
    ray_direction = ray.vec_norm
    relative_position = surface_com - start_point
    projection_distance = np.dot(relative_position, ray_direction)
    projection_point = start_point + np.clip(projection_distance, 0, np.linalg.norm(end_point - start_point)) * ray_direction
    
    
    neighbor_candidates = []
    for offset in np.linspace(-3,3,6):
        projection_distance = np.dot(relative_position, ray_direction) + offset
        projection_point = start_point + np.clip(projection_distance, 0, np.linalg.norm(end_point - start_point)) * ray_direction
        
        k_nearest_neighbors = tree.query(projection_point.reshape(1, -1), k=3)
        distances, indices = k_nearest_neighbors
            
        neighbor_points = input_points[indices][0]
        
        
        neighbor_candidates.append( np.concatenate((neighbor_points, distances[0].reshape(-1, 1)), axis=1))
    
        
    concatenated_matrix = np.concatenate(neighbor_candidates, axis=0)
    
    first_three_columns = concatenated_matrix[:, :3]
    
    # Find unique rows based on the first three columns
    unique_rows, unique_indices = np.unique(first_three_columns, axis=0, return_index=True)
    
    concatenated_matrix = concatenated_matrix[unique_indices]
    
    
    
    indices_of_smallest_values = np.argsort(concatenated_matrix[:, 3])[:3]
    
    concatenated_matrix = concatenated_matrix[indices_of_smallest_values]
    
    # if np.min(concatenated_matrix[0:,3]) > 1:
    #     #ray_intersect = ray.find_intersect_with_plane(ray.vec_norm, concatenated_matrix[0,0:3])
        
    #     ray_intersect = ray.find_intersect_with_plane(ray.vec_norm, np.asarray([0,0,skin_plane]))
           
        
    #     ray.surface_interaction_point = ray_intersect
    #     continue
    
    
    points_closest_to_axis = concatenated_matrix[0:,0:3]
    
    normal_vector, plane_point = define_plane(points_closest_to_axis)
    
    ray_intersect = ray.find_intersect_with_plane(normal_vector, plane_point)
    
    # ray.surface_interaction_point = ray_intersect

    cloud_surface_interaction_point = np.asarray(ray_intersect)
    
    return cloud_surface_interaction_point


class PointCloud(patient_surfaces.PatientSurface):
    def __init__(self, points, alpha=12. , irregular_cloud = False):
        self.type = "cloud"
        self.points = points

        self.cloud = pv.PolyData(self.points)
        self.volume = self.cloud.delaunay_3d(alpha=alpha)
        # volume.plot()
        self.shell = self.volume.extract_geometry()
        
        self.irregular_cloud = irregular_cloud
        
    def __repr__(self):
        return "Point cloud"

    def find_ray_interaction(self, ray):

        # points, ind = self.shell.ray_trace(ray.S, ray.T)

        # if len(points) == 0:
        #     return

        # shortest_distance = 1000000
        # first_intersect = []
        # for p in points:
        #     dist = math.dist(p, ray.S)
        #     if dist < shortest_distance:
        #         shortest_distance = dist
        #         first_intersect = p

        # cloud_surface_interaction_point = np.asarray(first_intersect)
        
        if not self.irregular_cloud:
            cloud_surface_interaction_point = shell_ray_first_intersect(self.shell, ray)
        else:
            cloud_surface_interaction_point = points_ray_intersect(self.points, ray)
        # if (cloud_surface_interaction_point == cloud_surface_interaction_point2).all():
        #     pass
        # else:
        #     print("--------------------------------------------")
        #     print(cloud_surface_interaction_point, cloud_surface_interaction_point2)
        #     print("--------------------------------------------")
        #     raise ValueError
        
        ray.surface_interaction_point = cloud_surface_interaction_point
        ray.surface_interaction_point_candidates.append(
            cloud_surface_interaction_point)
        
    def get_ray_interaction(self, ray):
        cloud_surface_interaction_point = shell_ray_first_intersect(self.shell, ray)
        
        return cloud_surface_interaction_point

    def plot_volume(self):
        self.volume.plot()

    def plot_shell(self):
        self.shell.plot()
