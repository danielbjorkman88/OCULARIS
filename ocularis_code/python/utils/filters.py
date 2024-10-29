# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 10:26:58 2024

@author: bjoerk_c
"""
from scipy.spatial import KDTree
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
import numpy as np
from scipy.ndimage import binary_erosion

from scipy.spatial import Delaunay



def make_spherical_object_hollow_thin(binary_array):
    """
    Convert a solid spherical object in a binary array to a hollow one with a 1-voxel thick shell.
    This method erodes the object and subtracts this erosion from the original to get a thin shell.

    Parameters:
    - binary_array (numpy.ndarray): A binary array with the spherical object represented by ones (1s).

    Returns:
    - numpy.ndarray: A binary array where only a 1-voxel thick outer shell of the spherical object remains.
    """
    # Perform erosion to remove one layer from the spherical object
    eroded_array = binary_erosion(binary_array, structure=np.ones((3, 3, 3)))

    # Calculate the difference between the original array and the eroded one
    hollow_array = binary_array ^ eroded_array

    return hollow_array.astype(np.bool)
def find_voxel_centers_vectorized(binary_mask, mesh_origin, resolution):
    """
    Find the centers of all voxels marked by 1s in a binary mask, using vectorized operations.

    Parameters:
    - binary_mask (numpy.ndarray): A 3D binary array representing the object.
    - mesh_origin (tuple or numpy.ndarray): The origin (x0, y0, z0) of the mesh grid.
    - resolution (tuple or numpy.ndarray): The resolution (dx, dy, dz) of the mesh grid.

    Returns:
    - numpy.ndarray: An array of voxel center points.
    """
    # Find indices of all active voxels
    x_indices, y_indices, z_indices = np.where(binary_mask == 1)
    
    # Calculate the centers of these voxels
    x_centers = x_indices * resolution[0] + mesh_origin[0] + resolution[0] / 2
    y_centers = y_indices * resolution[1] + mesh_origin[1] + resolution[1] / 2
    z_centers = z_indices * resolution[2] + mesh_origin[2] + resolution[2] / 2

    # Stack the coordinates to get an array of center points
    voxel_centers = np.vstack((x_centers, y_centers, z_centers)).T

    return voxel_centers


def segment_and_densify_point_cloud(points, num_segments, densify_resolution):
    """
    Segment a point cloud into regions defined by planes rotating around the z-axis,
    densify each segment individually, and then recombine them.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the point cloud.
    - num_segments (int): The number of segments to divide the point cloud into.
    - densify_function (function): A function that takes points and a resolution
                                   and returns a densified point cloud.
    - densify_resolution (float): The resolution parameter for the densification function.

    Returns:
    - numpy.ndarray: The recombined densified point cloud.
    """
    # Calculate the angle for each point in the XY plane relative to the positive X-axis
    angles = np.arctan2(points[:, 1], points[:, 0])
    
    # Normalize angles to the range [0, 2*pi)
    angles = np.mod(angles, 2*np.pi)
    
    # Determine the boundaries of the segments
    segment_boundaries = np.linspace(0, 2*np.pi, num_segments + 1)
    
    # Segment the points
    segmented_points = [None] * num_segments
    for i in range(num_segments):
        # Find points within the current segment
        in_segment = (angles >= segment_boundaries[i]) & (angles < segment_boundaries[i+1])
        segment_points = points[in_segment]
        
        # Densify the points within the current segment
        if len(segment_points) > 0:
            segmented_points[i] = densify_point_cloud(segment_points, densify_resolution)
        else:
            segmented_points[i] = np.array([]).reshape(0, 3)

    # Combine all densified segments
    densified_points = np.vstack(segmented_points)
    
    return densified_points


def filter_points_relative_to_plane(points, axis_point, axis_direction, distance_behind):
    """
    Filter points based on their position relative to a plane orthogonal to the axis direction,
    located a specified distance behind the axis point.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the point cloud.
    - axis_point (numpy.ndarray): The reference point (3D) on the axis.
    - axis_direction (numpy.ndarray): The direction vector (3D) of the axis.
    - distance_behind (float): The distance behind the axis point where the plane is located.

    Returns:
    - (numpy.ndarray, numpy.ndarray): A tuple of two arrays:
        - The first array contains points in front of the plane.
        - The second array contains points behind the plane.
    """
    # Normalize the axis direction vector
    axis_direction_normalized = axis_direction / np.linalg.norm(axis_direction)
    
    # Calculate the new axis point (threshold axis point) that is distance_behind behind the original axis point
    threshold_axis_point = axis_point - axis_direction_normalized * distance_behind
    
    # Calculate the vector from this new threshold axis point to each point
    vectors_to_points_from_threshold = points - threshold_axis_point
    
    # Project these vectors onto the axis direction to determine their position relative to the plane
    projections_from_threshold = np.dot(vectors_to_points_from_threshold, axis_direction_normalized)
    
    # Filter points based on their projection relative to the threshold axis point
    points_in_front_of_plane = points[projections_from_threshold > 0]  # Points in front of the plane
    points_behind_plane = points[projections_from_threshold <= 0]  # Points behind the plane

    return points_in_front_of_plane, points_behind_plane
    

def project_furthest_point_behind(points, axis_point, axis_direction):
    """
    Find the point in a set of points that is the furthest behind relative to an axis point and direction,
    then project this point onto the axis.

    Parameters:
    - points (numpy.ndarray): Nx3 array of points.
    - axis_point (numpy.ndarray): The reference point (3D) on the axis.
    - axis_direction (numpy.ndarray): The direction vector (3D) of the axis.

    Returns:
    - numpy.ndarray: The projected point (3D) on the axis of the point that is furthest behind.
    """
    # Normalize the axis direction vector
    axis_direction_normalized = axis_direction / np.linalg.norm(axis_direction)
    
    # Calculate vectors from the axis point to each point
    vectors_to_points = points - axis_point
    
    # Project vectors onto the axis direction
    projections = np.dot(vectors_to_points, axis_direction_normalized)
    
    # Identify the point furthest behind (minimum projection value)
    min_projection_index = np.argmin(projections)
    furthest_point_behind = points[min_projection_index]
    
    # Calculate the projection of this point onto the axis
    projection_distance = projections[min_projection_index]
    projected_point_on_axis = axis_point + projection_distance * axis_direction_normalized
    
    return projected_point_on_axis

def densify_point_cloud_3d(points, resolution=0.5):
    """
    Attempt to densify a 3D point cloud by adding new points within its volume.
    This approach uses 3D Delaunay triangulation to define the volume.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the original point cloud.
    - resolution (float): The desired density for adding points, in terms of spatial resolution.

    Returns:
    - numpy.ndarray: A densified version of the point cloud.
    """
    # Perform 3D Delaunay triangulation on the point cloud
    tri = Delaunay(points)
    
    # Generate potential points within the bounding box of the original point cloud
    x_min, y_min, z_min = np.min(points, axis=0) - resolution
    x_max, y_max, z_max = np.max(points, axis=0) + resolution
    x_range = np.arange(x_min, x_max, resolution)
    y_range = np.arange(y_min, y_max, resolution)
    z_range = np.arange(z_min, z_max, resolution)
    grid_x, grid_y, grid_z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
    potential_points = np.vstack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()]).T

    # Find which potential points are inside the convex hull defined by the triangulation
    inside_points = potential_points[tri.find_simplex(potential_points) >= 0]

    # Combine the original points with the new points inside the volume
    # Note: This might include a significant number of points; consider if resolution and volume are appropriate
    densified_points = np.vstack([points, inside_points])

    return densified_points

def densify_point_cloud(points, resolution = 0.5):
    """
    Densify a point cloud by interpolating new points.

    Parameters:
    - points (numpy.ndarray): Nx3 array representing the original point cloud.
    - resolution (float): The desired resolution for the densified point cloud.

    Returns:
    - numpy.ndarray: Mx3 array representing the densified point cloud.
    """
    
    assert len(points) > 2
    
    # Define the bounding box of the original point cloud
    x_min, y_min, z_min = np.min(points, axis=0)
    x_max, y_max, z_max = np.max(points, axis=0)

    # Generate a grid of points within the bounding box with the desired resolution
    x_grid, y_grid, z_grid = np.meshgrid(np.arange(x_min, x_max, resolution),
                                         np.arange(y_min, y_max, resolution),
                                         np.arange(z_min, z_max, resolution))

    # Interpolate the values of the new points based on the original points
    dense_points = griddata(points[:, :2], points[:, 2], (x_grid, y_grid), method='linear')

    # Filter out NaN values
    valid_indices = ~np.isnan(dense_points)
    dense_points = np.column_stack((x_grid.ravel()[valid_indices.ravel()],
                                    y_grid.ravel()[valid_indices.ravel()],
                                    dense_points.ravel()[valid_indices.ravel()]))


    densified_points = np.vstack([points, dense_points])

    return densified_points





def filter_for_neighbors(points_to_be_filtered, R, epsilon=1e-6):
    """
    Filter points that have at least 2 neighbors within distance R.

    Parameters:
        points_to_be_filtered (numpy.ndarray): Array of shape (n, 3) representing points in 3D space.
        R (float): Radius within which points are considered near each other.
        epsilon (float, optional): Small distance to exclude points that are too close. Defaults to 1e-6.

    Returns:
        numpy.ndarray: Subset of points that have at least 2 neighbors within distance R.
    """
    tree = KDTree(points_to_be_filtered)
    points_of_interest = np.zeros(len(points_to_be_filtered), dtype=bool)
    
    for i, query_point in enumerate(points_to_be_filtered):
        # Find the neighbors within distance R, excluding the point itself
        neighbors = tree.query_ball_point(query_point, r=R)
        
        # Exclude neighbors that are too close to the query point
        filtered_neighbors = [idx for idx in neighbors if np.linalg.norm(query_point - points_to_be_filtered[idx]) > epsilon]
        
        # Check if the point has at least 2 valid neighbors (including itself) within distance R
        if len(filtered_neighbors) >= 5:  # 2 because the point itself is included in neighbors
            points_of_interest[i] = True
    
    # Return the subset of points that have at least 2 neighbors within distance R
    return points_to_be_filtered[points_of_interest]




    
def sphere_func(coords, x0, y0, z0, r):
    """Defines a sphere function."""
    x, y, z = coords
    return (x - x0)**2 + (y - y0)**2 + (z - z0)**2 - r**2

def fit_sphere(points):
    """Fits a sphere to a set of points."""
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    # Initial guess for the sphere parameters (center and radius)
    initial_guess = (np.mean(x), np.mean(y), np.mean(z), np.std(x + y + z))
    # Fit the sphere parameters using curve_fit
    popt, pcov = curve_fit(sphere_func, (x, y, z), np.zeros_like(x), p0=initial_guess)
    # Extract the sphere parameters
    x0, y0, z0, r = popt
    return (x0, y0, z0), r

def filter_sphere(my_points, sphere_points):
    """
    Filter points in my_points that are outside the sphere defined by sphere_points.

    Parameters:
    - my_points: numpy array of shape (N, 3) representing points to be filtered
    - sphere_points: numpy array of shape (M, 3) representing points defining the sphere

    Returns:
    - outside_points_array: numpy array of shape (P, 3) containing points from my_points
                            that are outside the sphere
    """

    # Fit a sphere to the sphere_points
    sphere_center, sphere_radius = fit_sphere(sphere_points)

    # Compute the squared distance of each point from the sphere center
    squared_distances = np.sum((my_points - sphere_center)**2, axis=1)

    # Find points that are outside the sphere
    outside_points = squared_distances > sphere_radius**2

    # Get the points that are outside the sphere
    outside_points_array = my_points[outside_points]

    return outside_points_array
