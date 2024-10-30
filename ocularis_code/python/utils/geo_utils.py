# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np


def points_outside_cylinder(points, axis_point, axis_direction, radius, infront_of, behind):
    # Normalize the axis direction vector
    axis_direction = axis_direction / np.linalg.norm(axis_direction)
    
    # Calculate vectors from axis_point to each point in the points array
    vectors = points - axis_point
    
    # Project the vectors onto the axis direction to find their distance along the axis
    distances_along_axis = np.dot(vectors, axis_direction)
    
    # Project vectors onto the plane orthogonal to the axis to find their distance from the axis
    orthogonal_vectors = vectors - np.outer(distances_along_axis, axis_direction)
    distances_from_axis = np.linalg.norm(orthogonal_vectors, axis=1)
    
    # Check if points are inside the cylinder based on their distances
    inside_cylinder = (distances_along_axis >= infront_of) & (distances_along_axis <= behind) & (distances_from_axis <= radius)
    
    # Return points that are outside the cylinder
    return points[~inside_cylinder]

def find_maximum_radius(points, axis_point, axis_direction):
    """
    Finds the maximum radius of points in a numpy array from an axis defined by an axis point and direction.

    Parameters:
    points (numpy.ndarray): Nx3 array of points.
    axis_point (numpy.ndarray): The point on the axis (1x3 array).
    axis_direction (numpy.ndarray): The direction vector of the axis (1x3 array).

    Returns:
    float: The maximum radius of the points from the axis.
    """

    # Normalize the axis direction
    axis_direction = axis_direction / np.linalg.norm(axis_direction)

    # Find the vectors from the axis point to each point
    vectors_from_axis = points - axis_point

    # Project these vectors onto the plane perpendicular to the axis direction
    # by subtracting their projection on the axis direction from themselves
    vectors_on_plane = vectors_from_axis - np.outer(np.dot(vectors_from_axis, axis_direction), axis_direction)

    # Calculate the distance of each point to the axis (the length of the projection on the plane)
    radii = np.linalg.norm(vectors_on_plane, axis=1)

    # Return the maximum radius
    return np.max(radii)

def filter_points_within_radius(points, axis_point, axis_direction, radius):
    """
    Filters out points that are within a given radius from an axis defined by an axis point and direction.

    Parameters:
    points (numpy.ndarray): Nx3 array of points.
    axis_point (numpy.ndarray): The point on the axis (1x3 array).
    axis_direction (numpy.ndarray): The direction vector of the axis (1x3 array).
    radius (float): The radius within which points should be filtered out.

    Returns:
    numpy.ndarray: The filtered array of points.
    """

    # Normalize the axis direction
    axis_direction = axis_direction / np.linalg.norm(axis_direction)

    # Find the vectors from the axis point to each point
    vectors_from_axis = points - axis_point

    # Project these vectors onto the plane perpendicular to the axis direction
    vectors_on_plane = vectors_from_axis - np.outer(np.dot(vectors_from_axis, axis_direction), axis_direction)

    # Calculate the distance of each point to the axis (the length of the projection on the plane)
    radii = np.linalg.norm(vectors_on_plane, axis=1)

    # Filter out points where the radius is less than or equal to the given radius
    filtered_points = points[radii > radius]

    return filtered_points