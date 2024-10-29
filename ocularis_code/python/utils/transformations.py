# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:39:02 2022

@author: bjoerk_c
"""

import numpy as np
from scipy.spatial.transform import Rotation

def mapping(pointsIn, TMatrix):
    """ 
    function for mapping of 3d points by roto-translation Matrix 
    pointsIn = Nx3 points in three dimensional coordinates
    TMatrix = 4x4 rt matrix
    """
    if len(pointsIn.shape) == 1:
        pointsIn = np.expand_dims(pointsIn, axis=0)

    n, _ = pointsIn.shape
        
    pointsIn = np.c_[pointsIn, np.ones(n)]

    out = np.zeros(pointsIn.shape)

    for i in range(n):
        out[i, 0:] = TMatrix.dot(pointsIn[i].T)

    return out[0:, 0:3]



def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix



def rotation_matrix_from_vectors_scipy(vec1, vec2):
    """
    Calculate the rotation matrix that rotates vec1 to vec2.
    
    Parameters:
    - vec1: The initial vector.
    - vec2: The target vector.
    
    Returns:
    - rot_matrix: The rotation matrix.
    """
    # Normalize vectors
    vec1 = np.asarray(vec1) / np.linalg.norm(vec1)
    vec2 = np.asarray(vec2) / np.linalg.norm(vec2)
    
    # Calculate the rotation matrix
    rotation = Rotation.from_rotvec(np.cross(vec1, vec2))
    rot_matrix = rotation.as_matrix()
    
    return rot_matrix