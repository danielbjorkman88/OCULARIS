#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


from pyproton.volume.grid import Grid
import numpy as np
from pyproton.utils.constants import PKL_FILE_EXTENSION
from pyproton.utils.ioUtils import is_pickle, write_pickle


class Volume:
    """
    General class to deal with volumes.
    Currently only supports axis aligned scans (no rotation  in voxel to world matrix)
    Arrays are stored as C,X,Y,Z, C = 1 for grayscale images, C=3 for RGB images, C = 3 for vectorfields
    """

    def __init__(self, array, voxel_to_world_matrix):
        """
        Parameters
        ----------
        array: np.array
            Numpy array of which you want to create a volume.
        voxel_to_world_matrix:
            Matrix containing the information to convert voxel coordinates to world coordinates.
        """
        self._pixel_array = array
        self._voxel_to_world_matrix = voxel_to_world_matrix
        self._set_channel_maybe()
        self._grid = None

    @property
    def pixel_array(self) -> np.ndarray:
        return self._pixel_array.squeeze()

    @property
    def voxel_to_world_matrix(self) -> np.ndarray:
        return self._voxel_to_world_matrix

    @property
    def shape(self) -> tuple:
        return tuple(self.pixel_array.shape[-3:])

    @property
    def voxel_spacing(self) -> np.ndarray:
        return np.linalg.norm(self._voxel_to_world_matrix[:3,:3],axis=0)

    @property
    def origin(self) -> np.ndarray:
        return self._voxel_to_world_matrix[:3, 3]

    @property
    def orientation(self) -> np.ndarray:
        return self._voxel_to_world_matrix[:3, :3] / self.voxel_spacing

    @property
    def grid(self):
        if self._grid is None:        
            self._grid = Grid(self.origin, self.voxel_spacing, np.asarray(self.pixel_array.shape), "xx")
        return self._grid

    def _set_pixel_array(self, new_array):
        self._pixel_array = new_array

    def _set_voxel_to_world_matrix(self, new_voxel_to_world_matrix):
        self._voxel_to_world_matrix = new_voxel_to_world_matrix

    def _set_origin(self, new_origin):
        self._voxel_to_world_matrix[:3, 3] = new_origin

    def _set_voxel_spacing(self, new_voxel_spacing):
        self._voxel_to_world_matrix[:3, :3] = self.orientation*new_voxel_spacing

    def _set_channel_maybe(self):
        if len(self._pixel_array.shape) == 3:
            self._pixel_array = np.expand_dims(self._pixel_array, 0)

    def save(self, save_path):
        assert is_pickle(
            save_path
        ), f"Volume objects must be saved with a '{PKL_FILE_EXTENSION}' extension"
        write_pickle(save_path, self)
