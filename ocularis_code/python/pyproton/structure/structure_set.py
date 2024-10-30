#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from typing import Union
from collections.abc import Iterable
from scipy.spatial.transform import Rotation
from scipy import ndimage
import numpy as np
import copy

from pyproton.structure.structure import Structure
from pyproton.volume import Grid


class StructureSet:
    def __init__(self, structures: Union[None, Iterable[Structure]] = None):
        self._grid = None
        self._structures = {}

        if structures is not None:
            for s in structures:
                self.add_structure(s)

    def add_structure(self, struct: Structure):
        if self._grid is None:
            self._grid = struct.grid

        if struct.grid != self._grid:
            raise RuntimeError(
                "Trying to add structures on different grids to the "
                "same Structure Set"
            )

        self._structures[struct.organ.name] = struct

    @property
    def names(self) -> Iterable[str]:
        for s in self._structures.values():
            yield s.organ.name

    @property
    def name_list(self) -> list:
        name_list = []
        for s in self._structures.values():
            name_list.append(s.organ.name)
        return name_list

    @property
    def grid(self) -> Union[None, Grid]:
        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid
        for structure_name in self.name_list:
            self._structures[structure_name].grid = self._grid

    def translate_contour_coordinates(self, vector):
        for structure_name in self.name_list:
            for idx, _ in enumerate(self._structures[structure_name].contour_coordinates):
                self._structures[structure_name].contour_coordinates[idx] -= vector

    def centre_point_cloud_to_structure(self, structure_name):
        if structure_name not in self.name_list:
            print("Structure not part of structure set")
            return
        self.translate_contour_coordinates(
            self._structures[structure_name].com)

        print(f"Centered structure point clouds to {structure_name}")

    def centre_grid_to_structure(self, structure_name):
        if structure_name not in self.name_list:
            print("Structure not part of structure set")
            return
        com_of_voxels = np.asarray(ndimage.measurements.center_of_mass(
            self._structures[structure_name].binary_mask))
        target_com = com_of_voxels*self.grid.spacing + self.grid.origin

        new_grid = copy.deepcopy(self._grid)

        new_grid.origin -= target_com
        new_grid.meshgrids.update(new_grid)
        self.grid = new_grid

        print(
            f"Centered structure set voxel grids to {structure_name}")

    def rotate(self, rot_vector):
        r = Rotation.from_euler(
            'xyz', (rot_vector[0], rot_vector[1], rot_vector[2]), degrees=True)

        for structure_name in self.name_list:
            self._structures[structure_name].contour_coordinates = r.apply(
                self._structures[structure_name].contour_coordinates)

        print(f"Rotated point cloud by {rot_vector}")

    def resample_contour_points_in_grid(self, structure_names: list = None):

        if False not in [name in structure_names for name in self.name_list]:
            print("One or more structures not part of structure set")
            return

        if structure_names == None:
            structure_names = self.name_list

        for structure_name in structure_names:
            self._structures[structure_name].resample_contour_points_in_grid()

    def __getitem__(self, name: str) -> Structure:
        return self._structures[name]

    def __iter__(self) -> Structure:
        for name in self.names:
            yield self[name]

    def __repr__(self):
        title = "Structure set consisting of "
        for name in self.names:
            title += name + ", "
        return title
