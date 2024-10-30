#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from scipy import ndimage
import numpy as np

from pyproton.structure.structure_set import StructureSet


class EyeStructureSet(StructureSet):

    @property
    def gaze_vector(self):
        # Based on point cloud

        if "Sclera" in self.name_list:
            sclera_com = self["Sclera"].com

        elif "Sklera" in self.name_list:
            sclera_com = self["Sklera"].com
        else:
            return np.zeros(3)

        if "Lens" in self.name_list:
            lens_com = self["Lens"].com

        elif "Lens_L" in self.name_list:
            lens_com = self["Lens_L"].com
        else:
            return np.zeros(3)

        # sclera_com = self["Sclera"].com
        # lens_com = self["Lens"].com

        gaze_vector = lens_com - sclera_com
        gaze_vector = gaze_vector / np.linalg.norm(gaze_vector)

        return gaze_vector

    @property
    def gaze_vector_voxel(self):
        return np.zeros(3)
    #     # Based on voxel distribution
    #     struct = self["Sclera"]
    #     mask = struct.binary_mask
    #     com_bins = np.asarray(ndimage.measurements.center_of_mass(mask))
    #     sclera_com_voxel = com_bins*struct.grid.spacing + struct.grid.origin

    #     struct = self["Lens"]
    #     mask = struct.binary_mask
    #     com_bins = np.asarray(ndimage.measurements.center_of_mass(mask))
    #     lens_com_voxel = com_bins*struct.grid.spacing + struct.grid.origin

    #     gaze_vector_voxel = lens_com_voxel - sclera_com_voxel
    #     gaze_vector_voxel = gaze_vector_voxel / \
    #         np.linalg.norm(gaze_vector_voxel)

    #     return gaze_vector_voxel
