# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.
"""

from typing import Union
from collections.abc import Iterable
import copy
from scipy import ndimage
import numpy as np

import optis_tps


class BeamSet:
    def __init__(self, tps_config, patient_model):
        self._beams = []
        self.tps_config = tps_config
        self._patient_model = patient_model

        self.title = "Title string not defined yet"

    @property
    def beams(self) -> list:
        if not self._beams:
            for theta, phi in zip(self.tps_config["beam_angles_theta"], self.tps_config["beam_angles_phi"]):
                beam_config = copy.deepcopy(self.tps_config["config"])

                beam_config["data"] = self.tps_config["data"]

                # rotation in xz plane [Degrees]
                beam_config["Gantry rot theta"] = theta
                # rotation in yz plane [Degrees]
                beam_config["Gantry rot phi"] = phi

                beam_config["grid"] = self._patient_model.structure_set.grid
                beam_config["target"] = self._patient_model.structure_set["GTV"].binary_mask
                beam_config["Sclera"] = self._patient_model.structure_set["Sclera"].binary_mask
                beam_config["Medium"] = "Dicom"
                beam_config["Image"] = "3D" # 3D, xslice, yslice, zslice
                #beam_config["Patient_model"] = self._patient_model
                beam_config["eyeplan_collimator_points_in_plane"] = self._patient_model.eyeplan_model.collimator_points_in_plane
                #beam_config["aperture_expansion"] = self.tps_config["aperture_expansion"]
                
                
                
                target_voxel_com = np.asarray(ndimage.measurements.center_of_mass(beam_config["target"]))
                beam_config["Slice"] = [int(target_voxel_com[0]),int(target_voxel_com[1]),int(target_voxel_com[2])]

                grid = self._patient_model.structure_set.grid
                beam_config["grid"] = grid
                beam_config["Nvoxels"] = [
                    grid.size[0], grid.size[1], grid.size[2]]
                beam_config["Mesh dimensions"] = [grid.size[0]*grid.spacing[0],
                                                  grid.size[1]*grid.spacing[1],
                                                  grid.size[2]*grid.spacing[2]]
                

                self.add_beam(beam_config)

        return self._beams

    def add_beam(self, beam_config):
        self._beams.append(optis_tps.Beam(beam_config))

    def __iter__(self) -> optis_tps.Beam:
        for beam in self._beams:
            yield beam
            
    def __repr__(self):
        self.title = "Beam set of: \n"
        for beam in self.beams:
            self.title += beam.title + "\n"

        return self.title
