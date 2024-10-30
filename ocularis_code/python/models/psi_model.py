# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


import geometry
from pyproton.volume import Grid
from pyproton.structure.organ_data import OrganData, OrganType
from pyproton.structure.array_structure import ArrayStructure
from pyproton.handler import EyeplanStructureHandler
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import math

from scipy.interpolate import interp2d

from utils.vector_utils import angle
from pyproton.utils.concave_hull import ConcaveHull
from utils.transformations import rotation_matrix_from_vectors
from utils.inflections import find_inflexions

import models
import dose_engine

from scipy.spatial.transform import Rotation
import pyvista as pv
import sklearn.cluster
from utils.transformations import mapping
from configuration import constants


currentdir = Path(os.getcwd())
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))


class PSIModel(models.Model):

    def __init__(self, patient_config):
        super().__init__()
        self.patient_config = patient_config
        # self.config = patient_config["config"]
        self._config = None
        self._file_path = self.patient_config["patient_path"]

        self._structure_set = None
        # self._structure_set_gaze_centered = None
        # self._structure_set_clips_registered = None
        # self._structure_set_resampled_registered = None
        # self._doseplane_h = None
        # self._doseplane_v = None
        # self._dose_h = None
        # self._dose_v = None
        # self._dose_h_gaze_centered = None
        # self._dose_v_gaze_centered = None
        # self._dose_h_clips_registered = None
        # self._dose_v_clips_registered = None
        # self._polar_angle = None
        # self._azimuth_angle = None
        # self._target_range_from_ptd = None
        # self._modulation_range_from_ptd = None
        # self._gaze_vector_light = None
        # self._intersection_points = None
        # self._clips = None
        # self._clips_gaze_centered = None
        # self._grid = None
        # self._centre_of_model = None
        # self._relative_skin_pos = None
        # self._proximal_margin = 0
        # self._distal_margin = 0
        self._use_of_wedge = None
        self._treatment_eye = None
        self._fixation_eye = None
        self._pupil_to_centre_distance = None

        # self._clip_transformation = None
        # self._ptd_filename = self.patient_config["eyeplan ptd filename"]

        self.n_clips = 4

        # self._collimator_points_in_plane = []

    @property
    def config(self):
        return self._config
    
    @config.setter
    def config(self, conf):
        self._config = conf
        self.add_pupil_structure()