#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
import pydicom
from pyproton.handler.base_handler import BaseHandler
from pyproton.handler.base_handler import BaseHandler
from pyproton.handler.dicom_structure_handler import DicomStructureHandler

from pyproton.structure.array_structure import ArrayStructure
from pyproton.structure.organ_data import OrganData, OrganType
from pyproton.structure.eye_structure_set import EyeStructureSet


class DicomEyeStructureHandler(DicomStructureHandler):

    def _build_structure_set(self):
        structures = []
        for name, points in self._pointsets.items():
            organ = OrganData(name, OrganType.UNKNOWN)

            if type(points) == list:
                points = np.vstack(points)

            struct = ArrayStructure(organ, points)
            struct.color = self._structure_colors[name]
            structures.append(struct)
        return EyeStructureSet(structures)
