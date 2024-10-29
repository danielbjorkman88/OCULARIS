#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2022-04-02.

@author: Daniel Björkman
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
