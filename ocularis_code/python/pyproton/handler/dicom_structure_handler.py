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
from pyproton.handler.dicom import dicom_to_dict

from pyproton.structure.array_structure import ArrayStructure
from pyproton.structure.organ_data import OrganData, OrganType
from pyproton.structure.structure_set import StructureSet

class DicomStructureHandler(BaseHandler):
    def __init__(self, filename: str):
        """
        Parameters
        ----------
        filename: str
            Path to a DICOM file containing the structure set.
        """
        self._dicom_dict = None
        self._structure_names = None
        self._pointsets = None
        self._structure_colors = None

        self._filename = filename
        self._structureset = None

    @property
    def structure_set(self) -> StructureSet:
        if self._structureset is None:
            self._load_file()
            self._structureset = self._build_structure_set()
        return self._structureset

    def _load_file(self):
        self.dicom_dataset = pydicom.read_file(self._filename)
        self._dicom_dict = dicom_to_dict(self.dicom_dataset, recursive=True)

        # Collect names of the structures (regions of interest).
        self._structure_names = []
        for roi in self._dicom_dict["Structure Set ROI Sequence"]:
            self._structure_names.append(roi["ROI Name"])

        # Collect the points belonging to each structure and their color.
        self._pointsets = {}
        self._structure_colors = {}
        for ID, name in enumerate(self._structure_names):
            color = None
            all_points = []
            if "Contour Sequence" in self._dicom_dict["ROI Contour Sequence"][ID]:
                contour_sequence = self._dicom_dict["ROI Contour Sequence"][ID][
                    "Contour Sequence"
                ]

                color = self._dicom_dict["ROI Contour Sequence"][ID].get(
                    "ROI Display Color", ""
                )

                # Collect the points in each slice.
                for contour in contour_sequence:
                    n_points = contour["Number of Contour Points"]

                    points = np.empty((n_points, 3))
                    for i, val in enumerate(contour["Contour Data"]):
                        points[i // 3, i % 3] = val

                    all_points.append(points)

            self._structure_colors[name] = np.asarray(color)/255
            self._pointsets[name] = all_points

    def _build_structure_set(self):
        structures = []
        for name, points in self._pointsets.items():
            organ = OrganData(name, OrganType.UNKNOWN)
            
            if type(points) == list:
                points = np.vstack(points)

            struct = ArrayStructure(organ, points)
            struct.color = self._structure_colors[name]
            structures.append(struct)
        return StructureSet(structures)