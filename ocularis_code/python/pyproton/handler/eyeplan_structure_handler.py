#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2022-05-20.

@author:, Daniel Björkman
"""
import numpy as np
import os
import pandas as pd
from pathlib import Path

from pyproton.handler.base_handler import BaseHandler


from pyproton.structure.array_structure import ArrayStructure
from pyproton.structure.organ_data import OrganData, OrganType
from pyproton.structure.structure_set import StructureSet


class EyeplanStructureHandler(BaseHandler):
    def __init__(self, file_path: str):
        """
        Parameters
        ----------
        filename: str
            Path to a EYEPLAN files containing the structure set.
        """
        self._file_path = file_path
        self._structure_names = None
        self._pointsets = None
        self._structureset = None
        self._loaded_filenames = []

    @property
    def structure_set(self) -> StructureSet:
        if self._structureset is None:
            self._load_file()
            self._structureset = self._build_structure_set()
        return self._structureset

    def _load_file(self):

        files = os.listdir(self._file_path)

        # xlsx_files = list(filter(lambda x: x[-5:] == ".xlsx", files))

        csv_files = list(filter(lambda x: x[-4:] == ".csv", files))

        # for idx in range(len(csv_files)):
        #     if csv_files[idx].split("_")[1] == "scleralR.csv":
        #         csv_files.pop(idx)
        #         break

        # for idx in range(len(csv_files)):
        #     split_line = csv_files[idx].split("_")
        #     if len(split_line) == 3:
        #         if csv_files[idx].split("_")[1] == "trgt" and csv_files[idx].split("_")[2] == "baseArea.csv":
        #             csv_files.pop(idx)
        #             break

        self._pointsets = {}
        self._structure_names = []

        # for filename in xlsx_files:
        #     structure_name = filename.split("_")[1][0:-5]
        #     self._structure_names.append(structure_name)
        #     self._pointsets[structure_name] = np.array(pd.read_excel(
        #         Path(self._file_path) / filename))

        for filename in csv_files:
        
            if len(filename.split("_")) == 3:
                # for 'trgt_base.csv'
                structure_name = filename[7:-4].strip("_")
            else:
                structure_name = filename.split("_")[-1][0:-4].strip("_")
            self._structure_names.append(structure_name)
            self._pointsets[structure_name] = np.array(pd.read_csv(
                Path(self._file_path) / filename, header=None)).T
            self._loaded_filenames.append(Path(self._file_path) / filename)

        # for name in self._structure_names:
        #     print(name + ".xlsx")
        #     self._pointsets[name] = np.array(pd.read_excel(
        #         Path(self._file_path) / (name + ".xlsx")))

    def _build_structure_set(self):
        structures = []
        for name, points in self._pointsets.items():

            organ = OrganData(name, OrganType.UNKNOWN)

            if points.shape == (1, 1):
                points = np.asarray([points[0], points[0], points[0]])

            if type(points) == list:
                points = np.vstack(points)

            if points.shape[1] == 1:
                points = points.T

            if points.shape[1] == 4:
                points = np.delete(points, 3, 1)

            struct = ArrayStructure(organ, points)
            structures.append(struct)
        return StructureSet(structures)
