#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
from collections.abc import Iterable
import numpy as np
from pyproton.volume.grid import Grid
from pyproton.structure.organ_data import OrganData


class StructureSlice:
    @property
    def n_contours(self) -> int:
        pass

    def __getitem__(self, contour_index: int) -> np.ndarray:
        pass

    @property
    def slice_index(self) -> int:
        pass


class Structure(ABC):
    @property
    @abstractmethod
    def organ(self) -> OrganData:
        raise NotImplementedError()

    @property
    @abstractmethod
    def grid(self) -> Grid:
        raise NotImplementedError()

    @property
    @abstractmethod
    def contour_slices(self) -> Iterable[StructureSlice]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def contour_coordinates(self) -> np.ndarray:
        raise NotImplementedError()

    @property
    @abstractmethod
    def binary_mask(self) -> np.ndarray:
        raise NotImplementedError()

    @property
    @abstractmethod
    def volume(self) -> float:
        raise NotImplementedError()

    @property
    @abstractmethod
    def bounding_box(self):
        # TODO: Specify return type and format.
        raise NotImplementedError()

    @property
    @abstractmethod
    def convex_hull(self):
        # TODO: Specify return type and format.
        raise NotImplementedError()
