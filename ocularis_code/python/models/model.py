# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:44:14 2021

@author: bjoerk_c
"""

from abc import ABC

import os, sys
from pathlib import Path

currentdir = Path(os.getcwd())
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))

from pyproton.volume import Grid
from typing import Union

class Model(ABC):
    def __init__(self):
        self._grid = None
        self._structures = None
        self.surfaces = []
    
    def load_files(self):
        pass
    
    @property
    def grid(self) -> Union[None, Grid]:
        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid
        for structure_name in self.name_list:
            self._structures[structure_name].grid = self._grid
            
            