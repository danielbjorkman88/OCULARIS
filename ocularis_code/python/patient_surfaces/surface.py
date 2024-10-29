# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:44:14 2021

@author: bjoerk_c
"""

from abc import ABC


class PatientSurface(ABC):
    def __init__(self):
        self.type = None

    def find_ray_interaction(self, ray):
        pass
