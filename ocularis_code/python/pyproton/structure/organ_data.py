#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from enum import Enum


class OrganType(Enum):
    TARGET = "target"
    OAR = "oar"
    UNKNOWN = "unknown"


class OrganData:
    def __init__(self, name: str, organ_type: OrganType):
        self._name = name
        self._type = organ_type

        self._title = f"{self.name} of type {self.type}"

    @property
    def name(self):
        return self._name

    @property
    def type(self) -> OrganType:
        return self._type

    def __repr__(self):
        return self._title
