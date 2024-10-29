#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
