#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 08:50:42 2021

@author: Daniel Björkman
"""

import dose_engine


class BraggPeak(dose_engine.DepthProfile):
    def __init__(self, xes, yes, w, modulation):
        self.xes = xes
        self.yes = yes
        self.w = w
        self.modulation = modulation

        self.relative_contribution = []
        self.R90 = []
        self.R80 = []
        self.R20 = []
        self.R10 = []
