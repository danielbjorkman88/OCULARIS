#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 2022

@author: Daniel Björkman
"""

import dose_engine


class SOBP(dose_engine.DepthProfile):
    def __init__(self, xes, yes, modulation, depth_shallowest_peak, depth_deepest_peak):
        self.xes = xes
        self.yes = yes
        self.modulation = modulation

        self.depth_deepest_peak = depth_deepest_peak
        self.depth_shallowest_peak = depth_shallowest_peak

        self.R90 = None
        self.R80 = None
        self.R20 = None
        self.R10 = None

    @property
    def dist_R90_first_peak(self):
        if self.R90 == None:
            raise ValueError

        return self.R90 - self.depth_deepest_peak


    @property
    def dist_R90_last_peak(self):
        if self.R90 == None:
            raise ValueError

        return self.R90 - self.depth_shallowest_peak
