#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 08:50:42 2021

@author: Daniel Björkman

This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import os
import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path


class DepthProfile:
    def __init__(self):

        self.data_path = []
        self.filename = []
        self.R90 = []
        self.R80 = []
        self.R20 = []
        self.R10 = []

        # 90 and 80 % of maximum at proximal side of Bragg peak
        self.P90 = []
        self.P80 = []

        self.title = []

    def extrapolate(self):

        xes = self.xes

        yes = self.yes

        f = interp1d(xes, yes, fill_value='extrapolate')

        new_x = -3
        new_y = f(new_x)

        new_xes = np.append(new_x, xes)
        new_yes = np.append(new_y, yes)

        new_xes = new_xes - new_x

        self.xes = new_xes
        self.yes = new_yes

    def calc(self):

        f = interp1d(self.xes, self.yes)
        self.f = f

        f_rev = interp1d(self.yes, self.xes)

        self.P90 = f_rev(max(self.yes)*0.9).max()
        self.P80 = f_rev(max(self.yes)*0.8).max()

        subx = self.xes[self.xes > 0.99*f_rev(max(self.yes))]
        suby = self.f(subx)

        self.subx = subx
        self.suby = suby

        f_sub = interp1d(suby, subx)

        self.R100 = f_sub(max(suby)*1).max()
        self.R90 = f_sub(max(suby)*0.9).max()
        self.R80 = f_sub(max(suby)*0.8).max()
        self.R20 = f_sub(max(suby)*0.2).max()
        self.R10 = f_sub(max(suby)*0.1).max()

        self.title = "Falloff 90-10 = {} mm. Falloff 80-20 = {} mm. BP width 90 = {} mm. BP width 80 = {} mm".format(round(
            self.R10 - self.R90, 3), round(self.R20 - self.R80, 3), round(self.R90 - self.P90, 3), round(self.R80 - self.P80, 3))

    def load_scan(self, filename, DD_path=Path(os.getcwd())):

        self.data_path = DD_path
        self.filename = filename

        data = np.loadtxt(self.data_path / filename, skiprows=7)

        xes = data[::-1, 2]
        yes = data[0:, 5]
        yes[yes < 0] = 0

        yes = yes - yes[-1]

        self.label = filename[9:14]

        start = min(xes)
        xes = xes - start

        # Adding padding for simpler interpolation
        xes = np.append(xes, 1000)
        yes = np.append(yes, 0)

        max_value = max(yes)

        self.xes = xes
        self.yes_abs = yes
        self.yes = yes/max_value

        # self.extrapolate()
        self.calc()

    def load_npy(self, filenames, DD_path=Path(os.getcwd())):

        self.data_path = DD_path
        self.filename = filenames

        xes = np.load(self.data_path / filenames[0])
        yes = np.load(self.data_path / filenames[1])

        # Adding padding for simpler interpolation
        xes = np.append(xes, 1000)
        yes = np.append(yes, 0)

        self.xes = xes
        self.yes = yes

        # self.extrapolate()
        self.calc()
