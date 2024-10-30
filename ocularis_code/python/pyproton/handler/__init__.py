#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


from .base_handler import BaseHandler

# from .DicomModalityHandler import DicomModalityHandler
from .dicom_structure_handler import DicomStructureHandler
from .dicom_eye_structure_handler import DicomEyeStructureHandler
from .eyeplan_structure_handler import EyeplanStructureHandler