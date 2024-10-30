# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


import math
import numpy as np
import warnings

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2, degrees = False):
    
    angle = 0
      
    if (np.isnan(v1[0]) or np.isnan(v1[1])) or (np.isnan(v2[0]) or np.isnan(v2[1])):
        return angle
    
    # try block avoids accidental division by 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  
        angle = math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


    if(v1[0]*v2[1] - v2[1]*v2[0] < 0):   
      angle = -angle
      
      
    if degrees:
        return math.degrees(angle)
      
    return angle

