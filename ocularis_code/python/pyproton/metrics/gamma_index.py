
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
from numba import jit



class GammaIndex:

    def __init__(self, dose1: np.ndarray, dose2: np.ndarray, delta_d: float = 0.005, delta_r: float = 0.005):

        
        assert dose1.shape == dose2.shape
        

