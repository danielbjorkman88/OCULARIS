import numpy as np
from numba import jit



class GammaIndex:

    def __init__(self, dose1: np.ndarray, dose2: np.ndarray, delta_d: float = 0.005, delta_r: float = 0.005):

        
        assert dose1.shape == dose2.shape
        

