# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 13:42:52 2022

@author: bjoerk_c
"""


import os, sys, inspect
import logging
from pathlib import Path
from abc import ABC, abstractmethod


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import geometry

class Plotter(ABC):
    def __init__(self, algo):
    
            self.algo = algo    
    
            self.logger = logging.getLogger("Plotter.log")        
            if not self.logger.handlers:
                self.logger.setLevel(logging.INFO)
                self.logger.addHandler(logging.StreamHandler())
    
    
            self.curr_path = Path(os.getcwd())
            self.plot_path = self.curr_path.parent.parent / "plots"      
            if not os.path.isdir(self.plot_path):
                os.makedirs(self.plot_path)
            
            
            self.geo = geometry.OptisGeometry(self.algo.config)
            
    @abstractmethod
    def plot():
        pass
    
    
    
    
    