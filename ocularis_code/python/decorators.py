# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 09:55:23 2021

@author: bjoerk_c
"""

# Decorators


import time
import functools
from functools import lru_cache, wraps
import numpy as np

# def timer(func):
#     def wrapper(*arg, **kw):
#         before = time.time()
#         func(*arg, **kw)
#         print(f"Executed in {time.time() - before} seconds")
        
#     return wrapper


def timer(func):
    def wrapper(*arg, **kw):
        start = time.perf_counter()
        func(*arg, **kw)
        print(f"Executed in {time.perf_counter() - start} seconds")
        
    return wrapper    






def synchronized(lock):
    """ Synchronization decorator """
    def wrap(f):
        @functools.wraps(f)
        def newFunction(*args, **kw):
            with lock:
                return f(*args, **kw)
        return newFunction
    return wrap




def np_cache(function):
    @lru_cache()
    def cached_wrapper(hashable_array):
        array = np.array(hashable_array)
        return function(array)

    @wraps(function)
    def wrapper(array):
        return cached_wrapper(tuple(array))

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper


# import threading
# lock = threading.Lock()

# @synchronized(lock)
# def do_something():
#   # etc

# @synchronzied(lock)
# def do_something_else():
#   # etc