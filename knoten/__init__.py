import ctypes
from distutils import sysconfig
import os
import warnings

from csmapi import csmapi

# Register the usgscam plugin with the csmapi

lib = ctypes.CDLL(ctypes.util.find_library('usgscsm.so'))
if not lib:
    warnings.warn('Unable to load usgscsm shared library')
