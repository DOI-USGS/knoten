import ctypes
from ctypes.util import find_library
from distutils import sysconfig
import os
import warnings

from . import csm

from csmapi import csmapi

# Register the usgscam plugin with the csmapi

lib = ctypes.CDLL(find_library('usgscsm'))
if not lib:
    warnings.warn('Unable to load usgscsm shared library')
