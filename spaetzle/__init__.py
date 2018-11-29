import ctypes
from distutils import sysconfig
import os
import warnings

from csmapi import csmapi

# Register the usgscam plugin with the csmapi

lib = ctypes.CDLL(os.path.abspath(os.path.join(sysconfig.get_python_lib(), '../../libusgscsm.so')))
if not lib:
    warnings.warn('Unable to load usgscsm shared library')
