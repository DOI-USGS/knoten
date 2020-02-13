import ctypes
from ctypes.util import find_library
from distutils import sysconfig
import os
import warnings

from . import csm, bundle

from csmapi import csmapi

# Register the usgscam plugin with the csmapi
libusgscsm_path = find_library('usgscsm')

if not libusgscsm_path:
    warnings.warn('libusgscsm not installed, unable to load shared library.')

libusgscsm = ctypes.CDLL(libusgscsm_path)

if not libusgscsm._name:
    warnings.warn('Unable to load usgscsm shared library')
