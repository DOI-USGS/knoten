import ctypes
from ctypes.util import find_library
from glob import glob
import os
import warnings

# Register the usgscam plugin with the csmapi
libusgscsm_path = find_library('usgscsm')
    
if not libusgscsm_path:
    libcsmapi_path = find_library('csmapi')
    usgscsm_folder = os.path.join(os.path.split(libcsmapi_path)[0], "csmplugins")
    libusgscsm_path = ""
    if os.path.exists(usgscsm_folder):
        # Supports py < 3.10, if only supporting 3.10+ use: glob( "*[0-9].[0-9].[0-9].dylib", root_dir=usgscsm_folder)
        results = glob(os.path.join(usgscsm_folder, "*[0-9].[0-9].[0-9].dylib"))
        results.sort()
        libusgscsm_path = os.path.join(usgscsm_folder, results[-1])

if not os.path.exists(libusgscsm_path):
    warnings.warn('libusgscsm not installed, unable to find shared library.')

try:
    libusgscsm = ctypes.CDLL(libusgscsm_path)
except OSError:
    warnings.warn(f'Unable to load usgscsm shared library at {libusgscsm_path}')
