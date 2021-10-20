# -*- coding: utf-8 -*-

'''
@Description: setup for ZJU_AERO
@Author: Hejun Xie
@Date: 2020-04-06 20:54:22
LastEditors: Hejun Xie
LastEditTime: 2021-10-19 21:47:51
'''

from setuptools import setup, Extension, Command
import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

_interpolation_c = Extension( "ZJU_AERO/interp/_interpolation_c",
                   ["./ZJU_AERO/interp/interpolation_c.i","./ZJU_AERO/interp/interpolation_c.c"],
                   include_dirs = [numpy_include],
                  )

_spectrum_integ = Extension( "ZJU_AERO/core/_spectrum_integ",
                   ["./ZJU_AERO/core/spectrum_integ.i","./ZJU_AERO/core/spectrum_integ.c"],
                   include_dirs = [numpy_include],
                  )


DISTNAME = "ZJU_AERO"
AUTHOR = 'Hejun Xie - Zhejiang University'
AUTHOR_EMAIL = 'hejun.xie@zju.edu.cn'
URL = "https://github.com/Usami-Renko/ZJU_AERO"
LICENCE = 'MIT'
PYTHON_REQUIRES = ">=3.7"
# INSTALL_REQUIRES = ['numpy', 'scipy', 'xarray', 'netCDF4', 'h5py', 'pandas', 
#                     'pyproj', 'matplotlib', 'pyyaml', 'multiprocess', 'future']
INSTALL_REQUIRES = []
DESCRIPTION = "Accurate and Efficient Radar Operator"
LONG_DESCRIPTION = """Accurate and Efficient Radar Operator, developed by ZJU Lei Bi Lab.
Supports state of the art simulations of dual-polarimetric radar observable varibales.
"""
PLATFORMS = ["Linux"]
CLASSIFIERS = [
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Operating System :: POSIX :: Linux']
PACKAGES = ['ZJU_AERO', 'ZJU_AERO/interp', 'ZJU_AERO/radar', 'ZJU_AERO/utils', 'ZJU_AERO/const',
            'ZJU_AERO/db', 'ZJU_AERO/core', 'ZJU_AERO/hydro', 'ZJU_AERO/beam', 'ZJU_AERO/config', 'ZJU_AERO/nwp', 
            'ZJU_AERO/nwp/grapes', 'ZJU_AERO/nwp/wrf', 'ZJU_AERO/nwp/graph']
PACKAGE_DATA = {'ZJU_AERO/interp' : ['*.i','*.c'],
                'ZJU_AERO/core'   : ['*.i','*.c']}

# setup
setup(  name         = DISTNAME,
        version      = "0.1.4",
        author       = AUTHOR,
        license      = LICENCE,
        author_email = AUTHOR_EMAIL,
        description  = DESCRIPTION,
        long_description  = LONG_DESCRIPTION,
        python_requires   = PYTHON_REQUIRES,
        install_requires  = INSTALL_REQUIRES,
        url               = URL,
        platforms         = PLATFORMS,
        classifiers       = CLASSIFIERS,
        packages          = PACKAGES,
        package_data      = PACKAGE_DATA,
        ext_modules       =  [_interpolation_c, _spectrum_integ]
        )
