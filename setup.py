# -*- coding: utf-8 -*-

'''
@Description: setup for ZJU_AERO
@Author: Hejun Xie
@Date: 2020-04-06 20:54:22
LastEditors: Hejun Xie
LastEditTime: 2021-07-18 09:59:24
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

# setup
setup(  name        = "ZJU_AERO",
        description = "Accurate and Efficient Radar Operator of Zhejiang Univeristy",
        version     = "0.1",
        url         = 'https://gitee.com/zju_atmo_radiation/zju_-aero',
        author='Hejun Xie - Zhejiang University',
        author_email='hejun.xie@zju.edu.cn',
        license='GPL-3.0',
        packages=['ZJU_AERO', 'ZJU_AERO/interp', 'ZJU_AERO/radar', 'ZJU_AERO/utils', 'ZJU_AERO/const',
        'ZJU_AERO/db', 'ZJU_AERO/core', 'ZJU_AERO/hydro', 'ZJU_AERO/beam', 'ZJU_AERO/config', 'ZJU_AERO/nwp', 
        'ZJU_AERO/nwp/grapes', 'ZJU_AERO/nwp/wrf'],
        package_data   = {'ZJU_AERO/interpolation' : ['*.o','*.i','*.c','*.so']},
        install_requires=[
        'numpy',
        'scipy',
        'netCDF4',
        'h5py',
        'pyproj',
        'xarray'],
        zip_safe=False,
        ext_modules =  [_interpolation_c, _spectrum_integ]
        )
