# -*- coding: utf-8 -*-

'''
@Description: setup for ZJU_AERO
@Author: Hejun Xie
@Date: 2020-04-06 20:54:22
LastEditors: Hejun Xie
LastEditTime: 2020-11-20 17:39:35
'''

from setuptools import setup, Extension, Command
import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

_interpolation_c = Extension("_interpolation_c",
                   ["./ZJU_AERO/interp/interpolation_c.i","./ZJU_AERO/interp/interpolation_c.c"],
                   include_dirs = [numpy_include],
                  )

# setup
setup(  name        = "ZJU_AERO",
        description = "Accurate and Efficient Radar Operator of Zhejiang Univeristy",
        version     = "0.1",
        url         = 'https://gitee.com/zju_bilei_lab/zju_-aero',
        author='Hejun Xie - Zhejiang University',
        author_email='hejun.xie@zju.edu.cn',
        license='GPL-3.0',
        packages=['ZJU_AERO', 'ZJU_AERO/interp', 'ZJU_AERO/radar', 'ZJU_AERO/utils', 'ZJU_AERO/const',
        'ZJU_AERO/db', 'ZJU_AERO/core', 'ZJU_AERO/hydro', 'ZJU_AERO/beam', 'ZJU_AERO/config', 'ZJU_AERO/nwp', 
        'ZJU_AERO/nwp/grapes', 'ZJU_AERO/nwp/wrf'],
        package_data   = {'ZJU_AERO/interpolation' : ['*.o','*.i','*.c','*.so']},
        install_requires=['numpy'],
        zip_safe=False,
        ext_modules = [_interpolation_c]
        )
