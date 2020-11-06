# -*- coding: utf-8 -*-

'''
@Description: setup for pyCRO
@Author: Hejun Xie
@Date: 2020-04-06 20:54:22
LastEditors: Hejun Xie
LastEditTime: 2020-11-06 22:36:46
'''

from setuptools import setup, Extension, Command
import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

_interpolation_c = Extension("_interpolation_c",
                   ["./pyCRO/interpolation/interpolation_c.i","./pyCRO/interpolation/interpolation_c.c"],
                   include_dirs = [numpy_include],
                  )

# setup
setup(  name        = "pyCRO",
        description = "China Radar Operator written in python",
        version     = "0.1",
        url='https://github.com/Usami-Renko/pyCRO',
        author='Hejun Xie - Zhejiang University',
        author_email='hejun.xie@zju.edu.cn',
        license='GPL-3.0',
        packages=['pyCRO', 'pyCRO/interpolation', 'pyCRO/radar', 'pyCRO/utilities', 'pyCRO/constants',
        'pyCRO/db', 'pyCRO/core', 'pyCRO/hydrometeors', 'pyCRO/beam', 'pyCRO/config', 'pyCRO/nwp', 
        'pyCRO/nwp/grapes', 'pyCRO/nwp/wrf'],
        package_data   = {'pyCRO/interpolation' : ['*.o','*.i','*.c','*.so']},
        install_requires=['numpy'],
        zip_safe=False,
        ext_modules = [_interpolation_c]
        )
