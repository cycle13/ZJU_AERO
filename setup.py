# -*- coding: utf-8 -*-

'''
@Description: setup for pyCRO
@Author: Hejun Xie
@Date: 2020-04-06 20:54:22
@LastEditors: Hejun Xie
@LastEditTime: 2020-04-06 20:56:37
'''

from setuptools import setup, Extension, Command

# setup
setup(  name        = "pyCRO",
        description = "China Radar Operator written in python",
        version     = "0.1",
        url='https://github.com/Usami-Renko/pyCRO',
        author='Hejun Xie - Zhejiang University',
        author_email='hejun.xie@zju.edu.cn',
        license='GPL-3.0',
        packages=['pyCRO'],
        package_data   = {'pyCRO': []},
        install_requires=['numpy'],
        zip_safe=False
        )
