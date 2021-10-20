<!--
 * @Description: Readme for ZJU_AERO
 * @Author: Hejun Xie
 * @Date: 2020-04-06 20:52:07
 * @LastEditors: Hejun Xie
 * @LastEditTime: 2021-10-20 20:33:57
 -->
# ZJU_AERO
Accurate and Efficient Radar Operator, developed by ZJU - Lei-Bi Lab.
Supports state of the art simulations of dual-polarimetric radar observable varibales.

## Install
1. This project use swig to build C extensions. So check if you have installed swig.
If not, then install it with:
```
    $sudo apt-get install swig
```

2. With Anaconda or Miniconda install, it is recommended to create a new conda environment when using ZJU_AERO or even other packages. 
To create a new environment based on the environment.yml:

```
    $conda env create -f environment.yml
```

Install the radar ploting package to plot the simulation results, 
Currently ZJU_AERO support those two packages:
1) arm_pyart, please see https://github.com/ARM-DOE/pyart;
2) pycwr, please see https://github.com/YvZheng/pycwr;
For example, you can install arm_pyart with
```
    $conda install -c conda-forge arm_pyart
```

3.  Install the ZJU_AERO packge

1) Install from downloaded release package ZJU_AERO_x.y.z.tar.gz;
'''
    $tar -zxvf ZJU_AERO_<x.y.z>.tar.gz
    $cd ZJU_AERO_<x.y.z>
    $python setup.py install
'''

## Get started
1. To get started using ZJU_AERO, you need to first get your LUT and model file right in place

1) Hydrometeor back scattering look-up table

2) NWP Model output including hydrometeor 3D fields
