<!--
 * @Description: Readme for ZJU_AERO
 * @Author: Hejun Xie
 * @Date: 2020-04-06 20:52:07
 * @LastEditors: Hejun Xie
 * @LastEditTime: 2021-10-21 16:23:33
 -->
# ZJU_AERO
Accurate and Efficient Radar Operator, developed by ZJU - Lei-Bi Lab.  
Supports state of the art simulations of dual-polarimetric radar observable varibales.

## Install
### 1. Install SWIG
This project use `swig` to build C extensions. So check if you have installed swig.
If not, then install it with:
```
    $sudo apt-get install swig
```

### 2. Create conda environment for ZJU_AERO
With Anaconda or Miniconda install, it is recommended to create a new conda environment when using ZJU_AERO or even other packages.  
To create a new environment based on the environment.yml:

```
    $conda env create -f environment.yml
```

Install the radar ploting package to plot the simulation results, 
Currently ZJU_AERO supports those two packages:
* `arm_pyart`, please see https://github.com/ARM-DOE/pyart;
* `pycwr`, please see https://github.com/YvZheng/pycwr;  
For example, you can install arm_pyart with
```
    $conda install -c conda-forge arm_pyart
```

### 3. Install the ZJU_AERO packge
Install from downloaded release package ZJU_AERO_x.y.z.tar.gz with:

```
    $tar -zxvf ZJU_AERO_<x.y.z>.tar.gz
    $cd ZJU_AERO_<x.y.z>
    $python setup.py install
```

## Quick Start

**To get started using ZJU_AERO, you need to first get your LUT and model file right in place**

* **Hydrometeor Back Scattering Look-up Table**  
For test cases, we place the lookup-table under the following directory:  
You can specify the `db_name` in option_files: `microphysics - folder_lut`
```
    ZJU_AERO_<x.y.z>/pathos/lut/<db_name>/lut_SZ_<H>_<Freq>_<Freq>_<mp_scheme>_Level<A/B/C>.nc
```

* **NWP Model Output Including Hydrometeor 3D Fields**  
ZJU_AERO supports two kind of NWP model output: **GRAPES** and **WRF**  
You can specify the `model_name` in option_files: `nwp - name`  
Take WRF as example, we place the wrfout file for test cases under directory:  
```
    ZJU_AERO_<x.y.z>/pathos/WRF/thompson/wrfout_xxxx.nc
```
