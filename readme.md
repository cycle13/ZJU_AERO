<!--
 * @Description: Readme for ZJU_AERO
 * @Author: Hejun Xie
 * @Date: 2020-04-06 20:52:07
 * @LastEditors: Hejun Xie
 * @LastEditTime: 2020-11-20 17:09:46
 -->
# ZJU_AERO
China Radar Operator written in python

## Install
1.  Install arm_pyart, please see https://github.com/ARM-DOE/pyart;
    or install pycwr, please see https://github.com/YvZheng/pycwr to plot simulation results.

2.  Build your conda enviroment by provided rdop.yaml:
Modify the prefix to the direcory you want ZJU_AERO to install to:
<-- prefix: /home/xhj/software/miniconda3/envs/rdop
--> prefix: /some/directory/to/install/ZJU_AERO

```
    conda env create -f rdop.yaml
```

3.  build and install
```
    python setup.py build
    python setup.py install
```
