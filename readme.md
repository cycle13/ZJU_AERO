<!--
 * @Description: Readme for pyCRO
 * @Author: Hejun Xie
 * @Date: 2020-04-06 20:52:07
 * @LastEditors: Hejun Xie
 * @LastEditTime: 2020-09-27 21:33:55
 -->
# pyCRO
China Radar Operator written in python

## Install
1.  Install pyWRF, please see https://github.com/Usami-Renko/pyWRF.

2.  Install arm_pyart, please see https://github.com/ARM-DOE/pyart.

3.  Build your conda enviroment by provided rdop.yaml:
Modify the prefix to the direcory you want pyCRO to install to:
<-- prefix: /home/xhj/software/miniconda3/envs/rdop
--> prefix: /some/directory/to/install/pycro

```
    conda env create -f rdop.yaml
```

4.  build and install
```
    python setup.py build
    python setup.py install
```
