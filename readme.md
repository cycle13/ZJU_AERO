<!--
 * @Description: Readme for ZJU_AERO
 * @Author: Hejun Xie
 * @Date: 2020-04-06 20:52:07
 * @LastEditors: Hejun Xie
 * @LastEditTime: 2021-10-23 16:44:36
 -->
# ZJU_AERO (Accurate and Efficient Radar Operator)
Supports state-of-the-art simulations of dual-polarimetric radar observable varibales.

 - [Developers and contributors](CONTRIBUTORS.txt)

## Features
1. `Non-spherical Hydrometeor Particle Shapes`: Including spheroid, hexagonal column, superformula snowflake with different shape parameters... Nonspherical and inhomogeneous particle optical data base are computed with invariant-imbedding T-matrix.
2. `Orientation Perference`: Solid particle orientation perference data stems from snow camera observation (MASC)
3. `New Mass-Diameter Scheme`: Adapted for non-spherical shape hydrometeor database.
4. `Beam Broadening`: Implemented by beam-sampling method.
5. `Beam Bending`:  Use online beam propogation method, which accounts for atmosphere refraction index derived from moisture and temperature.
6. `(Partial) Terrain Block`: Terrian Blocking Effect can be simulated with ZJU_AERO by beam-sampling and propogation algorithm.
7. `Path Attenuation`: ZJU_AERO can simulate the attenuation effect in the path of radar beam.
8. `Spaceborne and Ground-Based Radar In One Operator`: Users can apply ZJU-AERO for spaceborne radar (like GPM-DPR) simulation and ground based radar simulation.


*References*:
* *Bi, Lei, and Ping Yang. "Accurate simulation of the optical properties of atmospheric ice crystals with the invariant imbedding T-matrix method." Journal of Quantitative Spectroscopy and Radiative Transfer 138 (2014): 17-35.*
* *Garrett, Timothy J., et al. "Orientations and aspect ratios of falling snow." Geophysical Research Letters 42.11 (2015): 4617-4622.*

## Install From Release Package
Currently ZJU_AERO only supports installing from release package. The more user-friendly method of conda package installation will be available soon.
### 1. Install GCC
First check you have [GCC](https://gcc.gnu.org/) installed on your machine. It is **mandatory** and should be used together with SWIG to build C-extentions for ZJU_AERO.
```
    $gcc -v
```
If not, try install it with:
```
    $sudo apt-get install gcc
```

### 2. Install SWIG
SWIG is a software development tool that connects programs written in C and C++ with a variety of high-level programming languages, like python. For more information on this topic, please view [http://www.swig.org/](http://www.swig.org/).  
ZJU_AERO builds C extensions with `SWIG`, to speed up the radar operator simulations.
So check the `SWIG` version on your machine.  
```
    $swig -version
```
Compilation tests have been performed for swig==3.x.x / 4.x.x. So if you do not have swig, or the swig on your machine is obsolete, install it with:
```
    $sudo apt-get install swig
```

### 3. Create Conda Environment For ZJU_AERO
We recommend [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for environment management of ZJU_AERO:

So if you do not have anaconda or miniconda installed on your machine, install it with:
```
    $wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $chmod +x Miniconda3-latest-Linux-x86_64.sh
    $./Miniconda3-latest-Linux-x86_64.sh
```


With Anaconda or Miniconda install, it is recommended to create a new conda environment when using ZJU_AERO or even other packages.  
To create a new environment based on the environment.yml:

```
    $conda env create -f environment.yml
```

Here we list the dependencies in environment.yml. You can check it by yourself:
* [numpy](https://github.com/numpy/numpy)
* [scipy](https://github.com/scipy/scipy)
* [xarray](https://github.com/pydata/xarray)
* [netcdf4-python](https://github.com/Unidata/netcdf4-python)
* [h5netcdf](https://github.com/h5netcdf/h5netcdf)
* [h5py](https://github.com/h5py/h5py)
* [pandas](https://github.com/pandas-dev/pandas)
* [pyproj](https://github.com/pyproj4/pyproj)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [pyyaml](https://github.com/yaml/pyyaml)
* [multiprocess](https://github.com/uqfoundation/multiprocess)
* [future](https://github.com/PythonCharmers/python-future)
* [basemap](https://github.com/matplotlib/basemap)

To plot the simulation results, you need to install a radar ploting package, 
Currently ZJU_AERO supports two packages:
* `arm_pyart`, please see https://github.com/ARM-DOE/pyart;
* `pycwr`, please see https://github.com/YvZheng/pycwr;

### 4. Install The ZJU_AERO Packge
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
ZJU_AERO supports two kinds of NWP model output: **GRAPES** and **WRF**  
You can specify the `model_name` in option_files: `nwp - name`  
Take WRF as example, we place the wrfout file for test cases under directory:  
```
    ZJU_AERO_<x.y.z>/pathos/WRF/thompson/wrfout_xxxx.nc
```

* **Example1: PPI Scan Simulation**  
    * Test Script: `example/ppi.py`
    * User Option File: `example/option_files/example.yml`
    * WRF Model File: `pathos/WRF/thompson/wrfout_d02_2021-08-08_00_00_00`
    * LUT for Grauel, Cloud Ice, Snow and Rain: 
      * `pathos/lut/tm_masc_release/lut_SZ_G_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_I_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_S_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_R_9_41_1mom_LevelB.nc`
    ```
    $ cd example
    $ python ppi.py
    ```

    * Plotted PPI Scan `Reflectivity`-`ZH`
    ![avatar](pictures/test_ppi_ZH.png)

    * Plotted PPI Scan `Differential Reflectivity`-`ZDR`
    ![avatar](pictures/test_ppi_ZDR.png)

    * Plotted PPI Scan `Radial Velocity`-`VR`
    ![avatar](pictures/test_ppi_RVEL.png)

* **Example2: RHI Scan Simulation**  
    * Test Script: `example/rhi.py`
    * User Option File: `example/option_files/example.yml`
    * WRF Model File: `pathos/WRF/thompson/wrfout_d02_2021-08-08_00_00_00`
    * LUT for Grauel, Cloud Ice, Snow and Rain: 
      * `pathos/lut/tm_masc_release/lut_SZ_G_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_I_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_S_9_41_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_R_9_41_1mom_LevelB.nc`
    ```
    $ cd example
    $ python rhi.py
    ```

    * Plotted RHI Scan `Reflectivity`-`ZH`
    ![avatar](pictures/test_rhi_ZH.png)

    * Plotted RHI Scan `Differential Reflectivity`-`ZDR`
    ![avatar](pictures/test_rhi_ZDR.png)

    * Plotted RHI Scan `Radial Velocity`-`VR`
    ![avatar](pictures/test_rhi_RVEL.png)

* **Example3: Spaceborne Scan Simulation**  
    * Test Script: `example/spaceborne.py`
    * User Option File: `example/option_files/example_spaceborne.yml`
    * WRF Model File: `pathos/GRAPES/typhoon_haishen_20200905/modelvar202009050000900.nc`
    * LUT for Grauel, Cloud Ice, Snow and Rain: 
      * `pathos/lut/tm_masc_release/lut_SZ_G_13_6_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_I_13_6_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_S_13_6_1mom_LevelB.nc`
      * `pathos/lut/tm_masc_release/lut_SZ_R_13_6_1mom_LevelB.nc`
    ```
    $ cd example
    $ python spaceborne.py
    ```

    * Simulated Spaceborne Scan `Reflectivity`-`ZH` at 0km altitude
    ![avatar](pictures/test_spaceborne_0km.png)

    * Observed Spaceborne Scan `Reflectivity`-`ZH` at 8km altitude
    ![avatar](pictures/gpm_dpr_swath_0km.png)

    * Simulated Spaceborne Scan `Reflectivity`-`ZH` at 0km altitude
    ![avatar](pictures/test_spaceborne_8km.png)

    * Observed Spaceborne Scan `Reflectivity`-`ZH` at 8km altitude
    ![avatar](pictures/gpm_dpr_swath_8km.png)

For more Details and examples, please read [User_Guide](doc/User_Guide-ZJU_AERO-0.1.4.pdf) 

Project development plan
----------

- [ ] Melted Ice simulation
- [x] Thompson Microphysics constants
- [ ] Online Documentation
- [ ] Upload ZJU_AERO as a conda package
 