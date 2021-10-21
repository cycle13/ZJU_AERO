#!/bin/bash

rm -rf build
mkdir build
cat install.log | xargs rm -rf
python setup.py install --record install.log

