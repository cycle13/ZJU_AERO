###
 # @Description: compile the dynamic lib
 # @Author: Hejun Xie
 # @Date: 2020-08-15 16:22:47
 # @LastEditors: Hejun Xie
 # @LastEditTime: 2020-08-15 16:36:49
### 
#!/bin/bash

fname=$1
swig -python ${fname}.i
gcc -fPIC -c ${fname}.c ${fname}_wrap.c -I/home/xhj/software/miniconda3/envs/rdop/lib/python3.7/site-packages/numpy/core/include/ -I/home/xhj/software/miniconda3/envs/rdop/include/python3.7m/ 
ld -shared ${fname}.o ${fname}_wrap.o -o _${fname}.so
