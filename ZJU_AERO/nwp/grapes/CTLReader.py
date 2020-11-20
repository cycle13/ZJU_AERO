# -*- coding: utf-8 -*- 
"""
Read Grads data
2018.03.11
@author: wenqiushi,modified by wanghao
"""
#import pandas as pd
import numpy as np
import datetime
import re
import os

NUMBER = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
'''
   此模块适用于所有时次的数据存储在一个文件中。
   若数据形式为多个时次数据存放于多个对应的文件中,需要先处理成一个文件
'''

class CTLReader(object):
    def __init__(self,ctlfilename,varname=-1):
        self.variables = {}
        self.dimensions = {}
        self.attributes = {}
        self.crement = {}
        self.ctlname = ctlfilename
        self.varname = varname
        
        with open(self.ctlname,'r') as f:
            self.ctl = f.read()

        p = re.compile("%s\s+(.*)" % ('dset'))
        m = p.search(self.ctl)
        if m.group(1)[0] == '^':
            path = os.path.dirname(ctlfilename)
            self.filename = path + os.sep + m.group(1)[1:]
        else:
            self.filename = m.group(1)[:]
        if not os.path.exists(self.filename):
            raise IOError('{} is not exist.'.format(self.filename))

        self._read_dimensions() #获取ctl中的维度信息
        self._read_data(self.varname)

    def _read_dimensions(self):
        if 'xdef' in self.ctl:
            p = re.compile("%s\s+(\d+)\s+linear\s+(%s)\s+(%s)" % ('xdef',NUMBER,NUMBER))
            m = p.search(self.ctl)
            self.variables['longitude'] = np.linspace(float(m.group(2)),
                                                      float(m.group(2))+float(m.group(3))*(int(m.group(1))-1),
                                                      int(m.group(1)))
            self.dimensions['longitude'] = int(m.group(1))
            self.crement['longitude'] = float(m.group(3))

        if 'ydef' in self.ctl:
            p = re.compile("%s\s+(\d+)\s+linear\s+(%s)\s+(%s)" % ('ydef',NUMBER,NUMBER))
            m = p.search(self.ctl)
            self.variables['latitude'] = np.linspace(float(m.group(2)),
                                                      float(m.group(2))+float(m.group(3))*(int(m.group(1))-1),
                                                      int(m.group(1)))
            self.dimensions['latitude'] = int(m.group(1))
            self.crement['latitude'] = float(m.group(3))

        if 'zdef' in self.ctl:
            p = re.compile("%s\s+(\d+)\s+(\w+)" % ('zdef'))
            m = p.search(self.ctl)
            if m.group(2) == 'levels':
                p = re.compile("%s\s+(\d+)\s+levels([\s\S]+)tdef" % ('zdef'))
                m = p.search(self.ctl)
                self.variables['levels'] = np.fromstring(m.group(2),sep='\n')
                self.dimensions['levels'] = len(self.variables['levels'])
            if m.group(2) == 'linear':
                p = re.compile("%s\s+(\d+)\s+(\w+)\s+(%s)\s+(%s)" % ('zdef',NUMBER,NUMBER))
                m = p.search(self.ctl)
                self.variables['levels'] = np.arange(int(m.group(3)),int(m.group(3))+int(m.group(1))*int(m.group(4)),int(m.group(4)))
                self.dimensions['levels'] = int(m.group(1))

        if 'tdef' in self.ctl:
            p = re.compile("%s\s+(\d+)\s+linear\s+(\w+)\s+(\w+)\s" % ('tdef'))
            m = p.search(self.ctl)

            times = []
            initime = datetime.datetime.strptime(m.group(2), '%Hz%d%b%Y')
            self.dimensions['time'] = int(m.group(1))
#           if 'mn' in m.group(3):
#               increment = (int(m.group(1))-1)*datetime.timedelta(minutes=int(re.sub("\D", "", m.group(3))))
#           else:
#               increment = (int(m.group(1))-1)*datetime.timedelta(hours=int(re.sub("\D", "", m.group(3))))
            if 'mn' in m.group(3):
                increment = datetime.timedelta(minutes=int(re.sub("\D", "", m.group(3))))
            else:
                increment = datetime.timedelta(hours=int(re.sub("\D", "", m.group(3))))
            
            self.crement['time'] = increment 

            for i in range(0,self.dimensions['time']):
                times.append(initime+increment*i)

            self.variables['time'] = times
            #print self.variables['time'], self.dimensions['time']
            #print 'initime',',','endtime',',','increment',',','dimensions'
            #print initime,',',endtime,',',increment,',',self.dimensions['time']

    def _read_data(self,varname):   #分析CTL中必要的说明变量
        undef = eval(re.search('undef (%s)' % NUMBER, self.ctl).group(1))  # 缺省值
        if not bool(re.search('options.*endian',self.ctl,flags=re.I)):  # CTL文件中需要包含endian信息
            big_endian = 'flase'
#           print "please check the ctl, you need add endian information."
#           exit()

        big_endian = bool(re.search('options.*big_endian',self.ctl,flags=re.I)) # big_endian or little_endian
        if big_endian: # 暂时解决办法，此种条件下byteswap()不起作用
            rflag = '>f4'  # 大端读法
        else:
            rflag = '<f4'  # 小端读法
        
        sequential = bool(re.search('options.*sequential',self.ctl,flags=re.I)) # sequential or direct
        if not sequential:
            place_hold = 0
        else:
            place_hold = 2

        allvar,dim,long_name = [],[],[]  # 生成所有变量及变量对应层次的列表

        read = False  #识别是否为目标变量的开关

        for line in self.ctl.split('\n'):
            if line.startswith('endvars'):
                read = False
            if read:
                p = re.compile('(\w+)\s+(\d+)\s+(\d+)\s+(.*)')    #目标变量行的正则范式
                m = p.match(line)
                allvar.append(m.group(1))
                dim.append(int(m.group(2)))
                long_name.append(m.group(4))

            if line.startswith('var'):
                read = True

        if self.varname == -1:
            self.varname = allvar

        for i in range(len(dim)): # 将层次为0的转化为1，后边计算需要
            if dim[i] == 0:
                dim[i] = 1

        for ivarname in self.varname:  # 此段代码来读取相应变量的数据
            var = self.variables[ivarname] = Variable(ivarname)       #生成特定的变量类并在本段方法中以"var"的别名进行描述
            index = allvar.index(ivarname)
            long_name_tmp = long_name[index]
            var.dimensions['levels'] = dim[index]
            var.dimensions['time'],var.dimensions['latitude'],var.dimensions['longitude'] = self.dimensions['time'],self.dimensions['latitude'],self.dimensions['longitude']
            if var.dimensions['levels'] == 0:  # 当读到CTL中对应的层次为0时，代表数据只有一层
                var.dimensions['levels'] = 1
            var.variables['levels'] = self.variables['levels'][0:var.dimensions['levels']]
            SPACE = self.dimensions['latitude']*self.dimensions['longitude']
            size = var.dimensions['levels']*(SPACE+place_hold)

            #var.dimensions_name = ('time','levels','latitude','longitude')    #当变量为四维数组或三维数组时变量的维度信息用同一方式表示，方便主程序取值写法统一
            var.dimensions_name = ('levels','latitude','longitude')

            data = np.zeros((self.dimensions['time'],var.dimensions['levels'],self.dimensions['latitude'],self.dimensions['longitude']))

            var.shape = tuple(var.dimensions[dim] for dim in var.dimensions_name)    #根据不同的维度信息创建维度宽度提示元组

            filename = open(self.filename,'rb')        # 打开文件
            
            for it in np.arange(0,self.dimensions['time']):
                for iz in np.arange(0,sum(dim[0:index])):  #'f4'
                    tmpdata = np.fromfile(filename,dtype=rflag,count=SPACE+place_hold)  ##跳过不需要的变量

                if sequential:        # 处理顺序存取文件
                    data[it,:,:,:] = np.fromfile(filename,dtype=rflag,count=size).reshape(-1,SPACE+place_hold)[:,
                                                                                                             int(place_hold/2):
                                                                                                            -int(place_hold/2)].reshape(var.shape)
                else:                 # 处理直接存取文件
                    data[it,:,:,:] = np.fromfile(filename,dtype=rflag,count=size).reshape(-1,SPACE+place_hold).reshape(var.shape)
                #   print data

                if self.dimensions['time'] != 1:  # 当时次不是1时，还需循环跳过后面的变量
                    for iz in np.arange(0,sum(dim[index+1:])):
                        tmpdata = np.fromfile(filename,dtype=rflag,count=SPACE+place_hold)  ##跳过不需要的变量

        #   if big_endian:  # 如果是big_endian数据，进行转换  # 目前情况下并不能转换，原因未知
        #       data = data.byteswap()
            var.attributes = {
                'long_name' : long_name_tmp
            }

            var.data = data

            filename.close()

class Variable(object):    #变量类定义
    def __init__(self,name,data=None):    #创世纪
        self.name = name                  #python说：“要有名字“！于是有了变量
        self.data = data                  #python说：”要有数据“！于是有了变量
        self.variables = {} 
        self.dimensions = {}    
    def __getitem__(self,index):
        return self.data[index]
    def __getattr__(self,key):
        return self.attributes[key]
