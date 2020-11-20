'''
Description: some utilities to dump data structures temporarily.
Author: Hejun Xie
Date: 2020-11-14 11:19:44
LastEditors: Hejun Xie
LastEditTime: 2020-11-14 11:21:48
'''

import os
import pickle
from functools import wraps 


class DATAdecorator(object):
    def __init__(self, workdir, pickle_speedup, pickle_filename):
        self.workdir = workdir
        self.pickle_speedup = pickle_speedup
        self.pickle_filename = pickle_filename
    
    def __call__(self, worker):
        @wraps(worker)
        def wrapped_worker(*args, **kwargs):
            if not self.pickle_speedup \
                or not os.path.exists(self.pickle_filename):
                cdir = os.getcwd()
                os.chdir(self.workdir)
                DATA = worker(*args, **kwargs)
                os.chdir(cdir)
                self.pickle_dump(DATA)
            else:
                DATA = self.pickle_load()           
            return DATA
        return wrapped_worker
            
    def pickle_dump(self, DATA):
        print('Dump data at {}'.format(self.pickle_filename))
        makenewdir(os.path.dirname(self.pickle_filename))
        with open(self.pickle_filename, "wb") as f:
            pickle.dump(DATA, f)

    def pickle_load(self):
        print('Load data at {}'.format(self.pickle_filename))
        with open(self.pickle_filename, "rb") as f:
            DATA = pickle.load(f)
        return DATA

def makenewdir(mydir):
    if not os.path.exists(mydir):
        os.system("mkdir {}".format(mydir))