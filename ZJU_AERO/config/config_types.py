'''
Description: Define some types of input from user config file
Author: Hejun Xie
Date: 2020-11-21 17:20:05
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 17:23:13
'''

# Global imports
import numpy as np
from textwrap import dedent

BASIC_TYPES = [float, int, str]

class Range(object):
    def __init__(self,x0,x1):
        """
        The Range class extends the range native python class to floats
        It can be used as a convenient way to test if a given values fall
        within a certain range: ex 2.3 in Range(1.2,5.3) returns True, but
        1.9 in Range(2.0, 10.0) return False.
        Args:
            x0: lower bound of the range
            x1: lower bound of the range
        Returns:
            A Range class instance
        """

        if type(x0) != type(x1):
            raise ValueError('range bounds are not of the same type!')
        if x1 <= x0:
            raise ValueError('Lower bound is larger than upper bound!')
        self.x0 = x0
        self.x1 = x1
        self.type = type(x0)
    def __contains__(self,z):
        return type(z) == self.type and z <= self.x1 and z >= self.x0
    def __str__(self):
        return 'Range of values from {:f} to {:f}'.format(self.x0,self.x1)

def generic_type(value):
    """
        Gets the type of any input, if the type is a numpy type (ex. np.float),
        the equivalent native python type is returned

        Args:
            value: the value whose type should be returned
        Returns:
            The type of the value
    """
    type_ = type(value)
    if type_ == np.int64:
        type_ = int
    elif type_ == np.float or type_ == np.float64 or type_ == np.float128:
        type_ = float
    elif type_ == np.str_:
        type_ = str
    return type_

# The TypeList class is used to check if a given array or list is of appropriate
# dimension and type

class TypeList(object):
    def __init__(self, types, dim = []):
        """
        Checks if a given array or list has the right type(s) and the right
            dimensions

        Args:
            types : a single python type or a list of types, to be checked for
            dim: a tuple e.g. (3,2), which specifies which shape the checked
                array must have (Optional). When not specified, the checked
                array can have any arbitrary dimensions
        Returns:
            A TypeList class instance, which can be used to check an array
                ex. np.array([1.,2.,'a']) == TypeList([float,str],[3,])
                will return True, but np.array(['a','b']) = TypeList([float])
                will return False
        """

        if type(types) != list:
            types = [types]
        # Note that dim = [], is used to indicate an arbitrary length
        if list(set(types)-set(BASIC_TYPES)):
            msg = ('''
            One of the specified types is invalid! Must be int, float, str'
            ''')
            raise ValueError(dedent(msg))
        if any([d<0 for d in dim]):
            raise(ValueError('Specified dimension is invalid (<0)!'))

        self.types = types
        self.dim = dim
    def __eq__(self, array):
        flag = False
        try:
            array = np.array(array)
            dim = array.shape
            # Check if dimensions are ok
            if len(self.dim): # Check only if dim != []
                flag = all([d1 == d2 for d1,d2 in zip(dim,self.dim)])
            else:
                flag = True
            # Check if all types are ok
            flag *= all([generic_type(v) in self.types for v in array.ravel()])
        except:
            pass
        return flag

    def __str__(self):
        if self.dim != []:
            msg = 'Array of {:s}, with dimensions {:s}'.format(self.types,self.dim)
        else:
            msg = 'Array of {:s}, with arbitrary dimensions'.format(self.type_)

        return dedent(msg)
