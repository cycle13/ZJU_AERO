'''
Description: Define some types of input from user config file
Author: Hejun Xie
Date: 2020-11-21 17:20:05
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 21:09:37
'''

# Global imports
import numpy as np
from textwrap import dedent

BASIC_TYPES = [float, int, str]

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
    if type_ in [np.int64]:
        classified_type = int
    elif type_ in [np.float, np.float64, np.float128]:
        classified_type = float
    elif type_ in [np.str_]:
        classified_type = str
    return classified_type

class Range(object):
    '''
    The Range class extends the range native python class to floats
    It can be used as a convenient way to test if a given values fall
    within a certain range: ex 2.3 in Range(1.2, 5.3) returns True, but
    1.9 in Range(2.0, 10.0) return False.
    '''
    def __init__(self,x0,x1):
        """
        Args:
            x0: lower bound of the range
            x1: lower bound of the range
        Returns:
            A Range class instance
        """

        # some check for x0 and x1
        if type(x0) != type(x1):
            raise ValueError('range bounds are not of the same type!')
        if x1 <= x0:
            raise ValueError('Lower bound is larger than upper bound!')
        
        # assign attributes
        self.x0 = x0
        self.x1 = x1
        self.type = type(x0)
    def __contains__(self,z):
        return type(z) == self.type and z <= self.x1 and z >= self.x0
    def __str__(self):
        return 'Range of values from {:f} to {:f}'.format(self.x0,self.x1)


class TypeList(object):
    '''
    The TypeList class is used to check if a given array or list is of appropriate
    dimension and type
    '''
    def __init__(self, types, dim=[]):
        """
        Checks if a given array or list has the right type(s) and the right
            dimensions

        Args:
            types: a single python type or a list of types, to be checked for
            dim: a tuple e.g. (3,2), which specifies which shape the checked
                array must have (Optional). When not specified, the checked
                array can have any arbitrary dimensions
        Returns:
            A TypeList class instance, which can be used to check an array
                ex. np.array([1.,2.,'a']) == TypeList([float,str],[3,])
                will return True, but np.array(['a','b']) = TypeList([float])
                will return False
        """

        # deal with single types
        if type(types) != list:
            types = [types]
        
        # check the input types and dims
        if not set(types).issubset(set(BASIC_TYPES)):
            msg = ('''
            One of the specified types is invalid! Must be int, float, str
            ''')
            raise ValueError(dedent(msg))
        if any([d < 0 for d in dim]):
            raise(ValueError('Specified dimension is invalid (<0)!'))

        # assign attributes
        self.types = types
        self.dim = dim
    
    def __eq__(self, array):
        '''
        Define the '==' operator
        ex. np.array([1.,2.,'a']) == TypeList([float,str],dim=[3]) will return True,
        but np.array(['a','b']) == TypeList([float]) will return False.
        '''

        flag = False
        
        # get the dimension of the array
        array = np.array(array)
        dim = array.shape

        # Check if dimensions are ok
        if len(self.dim): # if not arbitrary type is accepted
            flag = all([d1 == d2 for d1, d2 in zip(dim, self.dim)])
        else:
            flag = True
        
        # Check if all types are ok as defined in self.types
        flag *= all([generic_type(v) in self.types for v in array.ravel()])

        return flag

    def __str__(self):
        if self.dim != []:
            msg = 'Array of {}, dimensions = {}'.format(' or '.join([str(t) for t in self.types]), self.dim)
        else:
            msg = 'Array of {}, dimensions = arbitrary'.format(self.types)

        return dedent(msg)
