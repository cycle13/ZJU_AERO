'''
@Description: provides a set of convenient functions that are used
throughout the radar operator
@Author: Hejun Xie
@Date: 2020-07-16 10:06:10
LastEditors: Hejun Xie
LastEditTime: 2020-11-21 17:25:04
'''


# Global imports
import numpy as np
np.seterr(divide='ignore') # Disable divide by zero error
from textwrap import dedent

def get_earth_radius(latitude):
    '''
    Computes the radius of the earth at a specified latitude
    Args:
        latitude: latitude in degrees

    Returns:
        The earth radius in meters
    '''
    a = 6378.1370*1000 # Minimal earth radius (pole) in m
    b = 6356.7523*1000 # Maximal earth radius (equator) in m
    num = ((a** 2 * np.cos(latitude)) ** 2 + (b ** 2 * np.sin(latitude)) ** 2)
    den = ((a * np.cos(latitude)) ** 2+(b * np.sin(latitude)) ** 2)
    return np.sqrt(num / den)


def nansum_arr(x,y, cst = 0):
    """
    Sums up two arrays with possibly different shapes, by padding them before
    with zeros so they have the same dimensions. Ignores NaN values (they
    are treated as zeros)

    Args:
        x: first array
        y: second array
    Returns:
        The summed up array
    """
    x = np.array(x)
    y = np.array(y)

    diff = np.array(x.shape) - np.array(y.shape)
    pad_1 = []
    pad_2 = []
    for d in diff:
        if d < 0:
            pad_1.append((0,-d))
            pad_2.append((0,0))
        else:
            pad_2.append((0,d))
            pad_1.append((0,0))

    x = np.pad(x, pad_1, 'constant', constant_values = cst)
    y = np.pad(y, pad_2, 'constant', constant_values = cst)

    z = np.nansum([x,y],axis=0)
    return z

def sum_arr(x,y, cst = 0):
    """
    Sums up two arrays with possibly different shapes, by padding them before
    with zeros so they have the same dimensions

    Args:
        x: first array
        y: second array
    Returns:
        The summed up array
    """
    diff = np.array(x.shape) - np.array(y.shape)
    pad_1 = []
    pad_2 = []
    for d in diff:
        if d < 0:
            pad_1.append((0,-d))
            pad_2.append((0,0))
        else:
            pad_2.append((0,d))
            pad_1.append((0,0))


    x = np.pad(x,pad_1,'constant',constant_values=cst)
    y = np.pad(y,pad_2,'constant',constant_values=cst)

    z = np.sum([x,y],axis=0)

    return z

def vlinspace(a, b, N, endpoint=True):
    """
        Vectorized equivalent of numpy's linspace

        Args:
            a: list of starting points
            b: list of ending points
            N: number of linearly spaced values to compute between a and b
            endpoint: boolean (optional), if True, the endpoint will be included
                in the resulting series
        Returns:
            A matrix, where every column i is a linearly spaced vector between
                a[i] and b[i]
    """
    a, b = np.asanyarray(a), np.asanyarray(b)
    return a[..., None] + (b-a)[..., None]/(N-endpoint) * np.arange(N)

def nan_cumprod(x):
    """
    An equivalent of np.cumprod, where NaN are ignored

    Args:
        x: a 1D array
    Returns:
        The 1D array of cumulated products,
            ex: [1,2,NaN,5] will return [1,2,2,10]
    """
    x[np.isnan(x)]=1
    return np.cumprod(x)

def nan_cumsum(x):
    """
    An equivalent of np.cumsum, where NaN are ignored

    Args:
        x: a 1D array
    Returns:
        The 1D array of cumulated sums,
            ex: [1,2,NaN,5] will return [1,3,3,8]
    """
    x[np.isnan(x)]=0
    return np.cumsum(x)

def aliasing(v, nyquist):
    """
        Performs aliasing of a radial velocity

        Args:
            v: the velocity to fold
            nyquist: the nyquist velocity
        Returns:
            The folded radial velocity
    """
    theta = (v + nyquist)/(2*nyquist) * np.pi - np.pi/2.
    theta_fold = np.arctan(np.tan(theta))
    v_fold = (theta_fold + np.pi/2)*(2*nyquist)/np.pi - nyquist
    return v_fold

def combine_subradials(list_of_subradials):
    """
    Combines two subradials by combining their variables. The subradials
    should correspond to the same radar profile (theta and phi angles)

    Args:
        list_of_subradials: the list of subradials
    Returns:
        The combined subradial
    """
    x = list_of_subradials[0]
    for l in list_of_subradials:
        if (l.dist_profile == x.dist_profile).all():
            x.values.update(l.values)
        else:
            msg = '''
            Beams are not defined on the same set of coordinates, aborting
            '''
            print(dedent(msg))
            return
    return x

def row_stack(a1, a2):
    """
    Vertically stacks two rows with possibly different shapes,
    padding them with nan before, so they have the same shape

    Args:
        a1: the first row
        a2: the second row
    Returns:
        The array corresponding to the stacked rows
    """
    [N1, M1] = a1.shape
    [N2, M2] = a2.shape

    if N1 > N2:
        a2 = np.pad(a2,((0,0),(0,M1-M2)), mode='constant',
                    constant_values = np.nan)
    elif N2 < N1:
        a1=np.pad(a2,((0,0),(0,M2-M1)), mode='constant',
                  constant_values = np.nan)
    return np.vstack((a1, a2))

