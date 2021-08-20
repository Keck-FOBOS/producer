"""
Miscellaneous package utilities.

.. include:: ../include/links.rst
"""

from itertools import chain, combinations

from IPython import embed 

import numpy

def powerset(iterable, reverse=False):
    """"
    Construct an iterable that steps through all combinations of the
    provided iterable.

    This is pulled from the recipes provided by the itertools
    documentation.

    Examples:
        
        Get all unique combinations of the list [1,2,3]:
        >>> list(powerset([1,2,3]))
        [() (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)]

    Args:
        iterable (iterable):
            An iterable object
        reverse (:obj:`bool`, optional):
            Reverse the order (only roughly) of the iterable by placing
            the longer sequences first.
    
    Returns:
        `itertools.chain`: Iterable object that returns the sequence of
        combinations.
    """
    rng = range(len(iterable)+1)[::-1] if reverse else range(len(iterable)+1)
    return chain.from_iterable(combinations(iterable, r) for r in rng)


def polygon_winding_number(polygon, point):
    """
    Determine the winding number of a 2D polygon about a point.
    
    The code does **not** check if the polygon is simple (no interesecting line
    segments).  Algorithm taken from Numerical Recipes Section 21.4.

    Args:
        polygon (`numpy.ndarray`_):
            An Nx2 array containing the x,y coordinates of a polygon.
            The points should be ordered either counter-clockwise or
            clockwise.
        point (`numpy.ndarray`_):
            One or more points for the winding number calculation.
            Must be either a 2-element array for a single (x,y) pair,
            or an Nx2 array with N (x,y) points.

    Returns:
        :obj:`int`, `numpy.ndarray`_: The winding number of each point with
        respect to the provided polygon. Points inside the polygon have winding
        numbers of 1 or -1; see :func:`point_inside_polygon`.

    Raises:
        ValueError:
            Raised if ``polygon`` is not 2D, if ``polygon`` does not have two
            columns, or if the last axis of ``point`` does not have 2 and only 2
            elements.
    """
    # Check input shape is for 2D only
    if len(polygon.shape) != 2:
        raise ValueError('Polygon must be an Nx2 array.')
    if polygon.shape[1] != 2:
        raise ValueError('Polygon must be in two dimensions.')
    _point = numpy.atleast_2d(point)
    if _point.shape[1] != 2:
        raise ValueError('Point must contain two elements.')

    # Get the winding number
    nvert = polygon.shape[0]
    npnt = _point.shape[0]

    dl = numpy.roll(polygon, 1, axis=0)[None,:,:] - _point[:,None,:]
    dr = polygon[None,:,:] - point[:,None,:]
    dx = dl[...,0]*dr[...,1] - dl[...,1]*dr[...,0]

    indx_l = dl[...,1] > 0
    indx_r = dr[...,1] > 0

    wind = numpy.zeros((npnt, nvert), dtype=int)
    wind[indx_l & numpy.logical_not(indx_r) & (dx < 0)] = -1
    wind[numpy.logical_not(indx_l) & indx_r & (dx > 0)] = 1

    return numpy.sum(wind, axis=1)[0] if point.ndim == 1 else numpy.sum(wind, axis=1)


def point_inside_polygon(polygon, point):
    """
    Determine if one or more points is inside the provided polygon.

    Primarily a wrapper for :func:`polygon_winding_number`, that
    returns True for each point that is inside the polygon.

    Args:
        polygon (`numpy.ndarray`_):
            An Nx2 array containing the x,y coordinates of a polygon.
            The points should be ordered either counter-clockwise or
            clockwise.
        point (`numpy.ndarray`_):
            One or more points for the winding number calculation.
            Must be either a 2-element array for a single (x,y) pair,
            or an Nx2 array with N (x,y) points.

    Returns:
        :obj:`bool`, `numpy.ndarray`: Boolean indicating whether or not each
        point is within the polygon.
    """
    return numpy.absolute(polygon_winding_number(polygon, point)) == 1



