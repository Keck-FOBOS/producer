
import argparse

from IPython import embed

import numpy
from matplotlib import pyplot, patches


def triangle_number(n):
    """
    Return the nth triangle number; i.e., the value of
    ``1+2+3+...+(n-1)+n``.
    """
    return n*(n+1)/2


def nhex(n):
    """
    Return the number of grid points in a hexagonal grid given the
    number of points along its long axis.

    Args:
        n (:obj:`int`):
            Number of points along the long axis of the grid.  Must be
            odd.
    
    Returns:
        :obj:`int`: Number of grid points.
    
    Raises:
        ValueError:
            Raised if the provided integer is even.
    """
    if n % 2 != 1:
        raise ValueError('Number along axis must be odd.')
    return int(n + 2*(triangle_number(n-1) - triangle_number(n//2)))


def nearest_hexgrid(n, upper=False):
    """
    Provided an expected number of grid locations, perform a brute-force
    iteration to find the hexgrid parameters that yield at least ``n`` grid
    locations.

    Algorithm starts with a grid long axis with 3 points and just
    increases the number along the axis until the number of grid points
    is larger than the provided number.  If ``upper`` is True, the
    function treats the provided number as an upper limit on the allowed
    number of grid points.

    Args:
        n (:obj:`int`):
            Number of grid points to match.
        upper (:obj:`bool`, optional):
            Treat the provided number as an upper limit.  If False, the
            method returns the nearest hexagonal grid parameters that
            have at least the provided number of grid points.  If True,
            the hexagonal grid instead must have no more than the
            provided number of grid points.

    Returns:
        :obj:`tuple`: Two integers providing (1) the length of the long
        axis of the hexagonal grid and (2) the number of grid points.
    """
    i = 3
    ng = nhex(i)
    while ng < n:
        i += 2
        ng = nhex(i)
    return (i-1, nhex(i-1)) if upper else (i, ng)


def set_orientation(coo, orientation='horizontal'):
    """
    Set the orientation of the input set of coordinates.

    The input orientation is *assumed* to be horizontal or,
    equivalently, a 0 degree rotation.  I.e., if ``orientation`` is set
    to ``'horizontal'`` or ``0.``, the input coordinates are simply
    returned.  For non-zero rotation angles, repeated calls to this
    function result in additional rotations; i.e., there is no absolute
    reference frame that function can use.

    Args:
        coo (`numpy.ndarray`_):
            Array of Cartesian coordinates to rotate.
        orientation (:obj:`str`, :obj:`float`, optional):
            The orientation to set assuming that the input coordinates
            are provided at 0 rotation (horizontal).  The value of the
            orientation must be either ``'horizontal'``, ``'vertical'``,
            or a rotation angle in degrees.

    Returns:
        `numpy.ndarray`_: Rotated coordinates.
    """
    if isinstance(orientation, str) and orientation not in ['horizontal', 'vertical']:
        raise ValueError('Orientation must be horizontal, vertical, or a rotation in degrees.')
    if orientation == 'horizontal' or orientation == 0.:
        return coo
    return rotate(coo, 30. if orientation == 'vertical' else orientation)


def hexgrid_circle_convert(d, incircle=False):
    """
    Convert between the circumscribe circle and the incircle of a
    hexagon.

    Args:
        d (:obj:`float`):
            The reference hexagonal diameter.  By default this is the
            circumscribed circle (see ``incircle``).
        incircle (:obj:`bool`, optional):
            The provided diameter is the incircle instead of the
            circumscribed circle.

    Returns:
        :obj:`float`: Return the conversion of the input.  If the input
        is the circumscribed circle, this returns the incircle, and vice
        versa.
    """
    return d / numpy.cos(numpy.pi/6) if incircle else d * numpy.cos(numpy.pi/6)


def hexgrid_grid_convert(n, d, incircle=False, d_unit='cell', full=False):
    """
    Convert between the diameter of a full hexagonal grid and each
    hexagonal grid cell.

    This function returns the diameter of the full grid if the diameter
    of the grid cell is provided, and vice versa.  The reference circle
    returned coincides with the one provied; e.g., if the circumcircle
    of the full grid is provided, the circumcircle of an individual grid
    cell is returned.

    Args:
        n (:obj:`int`):
            The number of hexagonal grid cells along the grids long
            axis.
        d (:obj:`float`):
            The reference hexagonal diameter.  By default this is the
            circumscribed circle (see ``incircle``).
        incircle (:obj:`bool`, optional):
            The provided diameter is the incircle of the reference
            hexagon.
        d_unit (:obj:`str`, optional):
            The input unit of the diameter.  There are two options: (1)
            ``'cell'`` means the input diameter is for each grid cell;
            (2) ``'grid'`` means the input diameter is for the entire
            grid.
        full (:obj:`bool`, optional):
            When defining the size of the grid directly (i.e.,
            ``d_unit='grid'``), expand the size of the individual grid
            cells so that they fill the parent hexagon of the grid.
    """
    # Calculate the diameter of the full grid
    if d_unit == 'cell':
        # Calculate the circumcircle of the grid cell
        _d = hexgrid_circle_convert(d, incircle=True) if incircle else d
        # Incircle of the full grid
        grid_d = _d * (n - (n//2)/2)
        # Return the diameter with the definition that matches the input
        return grid_d if incircle else hexgrid_circle_convert(d, incircle=True)

    # Calculate the incircle of the full grid
    _d = d if incircle else hexgrid_circle_convert(d)
    # Get the circumcircle of the grid cel
    cell_d = _d/(n - (n//2)/2 - (numpy.cos(numpy.pi/3) if full else 0.0))
    # Return the diameter with the definition that matches the input
    return hexgrid_circle_convert(cell_d) if incircle else cell_d


def hexgrid(n, d, incircle=False, orientation='horizontal', d_unit='cell', vertices=False,
            close=False, full=False):
    r"""
    Generate the grid centers for a 2D hexagonal grid.

    The grid is generated to match the provided number of points along
    the long axis and diameter of each grid cell.  By default the
    diameter of each grid cell is defined as its circumscribed circle.
    This can be changed to its incircle --- :math:`d_i = d_c
    \cos(\pi/6)`, where :math:`d_i` is the incircle diameter and
    :math:`d_c` is the circumscribed circle diameter --- using
    ``incircle``.

    Each grid cell is separated by its *incircle* diameter; i.e., the
    orientation of each grid cell hexagon (see :func:`hexvert`) is
    vertical for a full grid orientation that is horizontal.  This
    ensures a 100% fill factor for the grid.

    Args:
        n (:obj:`int`):
            The number of hexagonal grid cells along the grids long
            axis.
        d (:obj:`float`):
            The diameter of each grid cell.  By default this is the
            circumscribed circle (see ``incircle``).
        incircle (:obj:`bool`, optional):
            Use the provided diameter to set the incircle of each
            hexagonal cell instead of its circumscribed circle.
        orientation (:obj:`str`, :obj:`float`, optional):
            Sets the orientation of the grid, must be either
            'horizontal', 'vertical', or a rotation angle in degrees
            relative to the horizontal orientation.  The horizontal and
            vertical orientations set the long axis of the grid along
            the Cartesian x and y axis, respectively.  The horizontal
            orientation is equivalent to a rotation angle of 0 and a
            vertical orientation is equivalent to a rotation angle of 30
            degrees.  While the polar-coordinate ordering of the
            vertices in the output array will change, note the shape
            degeneracy periodicity of 60 degrees.
        d_unit (:obj:`str`, optional):
            The input unit of the diameter.  There are two options: (1)
            ``'cell'`` means the input diameter is for each grid cell;
            (2) ``'grid'`` means the input diameter is for the entire
            grid.
        vertices (:obj:`bool`, optional):
            Return the vertices of each grid cell.  If False, only the
            grid centers are returned.
        close (:obj:`bool`, optional):
            If vertices are to be returned for each grid cell, this sets
            whether or not to "close" each cell hexagon by repeating the
            first set of coordinates at the end of the array of
            vertices.
        full (:obj:`bool`, optional):
            When defining the size of the grid directly (i.e.,
            ``d_unit='grid'``), expand the size of the individual grid
            cells so that they fill the parent hexagon of the grid.
    
    Returns:
        `numpy.ndarray`_: An array with shape :math:`(N_{\rm grid},2)`
        providing the x and y Cartesian coordinates for the center of
        each grid cell.
    """
    # Check the input
    if d_unit not in ['cell', 'grid']:
        raise ValueError('Provided diameter unit must be "cell" or "grid".')

    # WARNING: hexvert() depends on this current ordering of the
    # returned grid points in the loop below!!

    # Get the incircle diameter
    _d = d if incircle else hexgrid_circle_convert(d)
    if d_unit == 'grid':
        _d = hexgrid_grid_convert(n, _d, incircle=True, d_unit='grid', full=full)

    # Get the number of grid cells
    ng = nhex(n)

    # Generate the center of each grid cell
    grid = numpy.zeros((ng,2), dtype=float)
    ii = 0
    for i in range(n):
        for j in range(n-numpy.absolute(n//2-i)):
            nn = n//2
            grid[ii,0] = _d * numpy.absolute(nn-i)/2 + _d * (j-nn)
            grid[ii,1] = _d * numpy.sin(numpy.pi/3.0) * (nn-i)
            ii += 1

    if not vertices:
        # Just return the centers
        return set_orientation(grid, orientation=orientation)

    # Get the vertices of a grid cell
    v = hexvert(d=_d, incircle=True, close=close, orientation='vertical')
    return set_orientation(grid[:,None,:] + v[None,:,:], orientation=orientation)


def hexarea(d=1., incircle=False):
    r"""
    Calculate the hexagon area.

    If :math:`d` is the circumcircle diameter, the area is:

    .. math::

        A = \frac{3\sqrt{3}}{8} d^2

    Args:
        d (:obj:`float`, `numpy.ndarray`_, optional):
            The diameter of circumscribed circle.
        incircle (:obj:`bool`, optional):
            The provided diameter is for the incircle of the hexagon, instead of
            its circumscribed circle.

    Returns:
        :obj:`float`, `numpy.ndarray`_: The area of the hexagon(s).
    """
    _d = hexgrid_circle_convert(d, incircle=True) if incircle else d
    return 3*sqrt(3)*_d**2/8


def hexvert(d=1., incircle=False, close=False, orientation='horizontal'):
    r"""
    Construct the vertices of a hexagon.

    Args:
        d (:obj:`float`, optional):
            The diameter of circumscribed circle.
        incircle (:obj:`bool`, optional):
            Use the provided diameter to set the incircle of the hexagon
            instead of its circumscribed circle.
        close (:obj:`bool`, optional):
            "Close" the hexagon by repeating the first set of
            coordinates at the end of the array of vertices.
        orientation (:obj:`str`, :obj:`float`, optional):
            Sets the orientation of the grid, must be either
            'horizontal', 'vertical', or a rotation angle in degrees
            relative to the horizontal orientation.  The horizontal and
            vertical orientations set the long axis of the grid along
            the Cartesian x and y axis, respectively.  The horizontal
            orientation is equivalent to a rotation angle of 0 and a
            vertical orientation is equivalent to a rotation angle of 30
            degrees.  While the polar-coordinate ordering of the
            vertices in the output array will change, note the shape
            degeneracy periodicity of 60 degrees.

    Returns:
        `numpy.ndarray`_: An array with shape :math:`(6,2)`, or
        :math:`(y,2)` if ``close`` is True, providing the x and y
        Cartesian vertices of the hexagon.
    """
    # Generate the vertices using hexgrid.  The diameter passed to
    # hexgrid() has to be adjusted to appropriately place the vertices
    # and match the input diameter.
    v = hexgrid(3, d/2/(numpy.cos(numpy.pi/6)**2 if incircle else numpy.cos(numpy.pi/6)))
    # Resort into counter-clockwise order and remove the point at 0,0
    # WARNING: This requires a hard-coded order for the grid returned by
    # hexgrid!
    v = numpy.vstack((v[:3], v[4:]))
    v = v[numpy.argsort(numpy.arctan2(v[:,1], v[:,0]))]
    # Append the first point to close the shape
    if close:
        v = numpy.vstack((v,v[0]))
    return set_orientation(v, orientation=orientation)


def rotate(coo, rot, deg=True):
    """
    Rotate a set of coordinates.

    Args:
        coo (`numpy.ndarray`_):
            Array of Cartesian coordinates to rotate.
        rot (:obj:`float`):
            Counter-clockwise rotation angle.
        deg (:obj:`bool`, optional):
            If True, provided rotation angle is in degrees.  If False,
            the angle is in radians.

    Returns:
        `numpy.ndarray`_: An array with the same shape as ``coo`` with
        the rotated coordinates.
    """
    theta = numpy.radians(rot) if deg else rot
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    r = numpy.array([[cos, sin], [-sin, cos]])
    return numpy.dot(coo, r)


# def test_cell():

#     v = hexvert(orientation='horizontal') #vertical')

#     w,h = pyplot.figaspect(1)
#     fig = pyplot.figure(figsize=(1.5*w,1.5*h))
#     ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#     ax.plot(v[:,0], v[:,1], color='k', lw=0.5)
#     ax.scatter(v[:,0], v[:,1], marker='x', color='k', s=50)
#     ax.set_xlim(-1, 1)
#     ax.set_ylim(-1, 1)
#     ax.add_patch(patches.Circle((0.,0.), radius=0.5, facecolor='none', edgecolor='C3', zorder=3))
#     ax.add_patch(patches.Circle((0.,0.), radius=0.5 * numpy.cos(numpy.pi/6.),
#                  facecolor='none', edgecolor='C0', zorder=3))

#     pyplot.show()


# def test_grid():

# #    v = hexgrid(3, 1., incircle=True, vertices=True)
# #    t = hexvert(d=2.5/numpy.cos(numpy.pi/6), incircle=True)

# #    v = hexgrid(3, 2.5/numpy.cos(numpy.pi/6), incircle=True, vertices=True, d_unit='grid')
# #    t = hexvert(d=2.5/numpy.cos(numpy.pi/6), incircle=True)

#     v = hexgrid(5, 4/numpy.cos(numpy.pi/6), incircle=True, vertices=True, d_unit='grid', full=True)
#     t = hexvert(d=4/numpy.cos(numpy.pi/6), incircle=True)

#     w,h = pyplot.figaspect(1)
#     fig = pyplot.figure(figsize=(1.5*w,1.5*h))
#     ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#     for c in v:
#         ax.plot(c[:,0], c[:,1], color='k', lw=0.5)
#         ax.scatter(c[:,0], c[:,1], marker='x', color='k', s=50)
#     ax.plot(t[:,0], t[:,1], color='C3', lw=0.5)
#     ax.scatter(t[:,0], t[:,1], marker='x', color='C3', s=50)
#     ax.set_xlim(-5, 5)
#     ax.set_ylim(-5, 5)
# #    ax.add_patch(patches.Circle((0.,0.), radius=0.5, facecolor='none', edgecolor='C3', zorder=3))
# #    ax.add_patch(patches.Circle((0.,0.), radius=0.5 * numpy.cos(numpy.pi/6.),
# #                 facecolor='none', edgecolor='C0', zorder=3))

#     pyplot.show()


