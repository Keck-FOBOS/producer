"""
Functions used to tile a set of on-sky area with FOBOS pointings.

.. include:: ../include/links.rst
"""
import warnings

from IPython import embed

import numpy

from scipy import spatial

from sklearn.neighbors import KDTree

from astropy import units
from astropy import coordinates

from . import hexgrid
from . import util


# TODO: 
#   - Enable a 1/2 FOV offset; i.e., offset the "center" tile so that it isn't
#     pointing at the ``grid_center``
#   - Enable a scaling of the offset between pointing centers to increase
#     overlap between fields.
#   - Deal with tangent plane distortions for large on-sky distributions of
#     targets.  I.e., when do distortions become significant?  Or is there a
#     better coordinate projection to use?  Break the data up into chunks (see
#     pydl spheregroup)?
#   - Use BallTree instead of KDTree?
def uniform_tiles(ra, dec, grid_center=None, min_dens=None, min_n=None, fov=20./60.):
    """
    Construct a uniform grid of FOBOS pointings to observe a set of targets.

    Args:
        ra (`numpy.ndarray`_):
            Right ascension of objects to observe in decimal degrees.
        dec (`numpy.ndarray`_):
            Declination of objects to observe in decimal degrees.
        grid_center (:obj:`tuple`, optional):
            Tuple with a starting point (ra, dec) for the grid center in decimal
            degrees.  If None, starting point set at mean of provided
            coordinates.
        min_dens (:obj:`float`, optional):
            Minimum *density* of targets within a given pointing to allow before
            removing it from the grid.  Density calculations use the area that
            *exclude* the overlap between tiles.  If None and ``min_n`` is None,
            all fields with *any* targets are returned.  Ignored if ``min_n`` is
            provided.
        min_n (:obj:`int`, optional):
            Minimum *number* of targets within a given pointing to allow before
            removing it from the returned grid.  If None, set to 0.  If
            ``min_dens`` is provided, the area is used to set the minimum
            number.
        fov (:obj:`float`, optional):
            The *diameter* of the field of view in decimal degrees.  The uniform
            grid is constructed using hexagons with the circular field-of-view
            as the circumcircle of each hexagonal grid cell.  See
            :func:`~producer.hexgrid.hexgrid`.

    Returns:
        `numpy.ndarray`_: A 2D `numpy.ndarray`_ with field centers vectors with
        the center coordinates in decimal degrees for all grid pointings.
    """
    if min_dens is not None and min_n is not None:
        warnings.warn('Provided both min_dens and min_n.  Ignoring min_dens.')
    if min_n is None:
        min_n = 0 if min_dens is None else min_dens * hexgrid.hexarea(d=fov)

    if grid_center is None:
        grid_center = (numpy.mean(ra), numpy.mean(dec))

    # Project the coordinates to the tangent plane
    ntarg = ra.size
    center = coordinates.SkyCoord(grid_center[0]*units.deg, grid_center[1]*units.deg, frame='icrs')
    targets = coordinates.SkyCoord(ra*units.deg, dec*units.deg, frame='icrs')
    targ_offsets = targets.transform_to(coordinates.SkyOffsetFrame(origin=center))
    x = targ_offsets.lon.value
    y = targ_offsets.lat.value
    coo = numpy.column_stack((x,y))

    # Get the length of the long axis for the grid
    hull = spatial.ConvexHull(coo).vertices
    i, j = map(lambda x: hull[x], numpy.triu_indices(len(hull), k=1))
    sep = numpy.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2)
    ii = numpy.argmax(sep)
    # Set the grid width so that its short axis is the same length as the
    # longest span of the coordinates
    width = hexgrid.hexgrid_circle_convert(sep[ii], incircle=True)
    rot = numpy.arctan2(y[i[ii]] - y[j[ii]], x[i[ii]] - x[j[ii]])

    # Get the number of grid cells along the long axis of the target distribution
    n = int(numpy.ceil(width/fov)) - 2
    if n % 2 == 0:
        n += 1

    # Even with the above setting of the long axis of the grid to match the long
    # axis of the target distribution, may still miss targets.  This iteratively
    # constructs a uniform grid until all targets are covered by a grid cell.
    # Uses a KD-tree to speed up searching for which targets are in each
    # pointing.
    kdtree = KDTree(coo)
    in_grid = numpy.zeros(ntarg, dtype=bool)
    niter = 0
    while not numpy.all(in_grid):
        n += 2
        grid = hexgrid.hexgrid(n, fov, orientation=numpy.degrees(rot))
        groups = kdtree.query_radius(grid, fov/2.)
        # TODO: Could just check that the hull vertices are found...
        in_grid = numpy.isin(numpy.arange(ntarg), numpy.unique(numpy.concatenate(groups)))
        niter += 1

    # Only keep the grid points with sufficient targets
    keep_grid = numpy.array([g.size > min_n for g in groups])

    # Revert back to on-sky coordinates and return the grid centers
    grid_coo = coordinates.SkyCoord(grid[keep_grid,0]*units.deg, grid[keep_grid,1]*units.deg,
                                    frame=center.skyoffset_frame())
    grid_coo = coordinates.SkyCoord(grid_coo.frame.transform_to(coordinates.ICRS()))
    return numpy.column_stack((grid_coo.ra.value, grid_coo.dec.value))








