"""
Functions used to tile a set of on-sky area with FOBOS pointings.

.. include:: ../include/links.rst
"""
import warnings

from IPython import embed

import numpy

from scipy import spatial

from matplotlib import pyplot, patches

from sklearn.neighbors import KDTree

from astropy import units
from astropy import coordinates

from . import hexgrid
from . import util
from .astrometry import focal_plane_offsets


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
def uniform_tiles(ra, dec, grid_center=None, min_dens=None, min_n=None, fov=20./60.,
                  half_offset=False):
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
        half_offset (:obj:`bool`, optional):
            Offset the center by half the field of view.  If not grid center is
            provided, offset the default grid center (the mean of the provided
            target coordinates) by half of the field-of-view.  Ignored if
            ``grid_center`` is provided.

    Returns:
        `numpy.ndarray`_: A 2D `numpy.ndarray`_ with field centers vectors with
        the center coordinates in decimal degrees for all grid pointings.
    """
    if min_dens is not None and min_n is not None:
        warnings.warn('Provided both min_dens and min_n.  Ignoring min_dens.')
    if min_n is None:
        min_n = 0 if min_dens is None else min_dens * hexgrid.hexarea(d=fov)

    if grid_center is None:
        # TODO: Need to consider spherical geometry here!
        grid_center = (numpy.mean(ra) + (fov/2. if half_offset else 0.), numpy.mean(dec))

    # Project the coordinates to the tangent plane
    ntarg = ra.size
    coo = numpy.column_stack(focal_plane_offsets(ra, dec, grid_center))

    # Get the length of the long axis for the grid
    hull = spatial.ConvexHull(coo).vertices
    i, j = map(lambda x: hull[x], numpy.triu_indices(len(hull), k=1))
    sep = numpy.sqrt((coo[i,0] - coo[j,0])**2 + (coo[i,1] - coo[j,1])**2)
    ii = numpy.argmax(sep)
    # Set the grid width so that its short axis is the same length as the
    # longest span of the coordinates
    width = hexgrid.hexgrid_circle_convert(sep[ii], incircle=True)
    rot = numpy.arctan2(coo[i[ii],1] - coo[j[ii],1], coo[i[ii],0] - coo[j[ii],0])

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
    return numpy.column_stack(focal_plane_offsets(grid[keep_grid,0], grid[keep_grid,1],
                                                  grid_center, revert=True))


# TODO: Add projection type...
def show_tiles(tile_coo, ra=None, dec=None, fov=20./60., return_ax=False):
    """
    Args:
        tile_coo (`numpy.ndarray`_):
            Tile coordinates in RA (first column) and DEC (second column).
        ra (`numpy.ndarray`_, optional):
            Right ascension of objects to observe in decimal degrees.
        dec (`numpy.ndarray`_, optional):
            Declination of objects to observe in decimal degrees.
        fov (:obj:`float`, optional):
            The *diameter* of the field of view in decimal degrees.  The uniform
            grid is constructed using hexagons with the circular field-of-view
            as the circumcircle of each hexagonal grid cell.  See
            :func:`~producer.hexgrid.hexgrid`.
        return_ax (:obj:`bool`, optional):
            Instead of showing the plot, return the axis instance.

    Returns:
        Axis: Axis instance.  If ``return_ax`` is False, this is returned as
        None.
    """
    d = numpy.amax(numpy.amax(tile_coo, axis=0) - numpy.amin(tile_coo, axis=0)) + 2*fov
    xlim, ylim = numpy.mean(tile_coo, axis=0)[:,None] + numpy.array([-d/2, d/2])[None,:]

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(2.*w,2.*h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.minorticks_on()
    ax.grid(True, which='major', color='0.9', zorder=0, linestyle='-')
    ax.tick_params(which='major', direction='in', length=8, top=True, right=True)
    ax.tick_params(which='minor', direction='in', length=4, top=True, right=True)
    if ra is not None and dec is not None:
        ax.scatter(ra, dec, marker='.', s=30, lw=0, color='k', zorder=3, alpha=0.1)
    for i, c in enumerate(tile_coo):
        ax.add_patch(patches.Circle((c[0],c[1]), radius=fov/2, facecolor='C0', edgecolor='C0',
                                    zorder=5, alpha=0.1))
        ax.text(c[0], c[1], f'{i+1}', ha='center', va='center', color='C0', fontsize=16)
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('DEC [deg]')
    if return_ax:
        return ax
    
    pyplot.show()
    return None





