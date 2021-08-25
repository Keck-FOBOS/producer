"""
Generic methods for distributing starbug modules in a focal plane.

.. include:: ../include/links.rst
"""

import time

from IPython import embed

import numpy

from .. import hexgrid 

def module_layout(platescale, pitch, n_hex_module, n_hex_focal, module_pitch=None, max_dist=None):
    r"""
    Construct a hexagonal grid of Starbug modules for a focal plane.

    Args:
        platescale (:obj:`float`):
            Telescope plate scale in mm / arcsec.
        pitch (:obj:`float`):
            The incircle of each Starbug home hexagon in mm.  I.e., the pitch
            between each Starbug home position.
        n_hex_module (:obj:`int`):
            The number of Starbugs along the long axis of the module.  This sets
            the number of Starbugs per module.  Must be an odd number.  E.g., 3
            is 7, 5 is 19, etc.
        n_hex_focal (:obj:`int`):
            The number of modules along the long axis of the focal-plane.  This
            sets the number of modules in the focal-plane grid.  Must be an odd
            number. E.g., 3 is 7, 5 is 19, etc.
        module_pitch (:obj:`float`):
            The pitch in mm between the centers of each module.  If None, this
            is the minimum required to meet the needs of the number of modules
            and the Starbug pitch.  See
            :func:`~producer.hexgrid.hexgrid_grid_convert`.
        max_dist (:obj:`float`):
            The maximum distance in arcmin from the focal plane center to the
            module center allowed for including a module in the final layout.
            If None, no limit is imposed.

    Returns:
        :obj:`tuple`: Two 2D `numpy.ndarray`_ objects with (1) the coordinates
        for the Starbugs in a module relative to the module center, with shape
        :math:`(N_{\rm bug},2)`, and (2) the coordinates for the modules
        relative to the focal-plane center, with shape :math:`(N_{\rm mod},2)`.
    """
    # Get the Starbug center grid for an individual module (centered at
    # 0.0)
    starbug_grid = hexgrid.hexgrid(n_hex_module, pitch/platescale, incircle=True,
                                   orientation='vertical')
    
    # Get the pitch between modules
    if module_pitch is None:
        module_pitch = hexgrid.hexgrid_grid_convert(n_hex_module, pitch, incircle=True,
                                                    d_unit='cell', full=False)

    # Module grid
    module_grid = hexgrid.hexgrid(n_hex_focal, module_pitch/platescale, incircle=True)
    if max_dist is not None:
        module_dist = numpy.sqrt(module_grid[:,0]**2 + module_grid[:,1]**2)
        module_dist /= 60.
        keep_module = numpy.ones(module_dist.shape, dtype=bool) if max_dist is None \
                            else module_dist < max_dist
        module_grid = module_grid[keep_module]

    return starbug_grid, module_grid


