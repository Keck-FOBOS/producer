"""
Allocate apertures to targets.

Contains code originally written by 2021 Akamai intern, Brittany Ann
Ramos.

.. include:: ../include/links.rst
"""
from pathlib import Path
from IPython import embed 

import numpy


def random_targets(r, n=None, density=5., rng=None):
    r"""
    Draw a set of random x and y coordinates within a circle.

    Args:
        r (:obj:`float`):
            Radius of the circle.
        n (:obj:`int`, optional):
            The total number of points to draw.  If None, number drawn
            based on the density requested.
        density (:obj:`float`, optional):
            The average density of targets within the circle.  This is
            used to calculate the number of points to generate within
            the circle: ``n = int(numpy.ceil(density*numpy.pi*r**2))``.  Units
            must be appropriate match radius units.  
        rng (`numpy.random.Generator`_, optional):
            Random number generator to use.  If None, a new one is
            instantiated using `numpy.random.default_rng`_.

    Returns:
        :obj:`tuple`: Two vectors of length :math:`N_{\rm targ}` (the number of
        targets).  Cartesian x coordinates are in the first vector, y
        coordinates in the second.
    """
    # Calculate the number of points to match an expected density
    if n is None:
        n = int(numpy.ceil(density*numpy.pi*r**2))
    if rng is None:
        rng = numpy.random.default_rng()

    c = numpy.empty((0,2), dtype=float)
    overdraw = 1.5
    while c.shape[0] != n:
        # Draw more points than needed within the unit square until the correct
        # number is reached.
        # TODO: Probably a less brute-force way of doing this...
        c = rng.uniform(low=-1, high=1, size=(int(n*overdraw),2))
        # Find those within the r = 1
        indx = c[:,0]**2 + c[:,1]**2 < 1
        c = c[indx][:n]
        # Increase overdraw for next iteration
        overdraw *= 1.1
    
    return r*c[:,0], r*c[:,1]


def parse_targets(ifile, ra_c=1, dec_c=2, ap_c=None, default_ap=0):
    """
    Parse target coordinates and aperture types from an input file.

    Args:
        ifile (:obj:`str`):
            Columnated ascii file with the target coordinates.
        ra_c (:obj:`int`, optional):
            1-indexed column with the RA coordinates.  Assumed to be in decimal
            degrees.
        dec_c (:obj:`int`, optional):
            1-indexed column with the declination coordinates.  Assumed to be in
            decimal degrees.
        ap_c (:obj:`int`, optional):
            1-indexed column with the aperture type to assign to each target.
            If None, the type is not available in the input file and the
            ``default_types`` is used for all targets.  Apertures must be 0 for
            a single fiber or 1 for an IFU.
        default_ap (:obj:`int`, optional):
            If the aperture types are not provided in the file, this sets the
            type to assign to *all* apertures.  Apertures must be 0 for a single
            fiber or 1 for an IFU.

    Returns:
        :obj:`tuple`: Three numpy vectors with the coordinates and aperture type
        for each target.
    """
    # Check the input
    if default_ap not in [0, 1]:
        raise ValueError('Default aperture type must be 0 (single-fiber) or 1 (IFU).')
    # Instantiate the file path
    p = Path(ifile).resolve()
    # Check it exists
    if not p.exists():
        raise FileNotFoundError(f'{str(p)}')

    # Read the file
    db = numpy.genfromtxt(str(p), dtype=str).T
    # Check the requested columns exist
    if numpy.any(numpy.array([ra_c, dec_c, 1 if ap_c is None else ap_c]) > db.shape[0]):
        raise ValueError(f'{p.name} only contains {db.shape[0]} columns.  Check column requests.')
    # Collect the data and convert to the correct type
    ra = db[ra_c-1].astype(float)
    dec = db[dec_c-1].astype(float)
    ap = numpy.full(ra.size, default_ap, dtype=int) if ap_c is None else db[ap_c-1].astype(int)
    return ra, dec, ap






