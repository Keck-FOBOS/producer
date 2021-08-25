"""
Miscellaneous astrometry utilities.

.. include:: ../include/links.rst
"""

from astropy import units
from astropy import coordinates

def focal_plane_offsets(x, y, center, revert=False):
    """
    Convert between on-sky coordinates and focal-plane offset coordinates.

    Assumes ICRS frame.

    Args:
        x (`numpy.ndarray`_):
            Right ascension or longitude offset in decimal degrees.
        y (`numpy.ndarray`_):
            Declination or latitude offset in decimal degrees.
        center (:obj:`tuple`):
            Two-tuple with the center coordinate, right ascension and
            declination, for the offset computations in decimal degrees.
        revert (:obj:`bool`, optional):
            If True, the computation assumes the input coordinates (``x`` and
            ``y``) are offset coordinates and the function reverts them to RA
            and DEC given the provided center.

    Returns:
        :obj:`tuple`: Two `numpy.ndarray`_ objects with the converted
        coordinates in decimal degrees.  Coordinates are RA and DEC if
        ``revert`` is true; otherwise, they are longitude and latitude offsets
        relative to the center.
    """
    center = coordinates.SkyCoord(center[0]*units.deg, center[1]*units.deg, frame='icrs')

    if revert:
        coo = coordinates.SkyCoord(x*units.deg, y*units.deg, frame=center.skyoffset_frame())
        coo = coordinates.SkyCoord(coo.frame.transform_to(coordinates.ICRS()))
        return coo.ra.value, coo.dec.value

    coo = coordinates.SkyCoord(x*units.deg, y*units.deg, frame='icrs')
    offsets = coo.transform_to(coordinates.SkyOffsetFrame(origin=center))
    return offsets.lon.value, offsets.lat.value


