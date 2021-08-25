
.. include:: include/links.rst

.. _knownissues:

Known issues
============

- ``fobos_makeplan`` treats all tiles independently, meaning that single
  objects in tile overlap regions will be allocated to multiple tiles.  Need a
  better way to optimize this.

- ``fobos_makeplan`` currently does not allow for different exposure times per
  object. 

- As far as the plots from ``fobos_showcfg`` correctly illustrate this, it
  appears collisions are not completely removed during the aperture assignments.

- Collision radius currently set to 10 arcsec, but current design suggests this
  should be larger. 

- ``fobos_makeplan`` does not reposition unallocated apertures to avoid
  collisions with allocated fibers. 

- ``fobos_makeplan`` does not allocate sky apertures. 



