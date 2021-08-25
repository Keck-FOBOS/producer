"""
Allocate apertures to targets.

Contains code originally written by 2021 Akamai intern, Brittany Ann
Ramos.

.. include:: ../include/links.rst
"""

from IPython import embed 

import numpy
from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot, patches

from sklearn.neighbors import KDTree

from .collisions import remove_collisions


def assign_apertures_plot(objx, objy, apx, apy, assigned_obj, assigned_ap, collision=1/6.,
                          ignore_obj=None, ignore_ap=None, zoom=0.75):
    """
    Show a plot of the aperture assignments.  I.e.:

    .. code-block:: python

        assigned_obj, assigned_ap = assign_apertures(objx, objy, apx, apy)
        assign_apertures_plot(objx, objy, apx, apy, assigned_obj, assigned_ap)

    Args:
        objx (`numpy.ndarray`_):
            Cartesian x coordinate in tangent plane projection of object
            coordinates relative to the pointing center.  Shape must be
            1D and match ``objy``.
        objy (`numpy.ndarray`_):
            Cartesian y coordinate in tangent plane projection of object
            coordinates relative to the pointing center. Shape must be
            1D and match ``objx``.
        apx (`numpy.ndarray`_):
            Cartesian x coordinate in tangent plane projection of the
            "home" aperture coordinates relative to the pointing center.
            Shape must be 1D and match ``apy``.
        apy (`numpy.ndarray`_):
            Cartesian y coordinate in tangent plane projection of the
            "home" aperture coordinates relative to the pointing center.
            Shape must be 1D and match ``apx``.
        assigned_obj (`numpy.ndarray`_):
            Indices of the objects with assigned apertures.  This is the first
            object returned by :func:`~producer.allocate.assign_apertures`.
        assigned_ap (`numpy.ndarray`_):
            Indices of the apertures assigned to objects.  This is the second
            object returned by :func:`~producer.allocate.assign_apertures`.
        collision (:obj:`float`, optional):
            Aperture-to-aperture collision radius.
        ignore_obj (`numpy.ndarray`_, optional):
            An integer array with the indices of objects that were *ignored*
            during assignment.  If None, all objects were included in the
            assignment attempt.
        ignore_ap (`numpy.ndarray`_, optional):
            An integer array with the indices of apertures that were *ignored*
            during assignment.  If None, all apertures were included in the
            assignment attempt.
        zoom (:obj:`float`, optional):
            Factor to zoom in default limits of the plot relative to the nominal
            set by the extent of the allocated targets.  For example, ``zoom=2``
            means to *decrease* the plot range by a factor of two; to zoom out,
            use a zoom value that is less than 1.
    """
    # Get the required width in each direction
    Dx = (numpy.amax(objx[assigned_obj]) - numpy.amin(objx[assigned_obj]))/zoom
    Dy = (numpy.amax(objy[assigned_obj]) - numpy.amin(objy[assigned_obj]))/zoom
    # Force it to be square
    D = max(Dx, Dy)
    
    xlim = numpy.mean(objx[assigned_obj]) + numpy.array([-D/2, D/2])
    ylim = numpy.mean(objy[assigned_obj]) + numpy.array([-D/2, D/2])

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    ax.minorticks_on()
    ax.grid(True, which='major', color='0.9', zorder=0, linestyle='-')
    ax.tick_params(which='major', direction='in', length=8, top=True, right=True)
    ax.tick_params(which='minor', direction='in', length=4, top=True, right=True)

    for i, (_x, _y) in enumerate(zip(objx[assigned_obj], objy[assigned_obj])):
        kwargs = {} # if i > 0 else {'label': 'Assigned Aper.'}
        ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='C0',
                                    edgecolor='none', zorder=0, alpha=0.2, **kwargs))
    if len(assigned_ap) != apx.size:
        unassigned = numpy.setdiff1d(numpy.arange(apx.size), assigned_ap)
        if ignore_ap is not None:
            unassigned = numpy.setdiff1d(unassigned, ignore_ap)
        for i, (_x, _y) in enumerate(zip(apx[unassigned], apy[unassigned])):
            kwargs = {} # if i > 0 else {'label': 'Unassigned Aper.'}
            ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='C1',
                                        edgecolor='none', zorder=0, alpha=0.2, **kwargs))

    ax.scatter(objx[assigned_obj], objy[assigned_obj], marker='.', color='C0', s=30, lw=0,
               label='Alloc. Target')
    if len(assigned_obj) != objx.size:
        unassigned = numpy.setdiff1d(numpy.arange(objx.size), assigned_obj)
        if ignore_obj is not None:
            unassigned = numpy.setdiff1d(unassigned, ignore_obj)
        ax.scatter(objx[unassigned], objy[unassigned], marker='.', color='C1', s=30, lw=0,
                   label='Unalloc. Target')

    # FOBOS field-of-view
    ax.add_patch(patches.Circle((0.,0.), radius=10., facecolor='none', edgecolor='C3',
                                zorder=4))

    ax.set_xlabel(r'$\xi$ [arcmin]')
    ax.set_ylabel(r'$\eta$ [arcmin]')
    # TODO: The legend really slows down the plot...
#    ax.legend()
    
    pyplot.show()


def assign_apertures(objx, objy, apx, apy, collision=1/6., patrol=2.3, ignore_obj=None,
                     ignore_ap=None, allocated_obj=None, allocated_ap=None):
    """
    Assign free apertures to objects.

    Coordinate units must be the same for all arguments.  I.e., if the
    object and aperture coordinates are in arcmin, the collision and
    patrol radii must be also.

    To avoid collisions with previously allocated objects, use
    ``allocated_obj``.  Providing ``allocated_ap`` just removes these
    apertures from consideration; i.e., the information of where the
    apertures are located in the field is associated with
    ``allocated_obj``, not ``allocated_ap``.  Providing one does not
    mean you have to provide the other.

    Objects and apertures can also be ignored using ``ignore_obj`` and
    ``ignore_fib``.  Ignored objects/apertures are not considered during
    assignment.  These arguments are largely provided as book-keeping
    convenience; i.e., the returned indices are based on the full input
    arrays instead of a preselected subset passed to this function.

    Args:
        objx (`numpy.ndarray`_):
            Cartesian x coordinate in tangent plane projection of object
            coordinates relative to the pointing center.  Shape must be
            1D and match ``objy``.
        objy (`numpy.ndarray`_):
            Cartesian y coordinate in tangent plane projection of object
            coordinates relative to the pointing center. Shape must be
            1D and match ``objx``.
        apx (`numpy.ndarray`_):
            Cartesian x coordinate in tangent plane projection of the
            "home" aperture coordinates relative to the pointing center.
            Shape must be 1D and match ``apy``.
        apy (`numpy.ndarray`_):
            Cartesian y coordinate in tangent plane projection of the
            "home" aperture coordinates relative to the pointing center.
            Shape must be 1D and match ``apx``.
        collision (:obj:`float`, optional):
            Aperture-to-aperture collision radius.
        patrol (:obj:`float`, optional):
            Aperture patrol radius.
        ignore_obj (`numpy.ndarray`_, optional):
            An integer array with the indices of objects to *ignore*
            during assignment.  If None, all objects are included in
            assignment attempt.
        ignore_ap (`numpy.ndarray`_, optional):
            An integer array with the indices of apertures to *ignore*
            during assignment.  If None, all apertures are included in
            assignment attempt.
        allocated_obj (`numpy.ndarray`_, optional):
            An integer array with the indices of objects that have
            already been allocated to targets in the field, as needed
            to avoid collisions.
        allocated_ap (`numpy.ndarray`_, optional):
            An integer array with the indices of apertures that have
            already been allocated to targets in the field.  This is
            only used in book-keeping; i.e., allocated apertures cannot
            be reassigned and the returned assignment indices will match
            the indices of the provided ``apx`` and ``apy`` vectors.

    Returns:
        :obj:`tuple`: Two `numpy.ndarray`_ objects with the indices of
        the objects that have been assigned apertures, and the indices
        of the associated apertures.
    """

    # Check input
    if objx.ndim > 1:
        raise ValueError('Input object coordinates must be in 1D arrays.')
    if objx.shape != objy.shape:
        raise ValueError('Shape of input object coordinate arrays do not match.')
    if apx.ndim > 1:
        raise ValueError('Input aperture coordinates must be in 1D arrays.')
    if apy.shape != apy.shape:
        raise ValueError('Shape of input aperture coordinate arrays do not match.')
   
    # Number of objects and apertures
    nobj = objx.size
    nap = apx.size

    # Set the objects available for assignment and those already
    # allocated.
    _ignore = numpy.empty(0, dtype=int) if ignore_obj is None else ignore_obj.copy()
    if allocated_obj is None:
        prealloc_x = None
        prealloc_y = None
    else:
        _ignore = numpy.unique(numpy.append(_ignore, allocated_obj))
        prealloc_x = objx[allocated_obj]
        prealloc_y = objy[allocated_obj]
    _objx = objx.copy() if len(_ignore) == 0 else numpy.delete(objx, _ignore)
    _objy = objy.copy() if len(_ignore) == 0 else numpy.delete(objy, _ignore)
    obj_id = numpy.arange(nobj) if len(_ignore) == 0 else numpy.delete(numpy.arange(nobj), _ignore)

    # Set the apertures available for assignment
    _ignore = numpy.empty(0, dtype=int) if ignore_ap is None else ignore_ap.copy()
    if allocated_ap is not None:
        _ignore = numpy.unique(numpy.append(_ignore, allocated_ap))
    _apx = apx.copy() if len(_ignore) == 0 else numpy.delete(apx, _ignore)
    _apy = apy.copy() if len(_ignore) == 0 else numpy.delete(apy, _ignore)
    ap_id = numpy.arange(nap) if len(_ignore) == 0 else numpy.delete(numpy.arange(nap), _ignore)

    # Track the number of iterations required to avoid collisions
    niter = 0
  
    # Arrays to keep track of assigned apertures and objects
    assigned_ap = numpy.empty(0, dtype=int)
    assigned_obj = numpy.empty(0, dtype=int)

    # Iteratively assign apertures to objects until either no objects
    # are left or no apertures can be assigned
    while _objx.size > 0 and _apx.size > 0:

        # Build a KDTree and query the nearest nquery number of objects
        # for each aperture
        kdtree = KDTree(numpy.column_stack((_objx, _objy)))
        nquery = min(20, _objx.size)
        d, indx = kdtree.query(numpy.column_stack((_apx, _apy)), k=nquery)

        # Setup the distance array for these objects and assign each
        # aperture to a unique object.  Distances for object-aperture
        # pairs not returned by the KD-Tree query are set to a large
        # number, effectively forcing them to be ignored by the
        # linear-sum assignment function.
        api = numpy.tile(numpy.arange(_apx.size), (nquery,1)).T
        dist = numpy.full((_apx.size, _objx.size), 10*numpy.amax(d), dtype=float)
        dist[api.flat,indx.flat] = d.flat
        _a_ap, _a_obj = linear_sum_assignment(dist)

        # Select the assignments that are within the patrol radius of
        # the aperture and ignore the rest
        indx = dist[_a_ap,_a_obj] < patrol
        _a_ap = _a_ap[indx]
        _a_obj = _a_obj[indx]
        if _a_ap.size == 0:
            # No apertures can be assigned within the patrol radius, so
            # we're done!
            break

        # Thin the assigned coordinates of collided objects.  This needs
        # to account for all previously assigned apertures.
        if len(assigned_ap) == 0 and prealloc_x is None:
            not_collided = remove_collisions(_objx[_a_obj], _objy[_a_obj], collision)[1]
        else:
            if len(assigned_ap) == 0:
                _pre_x = prealloc_x
                _pre_y = prealloc_y
            elif prealloc_x is None:
                _pre_x = objx[assigned_obj]
                _pre_y = objy[assigned_obj]
            else:
                _pre_x = numpy.append(prealloc_x, objx[assigned_obj])
                _pre_y = numpy.append(prealloc_y, objy[assigned_obj])
            npre = _pre_x.size
            not_collided = remove_collisions(numpy.append(_pre_x, _objx[_a_obj]),
                                             numpy.append(_pre_y, _objy[_a_obj]),
                                             collision)[1][npre:]

        # Add uncollided assignments to list of assigned apertures and
        # objects
        assigned_ap = numpy.append(assigned_ap, ap_id[_a_ap[not_collided]])
        assigned_obj = numpy.append(assigned_obj, obj_id[_a_obj[not_collided]])

        # Remove the successfully assigned apertures from those available
        # to be reassigned
        _apx = numpy.delete(_apx, _a_ap[not_collided])
        _apy = numpy.delete(_apy, _a_ap[not_collided])
        ap_id = numpy.delete(ap_id, _a_ap[not_collided])

        # Remove the objects selected in the assignment, regardless of
        # whether or not they collide.  I.e., objects that collide with
        # currently assigned apertures will continue to collide!
        _objx = numpy.delete(_objx, _a_obj)
        _objy = numpy.delete(_objy, _a_obj)
        obj_id = numpy.delete(obj_id, _a_obj)

        niter += 1

    # TODO: Print a report?

    return assigned_obj, assigned_ap



