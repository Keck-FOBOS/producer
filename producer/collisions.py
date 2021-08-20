"""
Functions used to handle aperture collisions.

.. include:: ../include/links.rst
"""

from IPython import embed 

import numpy as np
from matplotlib import pyplot, ticker

from sklearn.neighbors import KDTree

from .util import powerset


def remove_collisions_plot(x, y, keep, collision, groups=None):
    """
    Diagnostic plot for the results of :func:`remove_collisions`.

    Args:
        x (`numpy.ndarray`_):
            Cartesian x coordinates
        y (`numpy.ndarray`_):
            Cartesian x coordinates
        keep (`numpy.ndarray`_):
            Boolean array selecting the cartesian coordinates that *do
            not* collide.
        collision (:obj:`float`):
            The collision distance.
        groups (array-like, optional):
            The list of friends-of-friends groups constructed by
            :func:`remove_collisions`.
    """
    xlim = max(np.amax(np.absolute(x)), np.amax(np.absolute(y))) + 2*collision
    ylim = [-xlim, xlim]
    xlim = [-xlim, xlim]

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if groups is not None:
        for i, g in enumerate(groups):
            ax.text(np.mean(x[g]), np.mean(y[g]), str(i), ha='center', va='center',
                    color=f'C{i}')
            ax.scatter(x[g], y[g], marker='x', s=30, color=f'C{i}')

    for i, (_x, _y) in enumerate(zip(x, y)):
        ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='k', edgecolor='none',
                                    zorder=0, alpha=0.05))
        if keep[i]:
            ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='none',
                                        edgecolor='k', zorder=1))
    pyplot.show()


# TODO: Allow KDTree to be prebuilt?
# TODO: Enable use of BallTree with hersine distance metric
def remove_collisions(x, y, collision, maxgrp=18, plot=False):
    r"""
    Remove "collisions" in a set of Cartesian coordinates.

    Two coordinates collide if they are closer than the provided minimum
    distance.

    The algorithm builds friends-of-friends groupings with a linking
    length set to the collision distance, and then determines the
    maximum number of group members that are all separated by at least
    the collision distance.

    Iteration through all combinations of group members quickly becomes
    intractable for large groups; for a group of size :math:`n`, the
    total number of combinations is :math:2^n`.  The ``maxgroup``
    parameter sets the maximum size allowed for iterating through all
    possible combinations.  Groups larger than this value are
    recursively "thinned" until its members are no greater than
    ``maxgrp``.  This thinning process can be very expensive for large
    groups.  For fields that are dense relative to the collision
    distance, first consider iteratively running this function by
    starting with a relatively small collision distance and gradually
    increasing it to the desired value.

    Args:
        x (`numpy.ndarray`_):
            Cartesian x coordinates
        y (`numpy.ndarray`_):
            Cartesian y coordinates
        collision (:obj:`float`):
            The minimum distance to impose between coordinates.
        maxgrp (:obj:`int`, optional):
            The maximum size of a friends-of-friends group for which to
            iterate through all possible combinations to find the
            largest set of uncollided coordinates.  See the function
            description.
        plot (:obj:`bool`, optional):
            Show a plot of the result.

    Returns:
        :obj:`tuple`: The function returns the `numpy.ndarray`_ (array
        of integer arrays) with the friends-of-friends groups (each
        group is set by an integer array providing the indices of the
        coordinate members) and a `numpy.ndarray`_ with the indices of
        the uncollided coordinates.
    """
    # Build a KDTree of the input coordinates
    kdtree = KDTree(np.column_stack((x,y)))

    # Use the tree to find the indices of points within the collision
    # distance of every other point
    groups = kdtree.query_radius(np.column_stack((x,y)), collision) #.tolist()
    grpn = np.array([g.size for g in groups])

    # This shortcut circumvents having to deal with the collisions below
    # by checking for the case when none of the coordinates actually
    # collide.
    n = x.size
    if np.all(grpn == 1):
        return groups, np.ones(n, dtype=bool)

    uniq = np.unique(np.concatenate(groups))

    merged = np.zeros(len(groups), dtype=bool)
    for i in range(groups.shape[0]):
        if len(groups[i]) == 1 or merged[i]:
            continue
        for j in range(i+1, groups.shape[0]):
            if len(groups[j]) == 1 or merged[j]:
                continue
            if len(np.intersect1d(groups[i], groups[j])) > 0:
                groups[i] = np.union1d(groups[i], groups[j])
                merged[j] = True

    groups = np.delete(groups, merged)

    if uniq.size != np.unique(np.concatenate(groups)).size:
        raise ValueError('CODING ERROR: Something wrong in group consolidation.')

    # Now iterate through each group:
    #   - Recursively thin large groups.  The recursion can lead to very
    #     long execution times for groups where the number of members is
    #     significantly larger than maxgrp.
    #   - Points in single-member groups are kept.
    #   - Exclude the second point in a pair groups.
    #   - For groups with three or more, start with the largest set (the
    #     length of the group - 1) and iterate through all combinations
    #     of subsets to find the first one (the one with the maximum
    #     number of members) without any collided coordinates.
    keep = np.zeros(n, dtype=bool)
    for i, _g in enumerate(groups):
        nmaxiter = 0
        g = _g.copy()
        # Deal with small groups
        if len(g) < 3:
            keep[g[0]] = True
            continue
        # Thin large groups to make the number of combinations tested
        # below tractable.
        while len(g) > maxgrp:
            # The change to the collision radius in this procedure is
            # completely ad hoc
            _groups, _keep = remove_collisions(x[_g], y[_g], (0.9 + 0.01*nmaxiter)*collision,
                                               maxgrp=maxgrp)
            g = _g[_keep]
            nmaxiter += 1
        # Starting with the largest subsets, iterate through all
        # combinations of group members until a set is found where all
        # the distances are larger than the collision distance.
        indx_combinations = powerset(np.arange(len(g)), reverse=True)
        for indx in indx_combinations:
            if len(indx) in [0, len(g)]:
                continue
            _indx = np.array(indx)
            if _indx.size < 2:
                # Only one point left, use it
                keep[g[_indx]] = True
                break
            # Calculate the upper triangle of the distance matrix (i.e.,
            # the distance from each point to every other point)
            di, dj = np.triu_indices(_indx.size, k=1)
            dist = (x[g[_indx[di]]] - x[g[_indx[dj]]])**2 + (y[g[_indx[di]]] - y[g[_indx[dj]]])**2

            # If the distances are all greater than the collision
            # distance, we're done and we can keep all of the
            # coordinates
            if np.all(dist > collision**2):
                keep[g[_indx]] = True
                break

        if nmaxiter > 0:
            # Recover well separated points from large groups
            indx = keep[_g]
            nindx = np.logical_not(indx)
            dist = (x[_g[nindx],None] - x[None,_g[indx]])**2 \
                        + (y[_g[nindx],None] - y[None,_g[indx]])**2
            keep[_g[nindx]] = np.all(dist > collision**2, axis=1)

    if plot:
        # Plot the groups and the thinned dataset
        remove_collisions_plot(x, y, keep, collision, groups=groups)

    return groups, keep



