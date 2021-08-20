"""
Allocate apertures to targets.

Contains code originally written by 2021 Akamai intern, Brittany Ann
Ramos.

.. include:: ../include/links.rst
"""

import time
from itertools import chain, combinations

from IPython import embed 

import numpy as np
from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot, patches, ticker

from sklearn.neighbors import KDTree

from astropy import units
from astropy.coordinates import SkyCoord, SkyOffsetFrame

from .collisions import remove_collisions


def assign_apertures(objx, objy, apx, apy, collision=1/6., patrol=2.3, ignore_obj=None,
                     ignore_ap=None, allocated_obj=None, allocated_ap=None):
    """
    Assign apertures to objects.

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
    _ignore = np.empty(0, dtype=int) if ignore_obj is None else ignore_obj.copy()
    if allocated_obj is None:
        prealloc_x = None
        prealloc_y = None
    else:
        _ignore = np.unique(np.append(_ignore, allocated_obj))
        prealloc_x = objx[allocated_obj]
        prealloc_y = objy[allocated_obj]
    _objx = objx.copy() if len(_ignore) == 0 else np.delete(objx, _ignore)
    _objy = objy.copy() if len(_ignore) == 0 else np.delete(objy, _ignore)
    obj_id = np.arange(nobj) if len(_ignore) == 0 else np.delete(np.arange(nobj), _ignore)

    # Set the apertures available for assignment
    _ignore = np.empty(0, dtype=int) if ignore_ap is None else ignore_ap.copy()
    if allocated_ap is not None:
        _ignore = np.unique(np.append(_ignore, allocated_ap))
    _apx = apx.copy() if len(_ignore) == 0 else np.delete(apx, _ignore)
    _apy = apy.copy() if len(_ignore) == 0 else np.delete(apy, _ignore)
    ap_id = np.arange(nap) if len(_ignore) == 0 else np.delete(np.arange(nap), _ignore)

    # Track the number of iterations required to avoid collisions
    niter = 0
  
    # Arrays to keep track of assigned apertures and objects
    assigned_ap = np.empty(0, dtype=int)
    assigned_obj = np.empty(0, dtype=int)

    # Iteratively assign apertures to objects until either no objects
    # are left or no apertures can be assigned
    while _objx.size > 0 and _apx.size > 0:

        # Build a KDTree and query the nearest nquery number of objects
        # for each aperture
        kdtree = KDTree(np.column_stack((_objx, _objy)))
        nquery = min(20, _objx.size)
        d, indx = kdtree.query(np.column_stack((_apx, _apy)), k=nquery)

        # Setup the distance array for these objects and assign each
        # aperture to a unique object.  Distances for object-aperture
        # pairs not returned by the KD-Tree query are set to a large
        # number, effectively forcing them to be ignored by the
        # linear-sum assignment function.
        api = np.tile(np.arange(_apx.size), (nquery,1)).T
        dist = np.full((_apx.size, _objx.size), 10*np.amax(d), dtype=float)
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
                _pre_x = np.append(prealloc_x, objx[assigned_obj])
                _pre_y = np.append(prealloc_y, objy[assigned_obj])
            npre = _pre_x.size
            not_collided = remove_collisions(np.append(_pre_x, _objx[_a_obj]),
                                             np.append(_pre_y, _objy[_a_obj]),
                                             collision)[1][npre:]

        # Add uncollided assignments to list of assigned apertures and
        # objects
        assigned_ap = np.append(assigned_ap, ap_id[_a_ap[not_collided]])
        assigned_obj = np.append(assigned_obj, obj_id[_a_obj[not_collided]])

        # Remove the successfully assigned apertures from those available
        # to be reassigned
        _apx = np.delete(_apx, _a_ap[not_collided])
        _apy = np.delete(_apy, _a_ap[not_collided])
        ap_id = np.delete(ap_id, _a_ap[not_collided])

        # Remove the objects selected in the assignment, regardless of
        # whether or not they collide.  I.e., objects that collide with
        # currently assigned apertures will continue to collide!
        _objx = np.delete(_objx, _a_obj)
        _objy = np.delete(_objy, _a_obj)
        obj_id = np.delete(obj_id, _a_obj)

        niter += 1

    # TODO: Print a report?

    return assigned_obj, assigned_ap


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
            the circle: ``n = int(np.ceil(density*np.pi*r**2))``.  Units
            must be appropriate match radius units.  
        rng (`numpy.random.Generator`_, optional):
            Random number generator to use.  If None, a new one is
            instantiated using `numpy.random.default_rng`_.

    Returns:
        `numpy.ndarray`_: Array of shape :math:`(N_{\rm targ}, 2)`,
        where :math:`N_{\rm targ}` is the number of targets.  Cartesian
        x coordinates are in the first column, y coordinates in the
        second.
    """
    # Calculate the number of points to match an expected density
    if n is None:
        n = int(np.ceil(density*np.pi*r**2))
    if rng is None:
        rng = np.random.default_rng()

    c = np.empty((0,2), dtype=float)
    overdraw = 1.5
    while c.shape[0] != n:
        # Draw 50% more points than needed within the unit square
        c = rng.uniform(low=-1, high=1, size=(int(n*overdraw),2))
        # Find those within the r = 1
        indx = c[:,0]**2 + c[:,1]**2 < 1
        c = c[indx][:n]
        # Increase overdraw for next iteration
        overdraw *= 1.1
    
    return r*c[:,0], r*c[:,1]


def parse_fobos_mode(mode):
    """
    Parse and validate the provided FOBOS spectrograph mode.

    Args:
        mode (:obj:`int`, array-like):
            The mode assignment for each spectrograph.  Can be 0 to turn
            off all apertures for a given spectrograph, but otherwise
            must be 1 for the single-fiber apertures, 2 for multi-IFU
            mode, or 3 for monolithic IFU mode.  If a single integer,
            all spectrographs are put in the same mode; otherwise, must
            provide the mode for each of the 3 spectrographs separately.
    
    Returns:
        :obj:`list`: A list with the integer mode identifier for each of
        the three spectrographs.
    """
    # Check the mode input
    _mode = np.atleast_1d(mode).astype(int)
    if _mode.size == 1:
        _mode = np.repeat(_mode, 3)
    if _mode.size != 3:
        raise ValueError(f'Mode not understood: {mode}.  Must enter single mode for all '
                         f'spectrographs, or the mode for each of the three spectrographs.')
    if any([m not in [0,1,2,3] for m in _mode]):
        raise ValueError('Mode must be either 1 for single-fiber apertures or 2 for IFU.')
    return _mode


# TODO:
#    - Convert this into a class to avoid needing to read the file
#      multiple times.
#    - In IFU mode, identify single fiber apertures as sky apertures.
#      In mixed mode observations, avoid using single fiber apertures to
#      observe single-fiber targets.
def fobos_aperture_grid(mode=1, payload='science', baseline=True, config=1):
    """
    Read the home positions for FOBOS apertures.

    .. warning::
        
        Filename with the data is hard-coded.

    Args:
        mode (:obj:`int`, array-like, optional):
            The mode assignment for each spectrograph.  Can be 0 to turn
            off all apertures for a given spectrograph, but otherwise
            must be 1 for the single-fiber apertures, 2 for multi-IFU
            mode, or 3 for monolithic IFU mode.  If a single integer,
            all spectrographs are put in the same mode; otherwise, must
            provide the mode for each of the 3 spectrographs separately.
        payload (:obj:`str`, optional):
            The payload coordinates to return.  Must be ``'all'`` for
            any bundle type, ``'science'`` for science apertures,
            ``'calib'`` for flux-calibration bundles, or ``'guide'`` for
            imaging guide bundles.
        baseline (:obj:`bool`, optional):
            Flag to only return the apertures in the CoDR baseline
            configuration.
        config (:obj:`int`, optional):
            Configuration selection for mapping aperture modules to each
            spectrograph.  Must be either: (1) modules mapped to each
            spectrograph fill coherent focal-plane zones; or (2) modules
            mapped to each spectrograph are spread over the entire focal
            plane.
    
    Returns:
        :obj:`tuple`: Four 1D `numpy.ndarray`_ objects with the
        spectrograph ID for each aperture, the aperture type (0=single
        fiber; 2=IFU; 3=imaging bunde; 4=flux-calibration bundle), and
        the relative x and y positions of the apertures in arcmin
        relative to the focal-plane center.
    """
    # Check payload
    if payload not in ['all', 'science', 'calib', 'guide']:
        raise ValueError(f'Payload {payload} unknown.  Must be science, calib, or guide.')

    # Parse the spectrograph mode
    _mode = parse_fobos_mode(mode)
    
    # Select the desired spectrograph configuration column
    if config not in [1, 2]:
        raise ValueError(f'Spectrograph configuration {config} unknown.  Must be 1 or 2.')
    sc = 2 if config == 1 else 4
    mc = sc + 1

    # Load data
    fp = np.loadtxt('fp_layout_codr_starbugs.txt', dtype=float)

    # Loop through the modes of each spectrograph
    indx = np.zeros(fp.shape[0], dtype=bool)
    spec_id = fp[:,sc].astype(int)
    spec_mode = fp[:,mc].astype(int)
    payload_id = fp[:,8].astype(int)
    for i in range(3):
        if _mode[i] == 0:
            continue
        selected_payloads = [0, 1, 2] if payload in ['all', 'science'] else []
        if payload in ['all', 'guide']:
            selected_payloads += [3]
        if payload in ['all', 'calib']:
            selected_payloads += [4]
        indx |= (spec_id == i+1) & np.isin(payload_id, selected_payloads) \
                    & (spec_mode & (1 << _mode[i]-1) != 0)

    # Restrict based on baseline
    if baseline:
        indx &= fp[:,9].astype(int) == 1

    return spec_id[indx], payload_id[indx], fp[indx,6], fp[indx,7]


#def get_cos(target_file, center_coo):
#    center = SkyCoord(center_coo[0]*units.deg, center_coo[1]*units.deg, frame='icrs')
#    
#
#    db = np.genfromtxt(target_file)
#    targets = SkyCoord(db[:,1]*units.deg, db[:,2]*units.deg, frame='icrs')
#    targ_offsets = targets.transform_to(SkyOffsetFrame(origin=center))
#
#    r = np.sqrt(targ_offsets.lon.to('arcmin').value**2 + targ_offsets.lat.to('arcmin').value**2)
#    indx = r < 10
#    #print(np.sum(r < 10))
#        
#    galx = targ_offsets.lon.to('arcmin').value[indx]
#    galy = targ_offsets.lat.to('arcmin').value[indx]
#        
#    return galx, galy

#def get_allcos(files, center_coo):
#    center = SkyCoord(center_coo[0]*units.deg, center_coo[1]*units.deg, frame='icrs')
#    all_db = None #insert code here
#    
#    for i,f in enumerate(files):
#        db = np.genfromtxt(f)
#        galid = np.full(db.shape[0],i,dtype=int)
#        
#        print(galid.shape)
#        print(db.shape)
#        db = np.hstack((galid.reshape(-1,1),db))
#       
#        if all_db is None:
#            all_db = db
#        else:
#            all_db = np.append(all_db,db, axis=0)
#    
#    print(all_db.shape)
#    targets = SkyCoord(all_db[:,2]*units.deg, all_db[:,3]*units.deg, frame='icrs')
#    targ_offsets = targets.transform_to(SkyOffsetFrame(origin=center))
#
#    r = np.sqrt(targ_offsets.lon.to('arcmin').value**2 + targ_offsets.lat.to('arcmin').value**2)
#    indx = r < 10
#
#    galx = targ_offsets.lon.to('arcmin').value[indx]
#    galy = targ_offsets.lat.to('arcmin').value[indx]
#        
#    return galx, galy, all_db[indx,0]

# TODO:
#   - Set fraction to assign to sky
#   - Write method that assigns sky fibers, guide bundles, calibration
#     bundles
def configure_observations(objx, objy, aptype, mode=None):
    """
    Construct a series of observations to observe a set of targets with
    a fixed field center.

    Args:
        objx (`numpy.ndarray`_):
            Cartesian x coordinate in tangent plane projection of object
            coordinates relative to the pointing center.  Shape must be
            1D and match ``objy``.
        objy (`numpy.ndarray`_):
            Cartesian y coordinate in tangent plane projection of object
            coordinates relative to the pointing center. Shape must be
            1D and match ``objx``.
        aptype (:obj:`int`, `numpy.ndarray`_):
            The aperture type for each object.  Can provide a single
            integer for all targets or an aperture type for each object.
            The aperture type must be 0 for a single-fiber aperture or 2
            for a 37-fiber IFU.
        mode (:obj:`int`, array-like, optional):
            The mode assignment for each spectrograph.    Can be 0 to turn
            off all apertures for a given spectrograph, but otherwise
            must be 1 for the single-fiber apertures, 2 for multi-IFU
            mode, or 3 for monolithic IFU mode.  If a single integer,
            all spectrographs are put in the same mode; otherwise must
            provide the mode for each of the 3 spectrographs separately.

    Returns:
        stuff
    """
    # Check input
    if objx.ndim > 1:
        raise ValueError('Input coordinate arrays must be 1D.')
    if objy.shape != objx.shape:
        raise ValueError('Input coordinate arrays must have the same shape.')
    nobj = objx.size
    _aptype = np.atleast_1d(aptype)
    if _aptype.size == 1:
        _aptype = np.full(nobj, aptype, dtype=int)
    if _aptype.shape != objx.shape:
        raise ValueError('For per-object aperture types, must have same shape as coordinates.')
    if not np.all(np.isin(_aptype, [0,2])):
        raise ValueError('Aperture types must be 0 (single-fiber) or 2 (IFU).')

    # If mode is None, first try to set it based on the aperture types
    if mode is None:
        if np.all(_aptype == 0):
            # All aperture are in single-fiber mode, so set the
            # instrument configuration accordingly
            mode = 1
        elif np.all(_aptype == 2):
            # All aperture are in IFU mode, so set the instrument
            # configuration accordingly
            mode = 2

    # If mode is still None, that means there is a mix of IFU and
    # single-fiber targets, and we need to find the best of 6 possible
    # mixed-mode configurations
    if mode is None:
        raise NotImplementedError('Cannot yet determine best mixed mode configuration.')

    # Parse the mode
    _mode = np.array(parse_fobos_mode(mode))

    # Objects to track assigned object and aperture IDs for each
    # observation
    obs_obj = []
    obs_ap = []
    obs_nap = []
    obs_mode = []
    _objx = objx.copy()
    _objy = objy.copy()
    _aptype = aptype.copy()
    obj_id = np.arange(nobj, dtype=int)
  
    # TODO: This can be moved inside the loop
    # Read the starbug coordinates
    spec_id, payload_id, apx, apy \
            = fobos_aperture_grid(mode=_mode, payload='science', config=1)
    ap_id = np.arange(apx.size, dtype=int)

    # Itertatively construct observations that each, as much as
    # possible, assign apertures to targets until no more targets can be
    # observed.
    while _objx.size > 0:

        # TODO: This ordering means single fibers are given priority
        # over IFUs

        _a_obj = None
        _a_ap = None

        obs_nap += [0]

        # Assign the single-fiber apertures to appropriate targets
        use_obj = _aptype == 0
        use_ap = payload_id == 0
        obs_nap[-1] = np.sum(use_ap)
        if np.any(use_obj) and np.any(use_ap):
            _a_obj, _a_ap = assign_apertures(_objx, _objy, apx, apy,
                                             ignore_obj=np.where(np.logical_not(use_obj))[0],
                                             ignore_ap=np.where(np.logical_not(use_ap))[0])
            if _a_obj.size == 0:
                _a_obj = None
                _a_ap = None

        # Assign the IFU apertures to appropriate targets
        use_obj = _aptype == 2
        use_ap = payload_id == 2
        obs_nap[-1] += np.sum(use_ap)
        if np.any(use_obj) and np.any(use_ap):
            _ifu_obj, _ifu_ap = assign_apertures(_objx, _objy, apx, apy,
                                                 ignore_obj=np.where(np.logical_not(use_obj))[0],
                                                 ignore_ap=np.where(np.logical_not(use_ap))[0],
                                                 allocated_obj=_a_obj, allocated_ap=_a_ap)
            if _ifu_obj.size > 0:
                _a_obj = np.append(_a_obj, _ifu_obj)
                _a_ap = np.append(_a_ap, _ifu_ap)

        if _a_obj is None or _a_obj.size == 0:
            # No apertures could be assigned so we're done!
            break

        # Add the object IDs, aperture IDs, and spectrograph mode for
        # this observation
        obs_obj += [obj_id[_a_obj]]
        obs_ap += [ap_id[_a_ap]]
        obs_mode += [_mode]

        # Remove the assigned objects from those available in the next
        # iteration
        _objx = np.delete(_objx, _a_obj)
        _objy = np.delete(_objy, _a_obj)
        _aptype = np.delete(_aptype, _a_obj)
        obj_id = np.delete(obj_id, _a_obj)

    return obs_obj, obs_nap, obs_ap, obs_mode


def main_old():
    
    #set files
    files = ['quiescent_sample.txt', 'starbursting_sample.txt', 'highz_sample.txt']
    
    
    center = [150.1,2.2]
    #n = [[100, 200, 400, 800], [2000, 4000, 8000, 16000]]
    
    #empty arrays for peak values
    
    #plots outside of loop
    #fig1,((ax4, ax5), (ax6, ax7)) = plt.subplots(2,2)
    fig1,(ax1, ax2, ax3) = pyplot.subplots(3,1, sharex=True)
    
    galx, galy, galid = get_allcos(files,center)
    indx = (galid==0) | (galid == 2)
    print(galx.shape, indx.shape)
    
    a_single, galx_s, galy_s, fibx, fiby, niter, SID1 = get_obs([1,0,1],galx[indx],galy[indx], J=galid[indx]) #mode 1 (1 & 3 files)
    
    indx = (galid ==1)
    a_IFU, galx_m, galy_m, fibx_m, fiby_m, niter_m, SID1_m = get_obs([0,2,0], galx[indx], galy[indx], J=galid[indx]) #mode 2  (2 file)
    
    print([len(a) for a in a_single])
    print([len(a) for a in a_IFU])
    print(galx.shape, galy.shape, galid.shape)
    
    total_obs = max(len(a_single),len(a_IFU)) #list of ids,multiply max used to use to create list thats as long as longest a_* list
    a = [None]*total_obs
    
    for i in range (total_obs):
        if i < len(a_single):
            if a[i] is None:
                a[i]=a_single[i]
            else:
                a[i]=np.append(a[i],a_single[i])
            print(len(a[i]))
 
        if i < len(a_IFU):
            if a[i] is None:
                a[i]=a_IFU[i]
            else:
                a[i]=np.append(a[i],a_IFU[i])
            print(len(a[i]))
    
    print([len(_a) for _a in a])
    #no of observations
    niter = total_obs
    arr = np.arange(niter)+1

    #cumulative sum of all observed stars in each observation
    c = [_a.size for _a in a]
    c = np.asarray(c)
    n_assign_stars = galx.size
    c_ = np.cumsum(c)
    
    #FOM variables
    n_avail_fib = fibx.size + fibx_m.size
    comp = c_/n_assign_stars
    eff = c/n_avail_fib
    mean_eff = np.cumsum(eff)/(arr)
    FOM = np.sqrt(comp)*mean_eff
    
    #completeness plot
    ax1.plot(arr, c_/n_assign_stars, marker = '.', label = f'{n_assign_stars}')
    ax1.set_title('Completeness of Targets')
    ax1.set_ylabel('Completeness (%)')
    ax1.legend()
    #ax1.text(0.05,0.01,f'n_assign = {n_assign_stars}', transform = ax1.transAxes)
        
    #efficency plot
    ax2.plot(arr, c/n_avail_fib, marker = '.', label = f'{n_avail_fib}') #need to change if modes change because fibx changes
    ax2.set_title('Efficency of Assigned Targets')
    ax2.set_ylabel('Efficency (%)')
    ax2.legend()
    #ax2.text(0.05,0.01,f'n_obs = {c.size}', transform = ax2.transAxes)
        
    #FOM plot
    ax3.plot(arr, FOM, marker = '.')
    ax3.set_title('Figure of Merit')
    ax3.set_xlabel('No. of Observations')
    ax3.set_ylabel('FOM')
         
    fig1.savefig(f'test_1')
    

    #FOM peak values
    FOM_ = np.max(FOM)
    arr_ = arr[np.argmax(FOM)]
    n_stars = galx.size
    n_fib = n_avail_fib
    
    #comp peak values
    comp_peak = np.max(comp)
    arr_c = comp[np.argmax(FOM)]
    
    #eff peak values
    eff_peak = np.max(eff)
    arr_e = eff[np.argmax(FOM)]
    arr_me = mean_eff[np.argmax(FOM)]
       
 
    #mode_name += [int(''.join([str(m) for m in mode))]

    np.savetxt(f'FOM_values1.txt',
              np.column_stack((121,
                               n_stars,
                               n_fib,
                               c_[-1],
                               FOM_,
                               arr_c,
                               arr_e, 
                               arr_me)), 
                               fmt = '%10d %7d %13d %7d %15.4f %10.4f %11.4f %11.4f',
                               header = 'Mode''\nNo. of Targets Observed\nNo. of Fibers Available\n' 
                               + f"{'Mode name':>8} {'# of Targets':>0} {'  # of Fibers':>10} {'Observed Targets':>0} {'  FOM Peak':>9} {'  Comp. Peak':>10}{'  Eff. Peak':>10} {'  Mean Eff. Peak':>10}")

def show_focal_plane(mode, config, baseline=True):

        
    spec_id, payload_id, apx, apy \
            = fobos_aperture_grid(mode=mode, payload='all', baseline=baseline, config=config)
    
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim([-14, 14])
    ax.set_ylim([-14, 14])
    # Single fibers
    indx = payload_id == 0
    ax.scatter(apx[indx], apy[indx], marker='.', s=30, lw=0, color='C1', zorder=3)
    # IFUs
    indx = payload_id == 2
    ax.scatter(apx[indx], apy[indx], marker='s', s=20, lw=0.5, color='C2', zorder=3)
    # Guide bundles
    indx = payload_id == 3
    ax.scatter(apx[indx], apy[indx], marker='o', s=20, lw=0.5, color='C0', zorder=3)
    # Flux-calibration bundles
    indx = payload_id == 4
    ax.scatter(apx[indx], apy[indx], marker='x', s=20, lw=0.5, color='C4', zorder=3)
    for x, y in zip(apx, apy):
        ax.add_patch(patches.Circle((x,y), radius=138/60., facecolor='k', edgecolor='none',
                                    zorder=0, alpha=0.05))
    ax.add_patch(patches.Circle((0.,0.), radius=10., facecolor='none', edgecolor='C3', zorder=1))
    pyplot.show()
    

def assign_apertures_plot(objx, objy, apx, apy, assigned_obj, assigned_ap, collision=1/6.):

    xlim = [min(np.amin(objx[assigned_obj])-3*collision,
                np.amin(objy[assigned_obj])-3*collision),
            max(np.amax(objy[assigned_obj])+3*collision,
                np.amax(objx[assigned_obj])+3*collision)]

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)

    for _x, _y in zip(objx[assigned_obj], objy[assigned_obj]):
        ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='C0',
                                    edgecolor='none', zorder=0, alpha=0.2))
    if len(assigned_ap) != apx.size:
        unassigned = np.setdiff1d(np.arange(apx.size), assigned_ap)
        print(apx[unassigned])
        print(apy[unassigned])
        for _x, _y in zip(apx[unassigned], apy[unassigned]):
            ax.add_patch(patches.Circle((_x,_y), radius=collision/2., facecolor='C3',
                                        edgecolor='none', zorder=0, alpha=0.2))

    ax.scatter(objx[assigned_obj], objy[assigned_obj], marker='.', color='C0', s=30, lw=0)
    if len(assigned_obj) != objx.size:
        unassigned = np.setdiff1d(np.arange(objx.size), assigned_obj)
        ax.scatter(objx[unassigned], objy[unassigned], marker='.', color='C3', s=30, lw=0)
    
    pyplot.show()


def main():

    rng = np.random.default_rng(99)

    ndens = 5
    nsim = 10
    density = np.geomspace(2.5, 40, ndens)
    nobj = np.zeros((nsim, ndens), dtype=int)
    nobs = np.zeros((nsim, ndens), dtype=int)
    nalloc = np.empty((nsim, ndens), dtype=object)
    nap = np.empty((nsim, ndens), dtype=object)
    completeness = np.empty((nsim, ndens), dtype=object)
    efficiency = np.empty((nsim, ndens), dtype=object)
    mean_efficiency = np.empty((nsim, ndens), dtype=object)

    for j in range(nsim):
        for i in range(ndens):
            print(f'Density {i+1}/{ndens}', end='\r')

            objx, objy = random_targets(10., density=density[i], rng=rng)
            nobj[j,i] = objx.size
            aptype = np.zeros(nobj[j,i], dtype=int)

            obs_obj, obs_nap, obs_ap, obs_mode = configure_observations(objx, objy, aptype)

            nobs[j,i] = len(obs_obj)
            nalloc[j,i] = np.array([len(o) for o in obs_obj])
            nap[j,i] = np.array(obs_nap)
            completeness[j,i] = np.cumsum(nalloc[j,i])/nobj[j,i]
            efficiency[j,i] = nalloc[j,i]/nap[j,i]
            mean_efficiency[j,i] = np.cumsum(efficiency[j,i])/(np.arange(nobs[j,i])+1)

        print(f'Density {ndens}/{ndens}')

    obs_mask = np.empty((nsim, ndens), dtype=object)
    max_nobs = np.amax(nobs, axis=0)
    for j in range(nsim):
        for i in range(ndens):
            obs_mask[j,i] = np.zeros(max_nobs[i], dtype=bool)
            if nobs[j,i] < max_nobs[i]:
                obs_mask[j,i][nobs[j,i]:] = True
                nalloc[j,i] = np.append(nalloc[j,i], np.zeros(max_nobs[i] - nobs[j,i]))
                nap[j,i] = np.append(nap[j,i], np.zeros(max_nobs[i] - nobs[j,i]))
                completeness[j,i] = np.append(completeness[j,i],
                                              np.zeros(max_nobs[i] - nobs[j,i]))
                efficiency[j,i] = np.append(efficiency[j,i], np.zeros(max_nobs[i] - nobs[j,i]))
                mean_efficiency[j,i] = np.append(mean_efficiency[j,i],
                                                 np.zeros(max_nobs[i] - nobs[j,i]))

    m_completeness = np.empty(ndens, dtype=object)
    m_efficiency = np.empty(ndens, dtype=object)
    m_fom = np.empty(ndens, dtype=object)
    
    s_completeness = np.empty(ndens, dtype=object)
    s_efficiency = np.empty(ndens, dtype=object)
    s_fom = np.empty(ndens, dtype=object)
    
    for i in range(ndens):
        comp = np.ma.MaskedArray(np.stack(completeness[:,i]), mask=np.stack(obs_mask[:,i]))
        m_completeness[i] = np.ma.mean(comp, axis=0).filled(0.0)
        s_completeness[i] = np.ma.std(comp, axis=0).filled(0.0)

        arr = np.ma.MaskedArray(np.stack(efficiency[:,i]), mask=np.stack(obs_mask[:,i]))
        m_efficiency[i] = np.ma.mean(arr, axis=0).filled(0.0)
        s_efficiency[i] = np.ma.std(arr, axis=0).filled(0.0)

        arr = np.ma.sqrt(comp) \
                * np.ma.MaskedArray(np.stack(mean_efficiency[:,i]), mask=np.stack(obs_mask[:,i]))
        m_fom[i] = np.ma.mean(arr, axis=0).filled(0.0)
        s_fom[i] = np.ma.std(arr, axis=0).filled(0.0)

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.2, 0.69, 0.6, 0.3])
    ax.minorticks_on()
    ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
    ax.set_xlim(0.8, np.amax(nobs) + 0.2)
    ax.set_ylim(0., 1.05)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    for i in range(ndens):
        ax.fill_between(np.arange(max_nobs[i])+1, m_completeness[i] + s_completeness[i],
                        y2=m_completeness[i] - s_completeness[i],
                        color=f'C{i}', lw=0, alpha=0.2, zorder=1)
        ax.scatter(np.arange(max_nobs[i])+1, m_completeness[i], color=f'C{i}',
                   marker='.',  s=30, lw=0, zorder=3)
        ax.plot(np.arange(max_nobs[i])+1, m_completeness[i], color=f'C{i}', zorder=3)

    ax = fig.add_axes([0.2, 0.38, 0.6, 0.3])
    ax.minorticks_on()
    ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
    ax.set_xlim(0.8, np.amax(nobs) + 0.2)
    ax.set_ylim(0., 1.05)
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    for i in range(ndens):
        ax.fill_between(np.arange(max_nobs[i])+1, m_efficiency[i] + s_efficiency[i],
                        y2=m_efficiency[i] - s_efficiency[i],
                        color=f'C{i}', lw=0, alpha=0.2, zorder=1)
        ax.scatter(np.arange(max_nobs[i])+1, m_efficiency[i], color=f'C{i}',
                   marker='.', s=30, lw=0, zorder=3)
        ax.plot(np.arange(max_nobs[i])+1, m_efficiency[i], color=f'C{i}', zorder=3)

    ax = fig.add_axes([0.2, 0.07, 0.6, 0.3])
    ax.minorticks_on()
    ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
    ax.set_xlim(0.8, np.amax(nobs) + 0.2)
    ax.set_ylim(0., 1.05)
    for i in range(ndens):
        ax.fill_between(np.arange(max_nobs[i])+1, m_fom[i] + s_fom[i],
                        y2=m_fom[i] - s_fom[i],
                        color=f'C{i}', lw=0, alpha=0.2, zorder=1)
        ax.scatter(np.arange(max_nobs[i])+1, m_fom[i],
                   color=f'C{i}', marker='.', s=30, lw=0, zorder=3)
        ax.plot(np.arange(max_nobs[i])+1, m_fom[i],
                color=f'C{i}', zorder=3)

    pyplot.show()




    embed()
    exit()

#    spec_id, payload_id, apx, apy \
#            = fobos_aperture_grid(mode=1, payload='science', config=1)
#
#    assigned_obj, assigned_ap = assign_apertures(objx, objy, apx, apy)
#
#    assign_apertures_plot(objx, objy, apx, apy, assigned_obj, assigned_ap)

 

if __name__ == '__main__':
    main()












