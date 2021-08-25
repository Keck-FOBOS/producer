"""
Construct a set of observations by assigning fibers to targetsAllocate apertures to targets.

Contains code originally written by 2021 Akamai intern, Brittany Ann
Ramos.

.. include:: ../include/links.rst
"""
import io
import sys
import time
import warnings
from pathlib import Path
from configparser import ConfigParser

from IPython import embed 

import numpy

from astropy import table

from .deploy import FOBOSApertures
from .allocate import assign_apertures
from .targets import parse_targets

# TODO:
#   - Set fraction to assign to sky
#   - Write method that assigns sky fibers, guide bundles, calibration
#     bundles
def configure_observations(objx, objy, aptype, mode=None, max_nobs=None):
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
            The aperture type must be 0 for a single-fiber aperture or 1
            for a 37-fiber IFU.
        mode (:obj:`int`, array-like, optional):
            The mode assignment for each spectrograph.    Can be 0 to turn
            off all apertures for a given spectrograph, but otherwise
            must be 1 for the single-fiber apertures, 2 for multi-IFU
            mode, or 3 for monolithic IFU mode.  If a single integer,
            all spectrographs are put in the same mode; otherwise must
            provide the mode for each of the 3 spectrographs separately.
        max_nobs (:obj:`int`, optional):
            Impose a maximum number of observations to configure.  If None,
            observations will be configured until no more objects are within the
            field of view.

    Returns:
        stuff
    """
    # Check input
    if objx.ndim > 1:
        raise ValueError('Input coordinate arrays must be 1D.')
    if objy.shape != objx.shape:
        raise ValueError('Input coordinate arrays must have the same shape.')
    nobj = objx.size
    _aptype = numpy.atleast_1d(aptype)
    if _aptype.size == 1:
        _aptype = numpy.full(nobj, aptype, dtype=int)
    if _aptype.shape != objx.shape:
        raise ValueError('For per-object aperture types, must have same shape as coordinates.')
    if not numpy.all(numpy.isin(_aptype, [0,1])):
        raise ValueError('Aperture types must be 0 (single-fiber) or 2 (IFU).')

    # If mode is None, first try to set it based on the aperture types
    if mode is None:
        if numpy.all(_aptype == 0):
            # All aperture are in single-fiber mode, so set the
            # instrument configuration accordingly
            mode = 1
        elif numpy.all(_aptype == 1):
            # All aperture are in IFU mode, so set the instrument
            # configuration accordingly
            mode = 2

    # If mode is still None, that means there is a mix of IFU and
    # single-fiber targets, and we need to find the best of 6 possible
    # mixed-mode configurations
    if mode is None:
        raise NotImplementedError('Cannot yet determine best mixed mode configuration.')

    # Setup the apertures
    ap = FOBOSApertures(mode=mode)

    # Objects to track assigned object and aperture IDs for each
    # observation
    obs_obj = []
    obs_ap = []
    obs_nap = []
    obs_mode = []
    _objx = objx.copy()
    _objy = objy.copy()
    _aptype = aptype.copy()
    obj_id = numpy.arange(nobj, dtype=int)

    # Get the aperture ID, coordinates, and payload type
    ap_id = ap.id.copy()
    apx, apy = ap.coo.T.copy()
    payload = ap.payload.copy()
    # TODO: This can be moved inside the loop
    # Read the starbug coordinates
    ap_indx = ap.select('science')

    # Total number of objects within the field-of-view
    n_in_fov = numpy.sum(objx**2 + objy**2 < (ap.fov*60/2)**2)

    # Itertatively construct observations that each, as much as
    # possible, assign apertures to targets until no more targets can be
    # observed.
    niter = 0
    while _objx.size > 0 and (max_nobs is None or niter < max_nobs):

        # TODO: This ordering means single fibers are given priority
        # over IFUs

        _a_obj = None
        _a_ap = None

        # Only consider those objects in the field of view
        in_fov = _objx**2 + _objy**2 < (ap.fov*60/2)**2

        # Assign the single-fiber apertures to appropriate targets
        use_obj = (_aptype == 0) & in_fov
        use_ap = (payload == 0) & ap_indx
        _nap = numpy.sum(use_ap)
        if numpy.any(use_obj) and numpy.any(use_ap):
            _a_obj, _a_ap = assign_apertures(_objx, _objy, apx, apy,
                                             ignore_obj=numpy.where(numpy.logical_not(use_obj))[0],
                                             ignore_ap=numpy.where(numpy.logical_not(use_ap))[0])
            if _a_obj.size == 0:
                _a_obj = None
                _a_ap = None

        # Assign the IFU apertures to appropriate targets
        use_obj = (_aptype == 1) & in_fov
        use_ap = (payload == 1) & ap_indx
        _nap += numpy.sum(use_ap)
        if numpy.any(use_obj) and numpy.any(use_ap):
            _ifu_obj, _ifu_ap \
                    = assign_apertures(_objx, _objy, apx, apy,
                                       ignore_obj=numpy.where(numpy.logical_not(use_obj))[0],
                                       ignore_ap=numpy.where(numpy.logical_not(use_ap))[0],
                                       allocated_obj=_a_obj, allocated_ap=_a_ap)
            if _ifu_obj.size > 0:
                _a_obj = numpy.append(_a_obj, _ifu_obj)
                _a_ap = numpy.append(_a_ap, _ifu_ap)

        if _a_obj is None or _a_obj.size == 0:
            # No apertures could be assigned so we're done!
            break

        # Add the object IDs, aperture IDs, and spectrograph mode for
        # this observation
        obs_obj += [obj_id[_a_obj]]
        obs_nap += [_nap]
        obs_ap += [ap_id[_a_ap]]
        obs_mode += [ap.mode]

        # Remove the assigned objects from those available in the next
        # iteration
        _objx = numpy.delete(_objx, _a_obj)
        _objy = numpy.delete(_objy, _a_obj)
        _aptype = numpy.delete(_aptype, _a_obj)
        obj_id = numpy.delete(obj_id, _a_obj)

        niter += 1

    return n_in_fov, obs_obj, obs_nap, obs_ap, obs_mode


def report_configurations(n_in_fov, obs_obj, obs_nap, obs_ap, obs_mode):
    """
    Construct a report of the observation configurations.
    """
    tab = table.Table()
    # Number of observations
    nobs = len(obs_obj)
    # Pointing number
    tab['Pointing'] = numpy.arange(nobs)+1
    # Number of available apertures in each observation
    tab['Avail'] = numpy.array(obs_nap)
    # Number of allocations in each observation
    tab['Alloc'] = numpy.array([len(o) for o in obs_obj])
    # Fractional growth of observed targets
    tab['Compl'] = numpy.cumsum(tab['Alloc'])/n_in_fov
    tab['Compl'].format = '.4f'
    # Fraction of all apertures assigned
    tab['Eff'] = tab['Alloc']/tab['Avail']
    tab['Eff'].format = '.4f'
    # Mean fraction of assigned apertures for this and all previous observations
    tab['MeanEff'] = numpy.cumsum(tab['Eff'])/(numpy.arange(nobs)+1)
    tab['MeanEff'].format = '.4f'

    print(f'Total number of targets available: {n_in_fov}')
    print(f'Total number of pointings: {nobs}')
    print('')
    tab.write(sys.stdout, format='ascii.fixed_width_two_line', delimiter=' ')


def write_configurations(root, ra, dec, center, obs_obj, obs_ap, obs_mode, objid=None, path=None,
                         ndig=None, tight=False, target_file=None, ra_c=None, dec_c=None):
    """
    Write a set of configuration files for each FOBOS observation.

    Args:
        root (:obj:`str`):
            The root name for all output files.
        ra (`numpy.ndarray`_):
            Right ascension coordinates for all considered objects.
        dec (`numpy.ndarray`_):
            Declination coordinates for all considered objects.  Shape must
            match ``ra``.
        center (:obj:`tuple`):
            RA and DEC coordinates for the FOBOS pointing center.
        obs_obj (:obj:`list`):
            List of `numpy.ndarray`_ objects identifying the indices of the
            objects observed from the provided list of coordinates.  The number
            of items in the list sets the number of revisits to the same
            pointing.  This is the same as the second object returned by
            :func:`~producer.plan.configure_observations`.
        obs_ap (:obj:`list`):
            List of `numpy.ndarray`_ objects identifying the indices of the
            FOBOS apertures used for each object observed.  List length must
            match ``obs_obj``.  This is the same as the fourth object returned
            by :func:`~producer.plan.configure_observations`.  The aperture
            indices must match indices when instantiating a
            :class:`~producer.deploy.FOBOSApertures` object in the specified
            mode (``obs_mode``).
        obs_mode (:obj:`list`):
            List of `numpy.ndarray`_ objects identifying the FOBOS mode; see 
            :class:`~producer.deploy.FOBOSApertures`.  List length must match
            ``obs_obj``.  This is the same as the last object returned by
            :func:`~producer.plan.configure_observations`.
        objid (`numpy.ndarray`_, optional):
            An array with identifiers for each object.  Each array element must
            convert directly to a string.  Uniqueness is not checked.  Shape
            must match ``ra``.  If None, just set to the 0-indexed array index.
        path (:obj:`str`, optional):
            Root path for all output files.  If None, either set to the parent
            path provided by ``root`` or set to the current directory.
        ndig (:obj:`int`, optional):
            Number of digits to use for the observation number in the output
            file names.  If None, this is set by the number of configurations to
            write.  E.g., 9 observations or less yield ``ndig=1``, 10-99
            observations yield ``ndig=2``, etc.
        tight (:obj:`bool`, optional):
            Output the configuration in "tight" format, where unallocated
            apertures are not included.
        target_file (:obj:`str`, optional):
            Name of the file with the original targets.  If provided, will be
            included in header of output configuration files.
        ra_c (:obj:`int`, optional):
            1-indexed column number with the RA coordinates in ``target_file``.
            Ignored if ``target_file`` is None.
        dec_c (:obj:`int`, optional):
            1-indexed column number with the DEC coordinates in ``target_file``.
            Ignored if ``target_file`` is None.
    """
    # Check the coordinate arrays
    if ra.ndim > 1:
        raise ValueError('Input coordinate arrays must be 1D.')
    if dec.shape != ra.shape:
        raise ValueError('Shape of coordinate arrays differ.')
    nobj = ra.size
    objid_type = 'name'
    if objid is None:
        objid = numpy.arange(nobj).astype(str)
        objid_type = 'index'
    if objid.shape != ra.shape:
        raise ValueError('Object ID array does not match shape of coordinate arrays.')
    objidlen = numpy.amax([len(o) for o in objid])

    _center = numpy.atleast_1d(center).ravel()
    if _center.size != 2:
        raise ValueError('Center coordinates must provide RA and DEC only.')

    # Check the observation lists
    nobs = len(obs_obj)
    if len(obs_ap) != nobs:
        raise ValueError(f'Incorrect number of aperture arrays; expected {nobs}, got '
                         f'{len(obs_ap)}.')
    if len(obs_mode) != nobs:
        raise ValueError(f'Incorrect number of instrument modes; expected {nobs}, got '
                         f'{len(obs_mode)}.')
    for indx in obs_obj:
        if numpy.any(indx >= nobj):
            raise ValueError('Object selection indices out of bounds of coordinate arrays.')

    # Set the file name number of digits
    if ndig is None:
        ndig = int(numpy.ceil(numpy.log10(nobs+1)))

    # Construct and check the output root
    _root = Path(root).resolve()
    if _root.is_dir():
        raise ValueError('Provided root is a directory.  Must provide leading name of files.')
    if path is None:
        _path = _root.parent
        _root = _root.name
    else:
        _path = Path(path).resolve()
        if not _path.is_dir():
            warnings.warn(f'Request path, {path}, currently does not exist and will be created.')
        _root = _root.name

    for i in range(nobs):
        ofile = str(_path / f'{_root}_{i+1:0{ndig}}.db')
        ap = FOBOSApertures(mode=obs_mode[i])
        if numpy.any(obs_ap[i] >= ap.nap):
            raise ValueError(f'Aperture selection indices for observation {i+1} are invalid.')

        cfg = ConfigParser()
        cfg['FOBOS_MODE'] = {'version': ap.version,
                             'mode': ','.join(ap.mode.astype(str)),
                             'baseline': ap.baseline,
                             'design': ap.config
                             }

        cfg['TARGETS'] = {'source': 'Unknown' if target_file is None else target_file,
                          'ra_col': 'None' if ra_c is None else ra_c,
                          'dec_col': 'None' if dec_c is None else dec_c,
                          'id_type': objid_type,
                          'allocated': len(obs_obj[i])
                          }

        cfg['POINTING'] = {'ra': f'{_center[0]:12.8f}',
                           'dec': f'{_center[1]:12.8f}'
                          }
        with io.StringIO() as f:
            cfg.write(f)
            cfg_str = f.getvalue()

        header_text = 'FOBOS configuration file\n\n' \
                      f'Generated: {time.strftime("%a %d %b %Y %H:%M:%S", time.localtime())}\n\n'
        header_text += '-'*70 + '\n\n'
        header_text += cfg_str
        header_text += '-'*70 + '\n\n'
        header_text += 'Columns are:\n' \
                       '    MID - Module ID number\n' \
                       '    SID - Spectrograph ID number\n' \
                       '    BID - Starbug ID number\n' \
                       '    Type - Starbug payload type: (0) single-fiber, (1) 37-fiber \n' \
                       '           IFU, (3) imaging bundle, (4) flux-calibration bundle.\n' \
                       '    OBJID - Object ID \n' \
                       '    RA - Target right ascension (decimal degrees)\n' \
                       '    Dec - Target declination (decimal degrees)\n\n' \
                       f'{"MID":>3}   {"SID":>3}   {"BID":>3}   {"Type":>4}   ' \
                       f'{"OBJID":>{objidlen}}   {"RA":>12}   {"DEC":>11}'

        tab = empty_configuration_table(nrows=ap.nap, objidlen=objidlen)
        # NOTE: The syntax "tab['MID'][:]" is needed here so that the table
        # column doesn't loose it's other attributes (format, description, etc.)
        tab['MID'][:] = ap.mid
        tab['SID'][:] = ap.spc
        tab['BID'][:] = ap.bid
        tab['Type'][:] = ap.payload
        tab['OBJID'][:] = 'None'
        tab['OBJID'][obs_ap[i]] = objid[obs_obj[i]]
        tab['RA'][obs_ap[i]] = ra[obs_obj[i]]
        tab['DEC'][obs_ap[i]] = dec[obs_obj[i]]
        if tight:
            indx = tab['OBJID'] != 'None'
            tab = tab[indx]

        print(f'Writing: {ofile}')
        with open(ofile, 'w') as f:
            for l in header_text.split('\n'):
                f.write(f'# {l}\n')
            tab.write(f, format='ascii.fixed_width_no_header', delimiter=' ')


def empty_configuration_table(nrows=None, objidlen=20):
    """
    Construct an empty FOBOS configuration table.

    Args:
        nrows (:obj:`int`, optional):
            Number of table rows.  If None, initialized to 0.

    Returns:
        `astropy.table.Table`_: Empty configuration table.
    """
    length = 0 if nrows is None else nrows
    _objidlen = max(5, objidlen)
    return table.Table([table.Column(name='MID', dtype=int, length=length,
                                     description='Module ID Number', format='>3'),
                        table.Column(name='SID', dtype=int, length=length,
                                     description='Spectrograph ID Number', format='>3'),
                        table.Column(name='BID', dtype=int, length=length,
                                     description='Starbug ID Number', format='>3'),
                        table.Column(name='Type', dtype=int, length=length,
                                     description='Starbug payload type', format='>4'),
                        table.Column(name='OBJID', dtype=f'<U{_objidlen}', length=length,
                                     description='Object identifier', format=f'>{_objidlen}'),
                        table.Column(name='RA', dtype=float, length=length,
                                     description='Object right ascension', format='12.8f'),
                        table.Column(name='DEC', dtype=float, length=length,
                                     description='Object declination', format='11.8f')
                       ])


def parse_configuration(config_file, parse_source=False):
    """
    Parse the FOBOS configuration from a configuration file.

    Parameters
    ----------
    config_file : :obj:`str`
        Name of the configuration file.
    parse_source : :obj:`bool`, optional
        If possible given the data within the configuration, re-read and return
        the full target list from the target source file, instead of only those
        with allocated apertures.
    
    Returns
    -------
    objid : `numpy.ndarray`_
        1d array with the target ID numbers
    ra : `numpy.ndarray`_
        Target right ascension in decimal degrees
    dec : `numpy.ndarray`_
        Target declination in decimal degrees
    center : :obj:`tuple`
        Two-tuple with the coordinates of the pointing center
    assigned_obj : `numpy.ndarray`_
        Integer array with the array index in, e.g., ``objid`` with the target
        that has been assigned a FOBOS aperture.
    ap : :class:`~producer.deploy.FOBOSApertures`
        Aperture deployment instance.
    assigned_ap : `numpy.ndarray`_
        Integer array with the indices of the apertures in ``ap`` that have been
        assigned to a target.
    """

    _cfg = Path(config_file).resolve()
    if not _cfg.is_file():
        raise FileNotFoundError(f'{str(_cfg)} does not exist!')

    # Read the data table
    db = numpy.genfromtxt(str(_cfg), dtype=str)
    mid, sid, bid, payload = db.T[:4].astype(int)
    objid = db[:,4]
    ra, dec = db.T[5:].astype(float)

    # All of this is to just parse the configuration...
    with open(str(_cfg), 'r') as f:
        lines = numpy.asarray(f.readlines())
    lines = lines[numpy.array([l[0] == '#' for l in lines])]
    lines = [l[2:].strip() for l in lines]
    for i in range(len(lines)):
        if len(lines[i]) > 0 and lines[i] == '-'*len(lines[i]):
            i += 1
            break
    for j in range(len(lines)-1, -1, -1):
        if len(lines[j]) > 0 and lines[j] == '-'*len(lines[j]):
            break

    # Read the configuration parameters
    cfg = ConfigParser()
    cfg.read_string('\n'.join(lines[i:j]))

    # Instantiate the aperture deployment object
    if cfg['FOBOS_MODE'].get('version') != FOBOSApertures.version:
        raise ValueError('FOBOS aperture-deployment version has changed!')
    mode = numpy.array(cfg['FOBOS_MODE'].get('mode').split(',')).astype(int)
    baseline = cfg['FOBOS_MODE'].getboolean('baseline')
    design = cfg['FOBOS_MODE'].getint('design')
    ap = FOBOSApertures(mode=mode, baseline=baseline, config=design)

    # Attempt to parse original target file?
    target_file = None
    if parse_source:
        target_file = Path(cfg['TARGETS'].get('source')).resolve()
        if not target_file.is_file():
            warnings.warn(f'Original target source file does not exist: {str(target_file)}')
            target_file = None
        ra_c = cfg['TARGETS'].get('ra_col')
        dec_c = cfg['TARGETS'].get('dec_col')
        if 'None' in [ra_c, dec_c]:
            warnings.warn('RA and/or DEC columns not defined.  Not re-reading source catalog.')
            target_file = None
        else:
            ra_c = int(ra_c)
            dec_c = int(dec_c)

    # Get the pointing center
    center = (cfg['POINTING'].getfloat('ra'), cfg['POINTING'].getfloat('dec'))

    # Find the assigned apertures indices in the 
    # TODO: There needs to be a better way to do this...
    assigned_ap = numpy.array([numpy.where((m == ap.mid) & (b == ap.bid))[0][0]
                                for m, b in zip(mid, bid)])

    if target_file is None:
        # Only returning the allocated list, so we're done
        return objid.astype(int), ra, dec, center, numpy.arange(assigned_ap.size), ap, assigned_ap

    if cfg['TARGETS'].get('id_type') != 'index':
        raise NotImplementedError('Can only handle object ID that are the index of the target '
                                  'in the source file.')

    # Re-read the target catalog
    full_ra, full_dec, _ = parse_targets(target_file, ra_c=ra_c, dec_c=dec_c)

    return numpy.arange(full_ra.size), full_ra, full_dec, center, objid.astype(int), \
                ap, assigned_ap




