"""
Hard-coded method to construct the FOBOS focal plane design files.

.. include:: ../include/links.rst
"""

import time

from IPython import embed

import numpy

from .focalplane import module_layout
from ..deploy import FOBOSModeBitMask


def baseline_module_type():
    """
    Define the type of each module.

    Module types are:
       0. single starbugs only
       1. mini-IFU module
       2. location of Big IFU
       3. imaging bundle module
       4. flux-calibration module

    """
    return numpy.array([0,   # 1
                        0,   # 2
                        0,   # 3
                        0,   # 4
                        3,   # 5
                        0,   # 6
                        0,   # 7
                        3,   # 8
                        0,   # 9
                        0,   # 10
                        3,   # 11
                        1,   # 12
                        1,   # 13
                        3,   # 14
                        1,   # 15
                        1,   # 16
                        3,   # 17
                        0,   # 18
                        0,   # 19
                        0,   # 20
                        1,   # 21
                        4,   # 22
                        1,   # 23
                        1,   # 24
                        4,   # 25
                        1,   # 26
                        0,   # 27
                        0,   # 28
                        0,   # 29
                        4,   # 30
                        3,   # 31
                        1,   # 32
                        1,   # 33
                        3,   # 34
                        1,   # 35
                        1,   # 36
                        3,   # 37
                        4,   # 38
                        0,   # 39
                        0,   # 40
                        1,   # 41
                        1,   # 42
                        4,   # 43
                        1,   # 44
                        1,   # 45
                        3,   # 46
                        1,   # 47
                        1,   # 48
                        0,   # 49
                        0,   # 50
                        1,   # 51
                        4,   # 52
                        1,   # 53
                        1,   # 54
                        0,   # 55
                        1,   # 56
                        1,   # 57
                        4,   # 58
                        1,   # 59
                        2,   # 60
                        0,   # 61
                        1,   # 62
                        1,   # 63
                        3,   # 64
                        1,   # 65
                        1,   # 66
                        4,   # 67
                        1,   # 68
                        1,   # 69
                        0,   # 70
                        0,   # 71
                        4,   # 72
                        3,   # 73
                        1,   # 74
                        1,   # 75
                        3,   # 76
                        1,   # 77
                        1,   # 78
                        3,   # 79
                        4,   # 80
                        0,   # 81
                        0,   # 82
                        0,   # 83
                        1,   # 84
                        4,   # 85
                        1,   # 86
                        1,   # 87
                        4,   # 88
                        1,   # 89
                        0,   # 90
                        0,   # 91
                        0,   # 92
                        3,   # 93
                        1,   # 94
                        1,   # 95
                        3,   # 96
                        1,   # 97
                        1,   # 98
                        3,   # 99
                        0,   # 100
                        0,   # 101
                        3,   # 102
                        0,   # 103
                        0,   # 104
                        3,   # 105
                        0,   # 106
                        0,   # 107
                        0,   # 108
                        0])  # 109


def module_mode(mod_type):
    """
    Provide the default mode and designated sky bits of the starbugs in each
    module type.

    There are 5 module types:

        0. single starbugs only
        1. mini-IFU module
        2. monolithic IFU
        3. imaging bundle module
        4. flux-calibration module
    
    The active starbugs depend on the spectrograph mode.  The modes are
    tracked using a bit value, where bits can be either:

        0. Single fiber mode.
        1. multi-IFU mode.
        2. Monolithic IFU mode.

    So, e.g., a starbug with a mode value of 7 is deployed for all 3 modes.
    Starbug ID (or array index + 1) per module is expected to be.

    |               03
    |           02      07
    |       01      06      12
    |           05      11
    |       04      10      16
    |           09      15
    |       08      14      19
    |           13      18
    |               17

    Only IFU modules have designated sky apertures, which are starbug IDs 3, 4,
    16, and 17.

    Args:
        mod_type (:obj:`int`):
            The module type.  Must be 0, 1, 2, 3, or 4, as listed above.

    Returns:
        :obj:`tuple`: Two integer `numpy.ndarray`_ objects providing (1) bits
        used to determine if the starbug is active during the spectrograph mode
        and (2) bits used to determine if the starbug houses a designated sky
        fiber during the spectrograph mode.
    """
    bm = FOBOSModeBitMask()
    active_bits = numpy.zeros(19, dtype=numpy.int16)
    sky_bits = numpy.zeros(19, dtype=numpy.int16)
    if mod_type == 0:
        # Single-fiber only module: All starbugs are active in MOS mode only.
        # None of the apertures are explicitly designated for sky observations.
        return bm.turn_on(active_bits, 'MOS'), sky_bits
    if mod_type == 1:
        # IFU module: All but the center starbug are active in MOS mode, none of
        # which are explicitly designated for sky.
        active_bits = bm.turn_on(active_bits, 'MOS')
        active_bits[9] = bm.turn_off(active_bits[9], 'MOS')
        # In IFU mode, the center starbug is active, and 4 single-fiber
        # apertures on the outer ring are designated for sky-only observations.
        active_bits[[2,3,9,15,16]] = bm.turn_on(active_bits[[2,3,9,15,16]], 'IFU')
        sky_bits[[2,3,15,16]] = bm.turn_on(sky_bits[[2,3,15,16]], 'IFU')
        return active_bits, sky_bits
    if mod_type == 2:
        # Monolithic IFU module: All starbugs are turn off.
        return active_bits, sky_bits
    if mod_type in [3, 4]:
        # Guide and flux-calibration modules: The center starbug is always on,
        # but the other bugs are only on in MOS mode.  None of the apertures are
        # explicitly designated for sky.
        active_bits = bm.turn_on(active_bits, 'MOS')
        active_bits[9] = bm.turn_on(active_bits[9], 'IFU')
        active_bits[9] = bm.turn_on(active_bits[9], 'MONO')
        return active_bits, sky_bits

    raise ValueError(f'Unknown module type: {mod_type}')


def _special_module_mode_config1(active_bits, sky_bits, module_id):
    bm = FOBOSModeBitMask()
    _active_bits = active_bits.copy()
    _sky_bits = sky_bits.copy()
    if module_id not in [38, 44, 49, 58, 59, 70, 80]:
        # No changes needed
        return _active_bits, _sky_bits
    if module_id == 44:
        # Turn off the current MOS mode bits
        _active_bits = bm.turn_off(_active_bits, 'MOS')
        # The remaining starbugs with non-zero bits should all be for IFU mode;
        # include these in MOS mode, as well, for both the active bits and the
        # designated sky fibers.
        indx = _active_bits > 0
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MOS')
        indx = _sky_bits > 0
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MOS')
        # The fiber budget allows us to turn on 5 more single fibers for MOS
        # mode only
        _active_bits[[0,4,10,11,13]] = bm.turn_on(_active_bits[[0,4,10,11,13]], 'MOS')
        return _active_bits, _sky_bits
    if module_id in [49, 59, 70]:
        # Turn on the outer ring of fibers for use as sky fibers during
        # the monolithic IFU mode
        indx = [0,1,2,3,6,7,8,10,11,12,15,16,17,18]
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MONO')
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MONO')
        return _active_bits, _sky_bits
    if module_id in [38, 58, 80]:
        # Turn on the outer ring of fibers for use as sky fibers during
        # the monolithic IFU mode
        indx = [0,1,2,3,6,7,8,11,12,15,16,17,18]
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MONO')
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MONO')
        return _active_bits, _sky_bits
        
    
def _special_module_mode_config2(active_bits, sky_bits, module_id):
    bm = FOBOSModeBitMask()
    _active_bits = active_bits.copy()
    _sky_bits = sky_bits.copy()
    if module_id not in [14, 15, 44, 54, 65, 70, 80]:
        # No changes needed
        return _active_bits, _sky_bits
    if module_id == 44:
        # Turn off the current MOS mode bits
        _active_bits = bm.turn_off(_active_bits, 'MOS')
        # The remaining starbugs with non-zero bits should all be for IFU mode;
        # include these in MOS mode, as well, for both the active bits and the
        # designated sky fibers.
        indx = _active_bits > 0
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MOS')
        indx = _sky_bits > 0
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MOS')
        # The fiber budget allows us to turn on 5 more single fibers for MOS
        # mode only
        _active_bits[[0,4,10,11,13]] = bm.turn_on(_active_bits[[0,4,10,11,13]], 'MOS')
        return _active_bits, _sky_bits
    if module_id in [14, 65, 70]:
        # Turn on the outer ring of fibers for use as sky fibers during
        # the monolithic IFU mode
        indx = [0,1,2,3,6,7,8,10,11,12,15,16,17,18]
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MONO')
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MONO')
        return _active_bits, _sky_bits
    if module_id in [15, 54, 80]:
        # Turn on the outer ring of fibers for use as sky fibers during
        # the monolithic IFU mode
        indx = [0,1,2,3,6,7,8,11,12,15,16,17,18]
        _active_bits[indx] = bm.turn_on(_active_bits[indx], 'MONO')
        _sky_bits[indx] = bm.turn_on(_sky_bits[indx], 'MONO')
        return _active_bits, _sky_bits
        
    
def special_module_mode(active_bits, sky_bits, module_id, spec_map):
    """
    Apply adjustments to active and sky bits for some modules:

        - Spectrograph 1 always deploys a single IFU for quick ToO follow-up.
          This means that module 44 always deploys an IFU and 4 sky fibers, and
          the other starbugs never do anything...

        - When using the monolithic IFU, the 12 starbugs in the outer hexagonal
          ring are activated for 2 modules in each spectrograph for sky-only
          observations.  The modules are 38, 49, 58, 59, 70, and 80 in
          spectrograph-to-module-mapping design 1, and 14, 15, 54, 65, 70, 80
          for design 2.

    """
    if spec_map == 1:
        return _special_module_mode_config1(active_bits, sky_bits, module_id)
    if spec_map == 2:
        return _special_module_mode_config2(active_bits, sky_bits, module_id)
    raise ValueError(f'Unknown spectrograph configuration {spec_map}.')


def module_spectrograph():
    """
    Define two possible module-to-spectrograph allocations.

    The first distributes the modules (roughly) evenly between the 3
    spectrographs across the full field of view.  The second divides the
    focal plane into three regions, one per spectrograph.

    The returned index gives the spectrograph number.
    """
    return numpy.array([[2,0],   # 1
                        [1,0],   # 2
                        [0,2],   # 3
                        [0,0],   # 4
                        [1,0],   # 5
                        [2,0],   # 6
                        [1,2],   # 7
                        [0,2],   # 8
                        [0,2],   # 9
                        [2,0],   # 10
                        [0,0],   # 11
                        [1,0],   # 12
                        [0,0],   # 13
                        [2,0],   # 14
                        [1,2],   # 15
                        [0,2],   # 16
                        [1,2],   # 17
                        [2,2],   # 18
                        [1,0],   # 19
                        [1,0],   # 20
                        [0,0],   # 21
                        [1,0],   # 22
                        [2,0],   # 23
                        [0,2],   # 24
                        [2,2],   # 25
                        [2,2],   # 26
                        [2,2],   # 27
                        [1,2],   # 28
                        [0,0],   # 29
                        [0,0],   # 30
                        [2,0],   # 31
                        [2,0],   # 32
                        [1,0],   # 33
                        [0,2],   # 34
                        [2,2],   # 35
                        [1,2],   # 36
                        [2,2],   # 37
                        [0,2],   # 38
                        [0,2],   # 39
                        [2,0],   # 40
                        [0,0],   # 41
                        [1,0],   # 42
                        [2,0],   # 43
                        [0,0],   # 44
                        [1,2],   # 45
                        [1,2],   # 46
                        [0,2],   # 47
                        [1,2],   # 48
                        [1,2],   # 49
                        [0,0],   # 50
                        [2,0],   # 51
                        [1,0],   # 52
                        [0,0],   # 53
                        [2,0],   # 54
                        [0,2],   # 55
                        [0,2],   # 56
                        [2,2],   # 57
                        [1,2],   # 58
                        [0,2],   # 59
                        [-1,-1],   # 60
                        [1,0],   # 61
                        [1,0],   # 62
                        [2,0],   # 63
                        [1,0],   # 64
                        [1,1],   # 65
                        [2,1],   # 66
                        [0,2],   # 67
                        [1,2],   # 68
                        [2,2],   # 69
                        [2,2],   # 70
                        [0,0],   # 71
                        [2,1],   # 72
                        [0,0],   # 73
                        [1,1],   # 74
                        [0,1],   # 75
                        [2,1],   # 76
                        [1,1],   # 77
                        [0,1],   # 78
                        [0,2],   # 79
                        [2,1],   # 80
                        [2,2],   # 81
                        [1,1],   # 82
                        [2,1],   # 83
                        [0,1],   # 84
                        [0,1],   # 85
                        [2,1],   # 86
                        [0,1],   # 87
                        [1,1],   # 88
                        [2,1],   # 89
                        [1,1],   # 90
                        [1,2],   # 91
                        [2,1],   # 92
                        [1,1],   # 93
                        [2,1],   # 94
                        [1,1],   # 95
                        [0,1],   # 96
                        [2,1],   # 97
                        [1,1],   # 98
                        [2,1],   # 99
                        [0,1],   # 100
                        [0,1],   # 101
                        [2,1],   # 102
                        [1,1],   # 103
                        [2,1],   # 104
                        [1,1],   # 105
                        [0,1],   # 106
                        [0,1],   # 107
                        [1,1],   # 108
                        [2,1]]).T  # 109


def fobos_module_layout(module_file, starbug_file, version='0.1'):

    starbug_grid, module_grid = module_layout(0.725, 19., 5, 13, module_pitch=90., max_dist=11.2)
    if starbug_grid.shape[0] != 19:
        raise ValueError('FOBOS modules should have 19 Starbugs.')
    if module_grid.shape[0] != 109:
        raise ValueError('Number of FOBOS modules changed.')

    module_dist = numpy.sqrt(module_grid[:,0]**2 + module_grid[:,1]**2)
    module_dist /= 60.

    # Write the file with the module coordinates
    module_type = baseline_module_type()
    module_spec = module_spectrograph()
    module_base = module_dist < 10.4
    #   0 - single starbugs only
    #   1 - Location of Big IFU
    #   2 - mini-IFU module
    #   3 - imaging bundle module
    #   4 - flux-calibration module
    header_text = f'FOBOS Starbug module layout (v{version})\n\n' \
                  f'Generated {time.strftime("%a %d %b %Y %H:%M:%S", time.localtime())}\n\n' \
                  'Columns are:\n' \
                  '    MID - Module ID number\n' \
                  '    XI, ETA - Module center coordinates in arcmin\n' \
                  '              relative to the focal plan center\n' \
                  '    T - The special payload type at the module\n' \
                  '        center: (0) single-fiber, (1) 37-fiber IFU, \n' \
                  '        (2) monolithic IFU, (3) imaging bundle,\n' \
                  '        (4) flux-calibration bundle.\n' \
                  '    B - Flag (1=True) that the module is in the\n' \
                  '        baseline CoDR design.\n\n'
    # Save the module positions and types
    numpy.savetxt(module_file,
                  numpy.column_stack((numpy.arange(module_grid.shape[0])+1,
                                      module_grid[:,0]/60., module_grid[:,1]/60.,
                                      module_type, module_base.astype(int))),
                  fmt='%5d %7.3f %7.3f %2d %2d',
                  header=header_text+f"{'MID':>3} {'XI':>7} {'ETA':>7} {'T':>2} {'B':>2}")

    # Build the coordinates, types, mode, etc, for all Starbugs
    module_id = numpy.tile(numpy.arange(module_grid.shape[0])+1, (starbug_grid.shape[0],1)).T
    active_bit_ran = numpy.zeros(module_id.shape, dtype=numpy.int16)
    active_bit_reg = numpy.zeros(module_id.shape, dtype=numpy.int16)
    sky_bit_ran = numpy.zeros(module_id.shape, dtype=numpy.int16)
    sky_bit_reg = numpy.zeros(module_id.shape, dtype=numpy.int16)
    for i in module_id[:,0]:
        active_bit_ran[i-1,:], sky_bit_ran[i-1,:] \
                = special_module_mode(*module_mode(module_type[i-1]), i, 1)
        active_bit_reg[i-1,:], sky_bit_reg[i-1,:] \
                = special_module_mode(*module_mode(module_type[i-1]), i, 2)
    starbug_id = numpy.tile(numpy.arange(starbug_grid.shape[0])+1, (module_grid.shape[0],1))
    starbug_type = numpy.zeros(module_id.shape, dtype=int)
    # "Special" starbug always at the module center
    starbug_type[:,9] = module_type
    starbug_spec_ran = numpy.tile(module_spec[0], (starbug_grid.shape[0],1)).T+1
    starbug_spec_reg = numpy.tile(module_spec[1], (starbug_grid.shape[0],1)).T+1
    starbug_base = numpy.tile(module_base, (starbug_grid.shape[0],1)).T

    starbug_coo = (module_grid[:,None,:] + starbug_grid[None,:,:])
    # Remove the module with the monolithic IFU
    indx = module_type != 2
    module_id = module_id[indx].ravel()
    starbug_id = starbug_id[indx].ravel()
    starbug_type = starbug_type[indx].ravel()
    active_bit_ran = active_bit_ran[indx].ravel()
    active_bit_reg = active_bit_reg[indx].ravel()
    sky_bit_ran = sky_bit_ran[indx].ravel()
    sky_bit_reg = sky_bit_reg[indx].ravel()
    starbug_spec_ran = starbug_spec_ran[indx].ravel()
    starbug_spec_reg = starbug_spec_reg[indx].ravel()
    starbug_base = starbug_base[indx].ravel()
    starbug_coo = starbug_coo[indx].reshape(-1,2)

    # Write the file with the Starbug coordinates
    header_text = f'FOBOS Starbug module layout (v{version})\n\n' \
                  f'Generated {time.strftime("%a %d %b %Y %H:%M:%S", time.localtime())}\n\n' \
                  'Provides bits to select which Starbugs are active (ON) and which \n' \
                  'are explicitly designated for sky observations (SKY) in each of 3 \n' \
                  'modes: (0) single-fiber aperture deployment, (1) 37-fiber IFU \n' \
                  'deployment, (2) monolithic IFU deployment.  These bits are based on \n' \
                  'two designs for the spectrograph-to-module mapping: modules for each \n' \
                  'spectrograph are distributed over (1) the full field-of-view or (2) a \n' \
                  'designated region covering 1/3 of the field-of-view.  For example, to \n' \
                  'select active starbugs in mode 0 in design 1, do the bitwise \n' \
                  'operation:\n\n' \
                  '    ON1 & (1 << 0) \n\n' \
                  'Columns are:\n' \
                  '    BID - Starbug ID number\n' \
                  '    MID - Module ID number\n' \
                  '    SPC1 - Spectrograph (design 1)\n' \
                  '    ON1 - Selection bit for active starbugs (design 1)\n' \
                  '    SKY1 - Selection bit for designated sky apertures (design 1)\n' \
                  '    SPC2 - Spectrograph (design 2)\n' \
                  '    ON2 - Selection bit for active starbugs (design 2)\n' \
                  '    SKY2 - Selection bit for designated sky apertures (design 2)\n' \
                  '    XI, ETA - Starbug "home" coordinates in arcmin\n' \
                  '              relative to the focal plane center\n' \
                  '    T - Payload type: (0) single-fiber, (1) 37-fiber \n' \
                  '        IFU, (2) monolithic IFU, (3) imaging bundle,\n' \
                  '        (4) flux-calibration bundle.\n' \
                  '    B - Flag (1=True) that the module/starbug is in the\n' \
                  '        baseline CoDR design.\n\n'
    # Save the starbug positions and types
    numpy.savetxt(starbug_file,
                  numpy.column_stack((starbug_id, module_id, starbug_spec_ran,
                                      active_bit_ran, sky_bit_ran, starbug_spec_reg,
                                      active_bit_reg, sky_bit_reg, starbug_coo[:,0]/60.,
                                      starbug_coo[:,1]/60., starbug_type,
                                      starbug_base.astype(int))),
                  fmt='%5d %4d %4d %4d %4d %4d %4d %4d %8.4f %8.4f %2d %2d',
                  header=header_text+f"{'BID':>3} {'MID':>4} {'SPC1':>4} {'ON1':>4} "
                                      f"{'SKY1':>4} {'SPC2':>4} {'ON2':>4} {'SKY2':>4} "
                                      f"{'XI':>8} {'ETA':>8} {'T':>2} {'B':>2}")



