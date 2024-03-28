"""
Provides class and functions to interact with FOBOS aperture deployment.

.. include:: ../include/links.rst
"""

import time
from itertools import chain, combinations

from IPython import embed 

import numpy
from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot, patches

from sklearn.neighbors import KDTree

from astropy import units
from astropy.coordinates import SkyCoord, SkyOffsetFrame

from . import data_file
from .collisions import remove_collisions
from .data.bitmask import BitMask


class FOBOSModeBitMask(BitMask):
    def __init__(self):
        mask_bits = dict(MOS='Single-fiber, multi-object spectroscopy mode',
                         IFU='Multiplexed IFU mode',
                         MONO='Monolithic IFU mode')
        super().__init__(list(mask_bits.keys()), descr=list(mask_bits.values()))

    def validate(self, mode):
        _mode = numpy.atleast_1d(mode)
        valid_modes = [0] + (numpy.array(list(self.bits.values()))+1).tolist()
        if any([m not in valid_modes for m in _mode]):
            raise ValueError('Modes invalide.  Must be 0 (off), 1 (single-fiber), 2 (multi-IFU), '
                             'or 3 (monolithic IFU).')


class FOBOSApertureBitMask(BitMask):
    def __init__(self):
        mask_bits = dict(SCIENCE='Science aperture',
                         SKY='Designated sky fiber',
                         GUIDE='Tracking and focus-offset imaging bundle',
                         CALIB='Flux-calibration bundle')
        super().__init__(list(mask_bits.keys()), descr=list(mask_bits.values()))


class FOBOSApertures:

    module_src = data_file('deploy/fobos_modules.db')
    """
    File providing the centers of all FOBOS starbug modules.  This is currently
    never used.
    """

    starbug_src = data_file('deploy/fobos_starbugs.db')
    """
    Main file with the distribution, types, etc of the FOBOS apertures.
    """

    version = '0.1'
    """
    Version of the aperture deployment design.
    """

    mode_bm = FOBOSModeBitMask()
    """
    Mode selection bit mask
    """

    ap_bm = FOBOSApertureBitMask()
    """
    Aperture selection bit mask
    """

    fov = 20./60.
    """
    Diameter of the field-of-view in decimal degrees.
    """

    def __init__(self, mode=1, baseline=True, config=1):
        """
        Args:
            mode (:obj:`int`, array-like, optional):
                The mode assignment for each spectrograph.  Must be 1 for the
                single-fiber apertures, 2 for multi-IFU mode, or 3 for
                monolithic IFU mode.  If a single integer, all spectrographs are
                put in the same mode; otherwise, must provide the mode for each
                of the 3 spectrographs separately.
            baseline (:obj:`bool`, optional):
                Flag to only return the apertures in the CoDR baseline
                configuration.
            config (:obj:`int`, optional):
                Configuration selection for mapping aperture modules to each
                spectrograph.  Must be either: (1) modules mapped to each
                spectrograph are spread over the entire focal plane; or (2)
                modules mapped to each spectrograph fill coherent focal-plane
                zones.
        """
        # TODO: Check version?
        _data = numpy.genfromtxt(str(self.starbug_src), dtype=str)
        self.bid = _data[:,0].astype(int)
        self.mid = _data[:,1].astype(int)
        self.nap = self.bid.size

        self.id = numpy.arange(self.nap, dtype=int)
        self.id_name = numpy.array([f'{m}-{b}' for m, b in zip(self.mid, self.bid)])

        self.spc = None
        self.on = None
        self.sky = None

        self._spc1 = _data[:,2].astype(int)
        self._on1 = _data[:,3].astype(int)
        self._sky1 = _data[:,4].astype(int)

        self._spc2 = _data[:,5].astype(int)
        self._on2 = _data[:,6].astype(int)
        self._sky2 = _data[:,7].astype(int)

        self.coo = _data[:,8:10].astype(float)
        self.in_baseline = _data[:,11].astype(int).astype(bool)

        self.payload = _data[:,10].astype(int)

        # Initialize the aperture types based only on the payload type.
        self.type = numpy.zeros(self.nap, dtype=self.ap_bm.minimum_dtype(asuint=True))
        indx = numpy.isin(self.payload, [0,1,2])
        self.type[indx] = self.ap_bm.turn_on(self.type[indx], 'SCIENCE')
        indx = self.payload == 3
        self.type[indx] = self.ap_bm.turn_on(self.type[indx], 'GUIDE')
        indx = self.payload == 4
        self.type[indx] = self.ap_bm.turn_on(self.type[indx], 'CALIB')

        self.mode = None
        self.baseline = None
        self.config = None
        self.active = None
        self.configure(mode=mode, baseline=baseline, config=config)

    @staticmethod
    def parse_mode(mode):
        """
        Parse and validate the provided FOBOS spectrograph mode.

        Args:
            mode (:obj:`int`, array-like):
                The mode assignment for each spectrograph.  Must be 0 to turn
                off all apertures, 1 to select single-fiber mode, 2 for
                multi-IFU mode, or 3 for monolithic IFU mode.  If a single
                integer, all spectrographs are put in the same mode; otherwise,
                must provide the mode for each of the 3 spectrographs
                separately.
    
        Returns:
            `numpy.ndarray`_: A 3-element array with the integer mode identifier
            for each of the three spectrographs.

        Raises:
            ValueError:
                Raised if the mode number is not understood, or if the mode is
                not provided as a single value or one for each spectrograph.
        """
        # Check the mode input
        _mode = numpy.atleast_1d(mode).astype(int)
        if _mode.size == 1:
            _mode = numpy.repeat(_mode, 3)
        if _mode.size != 3:
            raise ValueError(f'Mode not understood: {mode}.  Must enter single mode for all '
                            f'spectrographs, or the mode for each of the three spectrographs.')
        # Check the modes are valid        
        FOBOSApertures.mode_bm.validate(_mode)
        return _mode

    def configure(self, mode=None, baseline=None, config=None):
        """
        Set the active apertures based on the provided mode.

        Configuration sets internal attributes that can be accessed

        Args:
            mode (:obj:`int`, array-like):
                The mode assignment for each spectrograph.  Must be 1 for the
                single-fiber apertures, 2 for multi-IFU mode, or 3 for
                monolithic IFU mode.  If a single integer, all spectrographs are
                put in the same mode; otherwise, must provide the mode for each
                of the 3 spectrographs separately.
            baseline (:obj:`bool`, optional):
                Flag to only return the apertures in the CoDR baseline
                configuration.
            config (:obj:`int`, optional):
                Configuration selection for mapping aperture modules to each
                spectrograph.  Must be either: (1) modules mapped to each
                spectrograph are spread over the entire focal plane; or (2)
                modules mapped to each spectrograph fill coherent focal-plane
                zones.
        """
        if mode is None and self.mode is None:
            # Assume this is the first configuration
            self.mode = self.parse_mode(1)
        elif mode is not None:
            # Parse and check the mode
            self.mode = self.parse_mode(mode)

        if baseline is None and self.baseline is None:
            # Assume this is the first configuration
            self.baseline = True
        elif baseline is not None:
            self.baseline = baseline

        if config is None and self.config is None:
            # Assume this is the first configuration
            self.config = 1
        elif config is not None:
            self.config = config

        # Set the configuration
        if self.config == 1:
            self.spc = self._spc1
            self.on = self._on1
            self.sky = self._sky1
        elif self.config == 2:
            self.spc = self._spc2
            self.on = self._on2
            self.sky = self._sky2
        else:
            raise ValueError(f'Spectrograph configuration {config} unknown.  Must be 1 or 2.')

        # Reset all science and sky apertures
        indx = numpy.isin(self.payload, [0,1,2])
        self.type[indx] = self.ap_bm.turn_on(self.type[indx], 'SCIENCE')
        self.type = self.ap_bm.turn_off(self.type, 'SKY')

        # Loop through each spectrograph to find the active apertures and set
        # their type.
        self.active = numpy.zeros(self.nap, dtype=bool)
        for i in range(3):
            if self.mode[i] == 0:
                continue
            _mode = self.mode_bm.keys()[self.mode[i]-1]
            self.active |= (self.spc == i+1) & self.mode_bm.flagged(self.on, _mode)
            # Set sky fibers (if any)
            indx = self.mode_bm.flagged(self.sky, _mode)
            if numpy.any(indx):
                self.type[indx] = self.ap_bm.turn_off(self.type[indx], 'SCIENCE')
                self.type[indx] = self.ap_bm.turn_on(self.type[indx], 'SKY')

        # Restrict based on baseline
        if self.baseline:
            self.active &= self.in_baseline

    def select(self, payload):
        r"""
        Construct a boolean array that selects specific payload types.

        Type must be:
            - ``'science'`` for science apertures,
            - ``'sky'`` for designated sky apertures,
            - ``'guide'`` for imaging guide bundles, or
            - ``'calib'`` for flux-calibration bundles.

        Note that in MOS mode (mode=1), the only sky fibers are the ones
        designated for the always-ready IFU in spectrograph 1.

        Args:
            payload (:obj:`str`):
                The payload selection.
            return_id (:obj:`bool`, optional):
                Return the running (0-indexed) number for the aperture.

        Returns:
            `numpy.ndarray`_: A boolean array used to select the viable
            apertures of the requested type.
        """
        valid = list(self.ap_bm.bits.keys())
        if payload.upper() not in valid:
            raise ValueError(f"Unknown payload type: {payload.upper()}.  Options are: "
                             f"{', '.join(valid)}")
        return self.active & self.ap_bm.flagged(self.type, payload.upper())

    def show(self, include_patrol=False, by_spec=False, legend=True, ax=None, show=True,
             s=40, show_all=False):
        """
        Make a plot of the currently active apertures.

        Args:
            include_patrol (:obj:`bool`, optional):
                Use a light shaded region to show the patrol region of every
                aperture.
            by_spec (:obj:`bool`, optional):
                Instead of coloring the aperture locations by their payload
                type, color them by their spectrograph mapping.
            legend (:obj:`bool`, optional):
                Include the plot legend
        """

        if ax is None:
            w,h = pyplot.figaspect(1)
            fig = pyplot.figure(figsize=(1.5*w,1.5*h))
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            ax.set_xlim([-14, 14])
            ax.set_ylim([-14, 14])
            ax.minorticks_on()
            ax.grid(True, which='major', color='0.9', zorder=0, linestyle='-')
            ax.tick_params(which='major', direction='in', length=8, top=True, right=True)
            ax.tick_params(which='minor', direction='in', length=4, top=True, right=True)
            ax.set_xlabel(r'$\xi$ [arcmin]')
            ax.set_ylabel(r'$\eta$ [arcmin]')

        # Single fibers
        indx = self.payload == 0
        if self.baseline:
            indx &= self.in_baseline
        if not show_all:
            indx &= self.active
        if by_spec:
            _indx = indx & (self.spc == 1)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='.', s=s, lw=0, color='C0', zorder=5, label='Single Fiber (Spec 1)')
            _indx = indx & (self.spc == 2)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='.', s=s, lw=0, color='C2', zorder=5, label='Single Fiber (Spec 2)')
            _indx = indx & (self.spc == 3)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='.', s=s, lw=0, color='C4', zorder=5, label='Single Fiber (Spec 3)')
        else:
            ax.scatter(self.coo[indx,0], self.coo[indx,1],
                    marker='.', s=s, lw=0, color='C0', zorder=5, label='Single Fiber')
        # IFUs
        indx = self.payload == 1
        if self.baseline:
            indx &= self.in_baseline
        if not show_all:
            indx &= self.active
        if by_spec:
            _indx = indx & (self.spc == 1)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='s', s=s//2, lw=0.5, color='C0',zorder=5, label='IFU (Spec 1)')
            _indx = indx & (self.spc == 2)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='s', s=s//2, lw=0.5, color='C2',zorder=5, label='IFU (Spec 2)')
            _indx = indx & (self.spc == 3)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='s', s=s//2, lw=0.5, color='C4',zorder=5, label='IFU (Spec 3)')
        else:
            ax.scatter(self.coo[indx,0], self.coo[indx,1],
                       marker='s', s=s//2, lw=0.5, color='C1',zorder=5, label='IFU')

        # Guide bundles
        indx = self.payload == 3
        if self.baseline:
            indx &= self.in_baseline
        if not show_all:
            indx &= self.active
        ax.scatter(self.coo[indx,0], self.coo[indx,1],
                    marker='o', s=s//2, lw=0.5, color='C3', zorder=5, label='Guide Bundle')

        # Flux-calibration bundles
        indx = self.payload == 4
        if self.baseline:
            indx &= self.in_baseline
        if not show_all:
            indx &= self.active
        if by_spec:
            _indx = indx & (self.spc == 1)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='x', s=s//2, lw=1, color='C0', zorder=5,
                       label='Calib. Bundle (Spec 1)')
            _indx = indx & (self.spc == 2)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='x', s=s//2, lw=1, color='C2', zorder=5,
                       label='Calib. Bundle (Spec 1)')
            _indx = indx & (self.spc == 3)
            ax.scatter(self.coo[_indx,0], self.coo[_indx,1],
                       marker='x', s=s//2, lw=1, color='C4', zorder=5,
                       label='Calib. Bundle (Spec 1)')
        else:
            ax.scatter(self.coo[indx,0], self.coo[indx,1],
                       marker='x', s=s//2, lw=1, color='C4', zorder=5,
                       label='Calib. Bundle')

        # Include the patrol region (this in combination with the legend can
        # make the plot very slow)
        if include_patrol:
            indx = numpy.isin(self.payload, [0,1,3,4])
            if self.baseline:
                indx &= self.in_baseline
            if not show_all:
                indx &= self.active
            for x, y in self.coo[indx]:
                ax.add_patch(patches.Circle((x,y), radius=138/60., facecolor='k',
                                            edgecolor='none', zorder=3, alpha=0.05))

        # Big IFU
        if show_all:
            ax.add_patch(patches.RegularPolygon((9.6,0.), 6, radius=0.3, facecolor='C0',
                         edgecolor='C0', zorder=4, label='Monolithic IFU', alpha=0.3))

        # FOBOS field-of-view
        ax.add_patch(patches.Circle((0.,0.), radius=10., facecolor='none', edgecolor='C3',
                                    zorder=4, label='FOBOS FOV'))
        if legend:
            ax.legend(ncol=2)
        if show:
            pyplot.show()
        return None if show else ax



