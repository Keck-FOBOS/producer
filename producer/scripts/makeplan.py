"""
Main execution script for FOBOS Observation Planning.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from . import scriptbase

class MakePlan(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='FOBOS Observation Planning', width=width)
        parser.add_argument('target_file', type=str,
                            help='Name of a file with targets to be observed')

        # TODO:
        #   - Allow the user to provide their own tile pattern.
        #   - Allow user to provide ID column
        #   - Allow for fits table input?
        parser.add_argument('-m', '--mode', nargs='*', default=None, type=int,
                            help='The spectrograph mode.  If None, mode is determined/optimized '
                                 'based on the aperture types identified in the target file.  '
                                 'If a single integer, all three spectrographs are put in the '
                                 'same mode.  Otherwise, must be 3 integers.  Use mode=1 for '
                                 'single-fiber apertures, mode=2 for IFUs.')

        parser.add_argument('-c', '--columns', nargs='*', default=[1,2], type=int,
                            help='Provide the columns with the RA, DEC, and aperture type '
                                 '(optional) to use for each target.  If the aperture type is '
                                 'not provided in the file, it is assumed to be the same for all '
                                 'targets and set by the --ap_type option.')

        parser.add_argument('-a', '--ap_type', default=0, type=int,
                            help='If the aperture type is not available in the target file, use '
                                 'this type for all targets.  Aperture types must be 0 for '
                                 'single fibers and 1 for IFUs.')

        parser.add_argument('-r', '--root', default=None, type=str,
                            help='Root name for output field configuration files.  If None, '
                                 'based on the name of the provided target file.')

        parser.add_argument('-o', '--offset', default=False, action='store_true',
                            help='Offset the default tiling by half the FOV')

        parser.add_argument('--min_n', default=0, type=int,
                            help='Minimum number of objects in a tile to consider it viable.')

        parser.add_argument('--tile_only', default=False, action='store_true',
                            help='Only show the tiling of the input coordinates.')

        parser.add_argument('--max_nobs', default=None, type=int,
                            help='Maximum number of times to revisit a FOBOS pointing to '
                                 'acquire new targets.  If None, pointings will be revisited '
                                 'until there are no more objects to observe.')

        return parser

    @staticmethod
    def main(args):

        from pathlib import Path

        import numpy

        from astropy import units
        from astropy import coordinates

        from ..targets import parse_targets
        from .. import plan
        from ..tile import uniform_tiles, show_tiles
        from ..astrometry import focal_plane_offsets

        if len(args.columns) == 2:
            ra_c, dec_c = args.columns
            ap_c = None
        elif len(args.columns) == 3:
            ra_c, dec_c, ap_c = args.columns
        else:
            raise ValueError('Must provide either RA and DEC columns, or RA, DEC, and AP columns.')

        target_file = Path(args.target_file).resolve()
        if not target_file.is_file():
            raise FileNotFoundError(f'{str(target_file)} does not exist!')

        root = target_file.parent / target_file.stem \
                    if args.root is None else Path(args.root).resolve()
        if not root.parent.is_dir():
            root.parent.mkdir()
            
        # TODO: Allow for different exposure time per object...
        print('-'*70)
        ra, dec, ap = parse_targets(str(target_file), ra_c=ra_c, dec_c=dec_c, ap_c=ap_c,
                                    default_ap=args.ap_type)
        ntarg = ra.size
        print(f'Read {ntarg} targets from {target_file.name}.')

        # TODO: Also use this to identify which targets are in each tile?
        tiles = uniform_tiles(ra, dec, half_offset=args.offset, min_n=args.min_n)
        print(f'Constructed {tiles.shape[0]} tiles.')
        if args.tile_only:
            show_tiles(tiles, ra=ra, dec=dec)
            return

        ndig = int(numpy.ceil(numpy.log10(tiles.shape[0]+1)))

        # TODO:
        #   - Remove objects from the available list as they're observed in each
        #     tile
        #   - Start with highest density tiles?
        #   - Only perform one observation for each tile, then re-tile based on
        #   the remaining objects.
        #   - Allow for different exposure time per object...
        for i, tile in enumerate(tiles):
            print('-'*70)
            print(f'Allocating FOBOS apertures in tile {i+1}.')
            # Offset coordinates from tile center in arcmin
            x, y = map(lambda x : x*60, focal_plane_offsets(ra, dec, tuple(tile)))
            n_in_fov, obs_obj, obs_nap, obs_ap, obs_mode \
                    = plan.configure_observations(x, y, ap, max_nobs=args.max_nobs)
            plan.report_configurations(n_in_fov, obs_obj, obs_nap, obs_ap, obs_mode)
            plan.write_configurations(f'{str(root)}_{i+1:0{ndig}}', ra, dec, tile, obs_obj, obs_ap,
                                      obs_mode, tight=True, target_file=str(target_file),
                                      ra_c=ra_c, dec_c=dec_c)
            print('-'*70)



