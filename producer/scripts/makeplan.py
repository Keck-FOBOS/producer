"""
Main execution script for FOBOS Observation Planning.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from . import scriptbase


class MakePlan(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='FOBOS Observation Planning', width=width)
        parser.add_argument('target_file', type=str,
                            help='Name of a file with targets to be observed')

        parser.add_argument('-m', '--mode', nargs='*', default=None,
                            help='The spectrograph mode.  If None, mode is determined/optimized '
                                 'based on the aperture types identified in the target file.  '
                                 'If a single integer, all three spectrographs are put in the '
                                 'same mode.  Otherwise, must be 3 integers.  Use mode=1 for '
                                 'single-fiber apertures, mode=2 for IFUs.')

        parser.add_argument('-c', '--columns', nargs='*', default=[1,2],
                            help='Provide the columns with the RA, DEC, and aperture type '
                                 '(optional) to use for each target.  If the aperture type is '
                                 'not provided in the file, it is assumed to be the same for all '
                                 'targets and set by the --ap_type option.')

        parser.add_argument('-a', '--ap_type', default=1, type=int,
                            help='If the aperture type is not available in the target file, use '
                                 'this type for all targets.  Aperture types must be 1 for '
                                 'single fibers and 2 for IFUs.')

        parser.add_argument('-r', '--root', default=None, type=str,
                            help='Root name for output field configuration files.  If None, '
                                 'based on the name of the provided target file.')

        parser.add_argument('-s', '--show', default=False, action='store_true',
                            help='Show field configurations and observing sequence summary plots.')

        return parser

    @staticmethod
    def main(args):

        from .targets import parse_targets

        if len(args.columns) == 2:
            ra_c, dec_c = args.columns
        elif len(args.columns) == 3:
            ra_c, dec_c, ap_c = args.columns
        else:
            raise ValueError('Must provide either RA and DEC columns, or RA, DEC, and AP columns.')

        ra, dec, ap = parse_targets(args.target_file, ra_c=ra_c, dec_c=dec_c, ap_c=ap_c,
                                    default_ap=args.ap_type)


