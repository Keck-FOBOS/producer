"""
Main execution script for FOBOS Observation Planning.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from . import scriptbase

class ShowCfg(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='Show FOBOS observation configuration',
                                    width=width)
        parser.add_argument('config_file', type=str, help='Name of the configuration file')

        parser.add_argument('-f', '--full', default=False, action='store_true',
                            help='If available, re-read the file with the original targets used '
                                 'to generate the configuration and show them in the plot.')
        return parser

    @staticmethod
    def main(args):

        import numpy

        from .. import plan
        from ..astrometry import focal_plane_offsets
        from ..allocate import assign_apertures_plot

        # NOTE: parse_configuration checks that the file exists
        objid, ra, dec, center, assigned_obj, ap, assigned_ap \
                = plan.parse_configuration(args.config_file, parse_source=args.full)

        ignore_ap = numpy.where(numpy.logical_not(ap.select('science')))

        x, y = map(lambda x : x*60, focal_plane_offsets(ra, dec, center))
        assign_apertures_plot(x, y, ap.coo[:,0], ap.coo[:,1], assigned_obj, assigned_ap,
                              ignore_ap=ignore_ap)


