"""
Main execution script for FOBOS Observation Planning.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

from . import scriptbase

class Layout(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show FOBOS aperture deployment for a given mode',
                                    width=width)
        parser.add_argument('--mode', nargs='*', default=1, type=int,
                            help='FOBOS spectrograph mode(s).  Can be a single integer, defining '
                                 'a single mode for all spectrographs, or a set of three integers '
                                 'that set the mode for each spectrograph individually.  '
                                 'Spectrograph modes are: (1) single-fiber apertures, (2) '
                                 '37-fiber IFUs, (3) monolithic IFU.')
        parser.add_argument('--design', type=int, default=1,
                            help='Focal-plane layout design.  Currently testing two designs.  '
                                 'One where the spectrograph-to-module mapping distributes all '
                                 'three spectrographs over the full FOBOS field-of-view (design '
                                 '1) and one where the modules for a given spectrograph are '
                                 'confined to a coherent region in the focal plane.')
        return parser

    @staticmethod
    def main(args):

        from ..deploy import FOBOSApertures

        ap = FOBOSApertures(mode=args.mode, baseline=True, config=args.design)
        ap.show()


