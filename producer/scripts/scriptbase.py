"""
Implements base classes for use with scripts.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import argparse
import textwrap
from functools import reduce

class SmartFormatter(argparse.HelpFormatter):
    r"""
    Enable a combination of both fixed-format and wrappable lines to be
    formatted for the help statements for command-line arguments used with
    `argparse.ArgumentParser`_.

    Borrows from
    https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text

    Help strings that use this formatter *must* begin with "R|".  If not, the
    help string is parsed by the base class.

    When parsed by this formatter, the leading "R|" characters are stripped and
    the lines to be printed are parsed using `str.splitlines`_.  Each resulting
    line is wrapped using `textwrap.wrap`_, unless it begins with the characters
    "F|", which forces the line to remain unaltered (except for stripping the
    leading characters).

    For example, if you add an argument like this:

    .. code-block:: python

        parser.add_argument('-t', '--tell_file', type=str,
                            help='R|Configuration file to change default telluric parameters.  '
                                 'Note that the parameters in this file will be overwritten if '
                                 'you set argument in your terminal.  The --tell_file option '
                                 'requires a .tell file with the following format:\n'
                                 '\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = qso\n'
                                 'F|         redshift = 7.6\n'
                                 'F|         bal_wv_min_max = 10825,12060\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = star\n'
                                 'F|         star_type = A0\n'
                                 'F|         star_mag = 8.\n'
                                 'OR\n'
                                 'F|    [tellfit]\n'
                                 'F|         objmodel = poly\n'
                                 'F|         polyorder = 3\n'
                                 'F|         fit_wv_min_max = 9000.,9500.\n'
                                 '\n')

    The result will be (depending on the width of your console):

    .. code-block:: console

        -t TELL_FILE, --tell_file TELL_FILE
                          Configuration file to change default telluric
                          parameters.  Note that the parameters in this file
                          will be overwritten if you set argument in your
                          terminal.  The --tell_file option requires a .tell
                          file with the following format:

                              [tellfit]
                                   objmodel = qso
                                   redshift = 7.6
                                   bal_wv_min_max = 10825,12060
                          OR
                              [tellfit]
                                   objmodel = star
                                   star_type = A0
                                   star_mag = 8.
                          OR
                              [tellfit]
                                   objmodel = poly
                                   polyorder = 3
                                   fit_wv_min_max = 9000.,9500.
    """
    def _split_lines(self, text, width):
        """
        Split the provided text into width constrained lines.

        See the class description for formatting instructions.
        """
        if text.startswith('R|'):
            lines = text[2:].splitlines()
            for i in range(len(lines)):
                if lines[i].startswith('F|'):
                    lines[i] = [lines[i][2:]]
                elif len(lines[i]) == 0:
                    lines[i] = [' ']
                else:
                    lines[i] = textwrap.wrap(lines[i], width)
            return reduce(list.__add__, lines)
        return super()._split_lines(text, width)


class ScriptBase:
    """
    Provides a base class for all scripts.
    """
    @classmethod
    def entry_point(cls):
        """
        Defines the main script entry point.
        """
        cls.main(cls.parse_args())

    # TODO: Combining classmethod and property works in python 3.9 and later
    # only: https://docs.python.org/3.9/library/functions.html#classmethod
    # Order matters.  In python 3.9, it would be:
    #
    # @classmethod
    # @property
    #
    # Because we're not requiring python 3.9 yet, we have to leave this as a
    # classmethod only:
    @classmethod
    def name(cls):
        """
        Provide the name of the script.  By default, this is the name of the
        module with "fobos" prepended.
        """
        return f"fobos_{cls.__module__.split('.')[-1]}"

    @classmethod
    def parse_args(cls, options=None):
        """
        Parse the command-line arguments.
        """
        parser = cls.get_parser()
        return parser.parse_args() if options is None else parser.parse_args(options)

    # Base classes should override this
    @staticmethod
    def main(args):
        """
        Execute the script.
        """
        pass

    # Base classes should override this.  Ideally they should use this
    # base-class method to intantiate the ArgumentParser object and then fill in
    # the relevant parser arguments
    @classmethod
    def get_parser(cls, description=None, width=None,
                   formatter=argparse.ArgumentDefaultsHelpFormatter):
        """
        Construct the command-line argument parser.

        Args:
            description (:obj:`str`, optional):
                A short description of the purpose of the script.
            width (:obj:`int`, optional):
                Restrict the width of the formatted help output to be no longer
                than this number of characters, if possible given the help
                formatter.  If None, the width is the same as the terminal
                width.
            formatter (`argparse.HelpFormatter`_):
                Class used to format the help output.

        Returns:
            `argparse.ArgumentParser`_: Command-line interpreter.
        """
        return argparse.ArgumentParser(description=description,
                                       formatter_class=lambda prog: formatter(prog, width=width))



