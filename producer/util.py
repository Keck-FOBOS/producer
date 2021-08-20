"""
Miscellaneous utilities.

.. include:: ../include/links.rst
"""

from itertools import chain, combinations

from IPython import embed 

def powerset(iterable, reverse=False):
    """"
    Construct an iterable that steps through all combinations of the
    provided iterable.

    This is pulled from the recipes provided by the itertools
    documentation.

    Examples:
        
        Get all unique combinations of the list [1,2,3]:
        >>> list(powerset([1,2,3]))
        [() (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)]

    Args:
        iterable (iterable):
            An iterable object
        reverse (:obj:`bool`, optional):
            Reverse the order (only roughly) of the iterable by placing
            the longer sequences first.
    
    Returns:
        `itertools.chain`: Iterable object that returns the sequence of
        combinations.
    """
    rng = range(len(iterable)+1)[::-1] if reverse else range(len(iterable)+1)
    return chain.from_iterable(combinations(iterable, r) for r in rng)



