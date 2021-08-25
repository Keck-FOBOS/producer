"""
Testing utilities
"""

from producer import data_file

def test_data_file(filename=None):
    root = data_file() / 'tests'
    return root if filename is None else root / filename



