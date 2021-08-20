"""
Testing utilities
"""

from pathlib import Path
from pkg_resources import resource_filename

def data_file(filename=None):
    root = Path(resource_filename('producer', 'data')).resolve()
    return root if filename is None else root / filename

def test_data_file(filename=None):
    root = data_file() / 'tests'
    return root if filename is None else root / filename



