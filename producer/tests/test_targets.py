
from IPython import embed

import numpy

from producer import targets
from producer.tests.util import test_data_file


def test_random():
    n = 100
    r = 10.
    x, y = targets.random_targets(r, n=n)
    assert x.size == n, 'Incorrect number of targets'
    assert numpy.all(x**2 + y**2 < r**2), 'Not limited to circle'


def test_parse():
    n = 100
    r = 10.
    rng = numpy.random.default_rng(99)
    x, y = targets.random_targets(r, n=n, rng=rng)
    opath = test_data_file('test_parse.tmp')
    opath.parent.mkdir(exist_ok=True)
    numpy.savetxt(str(opath), numpy.column_stack((x,y)), fmt=' %8.4f %8.4f')
    _x, _y, _a = targets.parse_targets(str(opath), default_ap=1)

    assert numpy.allclose(numpy.around(x, decimals=4), _x), 'Bad read'
    assert numpy.all(_a == 1), 'Aperture bad'

    a = rng.integers(low=1, high=3, size=n)
    numpy.savetxt(str(opath), numpy.column_stack((x,y,a)), fmt=' %8.4f %8.4f %2d')
    _x, _y, _a = targets.parse_targets(str(opath), ap_c=3, default_ap=1)
    assert numpy.array_equal(a, _a), 'Aperture bad'

    # Clean up
    opath.unlink()


