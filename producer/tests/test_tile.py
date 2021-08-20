
from IPython import embed

import numpy

from sklearn.neighbors import KDTree

from producer import tile

def test_uniform():
    rng = numpy.random.default_rng(99)
    ntarg = 80000
    x, y = rng.uniform(low=1, high=3, size=(2,ntarg))
    grid_coo = tile.uniform_tiles(x, y)

    assert grid_coo.shape[0] == 73, 'Different number of pointings'

    kdtree = KDTree(numpy.column_stack((x, y)))
    groups = kdtree.query_radius(grid_coo, 10/60)
    assert numpy.all(numpy.isin(numpy.arange(ntarg), numpy.unique(numpy.concatenate(groups)))), \
                'Did not cover all targets'

    # from matplotlib import pyplot, patches
    # w,h = pyplot.figaspect(1)
    # fig = pyplot.figure(figsize=(2.*w,2.*h))
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # ax.set_xlim([0.5, 3.5])
    # ax.set_ylim([0.5, 3.5])
    # ax.scatter(x, y, marker='.', s=30, lw=0, color='k')
    # for i, c in enumerate(grid_coo):
    #     ax.add_patch(patches.Circle((c[0],c[1]), radius=10./60, facecolor='C0', edgecolor='none',
    #                                 zorder=0, alpha=0.2))
    # pyplot.show()



