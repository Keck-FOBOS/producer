
from pathlib import Path
from pkg_resources import resource_filename

from IPython import embed

import numpy

from producer.deploy import FOBOSApertures
from producer.util import string_table

def write_deployment_table(ofile):
    data_table = numpy.empty((3, 7, 4), dtype=object)
    data_table[0, 0,:] = ['Payload', 'MOS', 'IFU', 'MONO']
    data_table[1, 0,:] = ['Payload', 'MOS', 'IFU', 'MONO']
    data_table[2, 0,:] = ['Payload', 'MOS', 'IFU', 'MONO']
    data_table[0,1:,0] = ['Designated Sky', 'Single-Fiber', '37-fiber IFU', '7-fiber Flux Cal. IFU',
                         '547-fiber Mono IFU', 'Total Fibers']
    data_table[1,1:,0] = ['Designated Sky', 'Single-Fiber', '37-fiber IFU', '7-fiber Flux Cal. IFU',
                         '547-fiber Mono IFU', 'Total Fibers']
    data_table[2,1:,0] = ['Designated Sky', 'Single-Fiber', '37-fiber IFU', '7-fiber Flux Cal. IFU',
                         '547-fiber Mono IFU', 'Total Fibers']

    nap = numpy.zeros((3,5), dtype=int)
    nfib = numpy.array([1, 1, 37, 7, 547])

    ap = [FOBOSApertures(mode=1), FOBOSApertures(mode=2), FOBOSApertures(mode=3)]

    # Loop through the spectrographs
    for i in range(3):
        # Loop through the modes
        for j in range(3):
            spc_indx = ap[j].spc == i+1
            # Number of sky apertures
            nap[j,0] = numpy.sum(ap[j].select('sky') & spc_indx)
            # Number of single-fiber apertures
            nap[j,1] = numpy.sum(ap[j].select('science') & (ap[j].payload == 0) & spc_indx)
            # Number of 37-fiber IFUs
            nap[j,2] = numpy.sum(ap[j].select('science') & (ap[j].payload == 1) & spc_indx)
            # Number of 7-fiber calibration bundles
            nap[j,3] = numpy.sum(ap[j].active & (ap[j].payload == 4) & spc_indx)
            # Number of 547-fiber monolithic bundles
            nap[j,4] = int(j == 2)

        data_table[i,1:6,1:] = nap.T.astype(str)
        data_table[i,6,1:] = numpy.sum(nap * nfib[None,:], axis=1).astype(str)

    lines = '\n**Spectrograph 1**\n\n' + string_table(data_table[0], delimeter='rst') \
            + '\n\n**Spectrograph 2**\n\n' + string_table(data_table[1], delimeter='rst') \
            + '\n\n**Spectrograph 3**\n\n' + string_table(data_table[2], delimeter='rst') + '\n\n'

    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


def main():

    path = Path(resource_filename('producer', '')).parent / 'docs' / 'include'
    if not path.is_dir():
        path.mkdir()

    ofile = path / 'fobos_deployments.rst'
    write_deployment_table(str(ofile))


if __name__ == '__main__':
    main()
