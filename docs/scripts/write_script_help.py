"""
Dynamically build the rst documentation with the script help text.
"""

import time
from importlib import resources

from producer.scripts import script_classes


#-----------------------------------------------------------------------------

def write_help(script_cls, opath, width=80):
    exe = script_cls.name()
    ofile = opath / f'{exe}.rst'
    lines = ['.. code-block:: console', '']
    lines += ['    $ {0} -h'.format(exe)]
    parser = script_cls.get_parser(width=80)
    parser.prog = exe
    lines += ['    ' + l for l in parser.format_help().split('\n')]
    print('Writing: {0}'.format(ofile))
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    root = resources.files('producer').parent
    path = root / 'docs' / 'help'
    if not path.is_dir():
        path.mkdir(parents=True)

    # Get the list of script names and script classes
    scr_clss = script_classes()
    for name, script_cls in scr_clss.items():
        write_help(script_cls, path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


