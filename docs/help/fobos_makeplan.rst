.. code-block:: console

    $ fobos_makeplan -h
    usage: fobos_makeplan [-h] [-m [MODE ...]] [-c [COLUMNS ...]] [-a AP_TYPE]
                          [-r ROOT] [-o] [--min_n MIN_N] [--tile_only]
                          [--max_nobs MAX_NOBS]
                          target_file
    
    FOBOS Observation Planning
    
    positional arguments:
      target_file           Name of a file with targets to be observed
    
    optional arguments:
      -h, --help            show this help message and exit
      -m [MODE ...], --mode [MODE ...]
                            The spectrograph mode. If None, mode is
                            determined/optimized based on the aperture types
                            identified in the target file. If a single integer, all
                            three spectrographs are put in the same mode. Otherwise,
                            must be 3 integers. Use mode=1 for single-fiber
                            apertures, mode=2 for IFUs. (default: None)
      -c [COLUMNS ...], --columns [COLUMNS ...]
                            Provide the columns with the RA, DEC, and aperture type
                            (optional) to use for each target. If the aperture type
                            is not provided in the file, it is assumed to be the
                            same for all targets and set by the --ap_type option.
                            (default: [1, 2])
      -a AP_TYPE, --ap_type AP_TYPE
                            If the aperture type is not available in the target
                            file, use this type for all targets. Aperture types must
                            be 0 for single fibers and 1 for IFUs. (default: 0)
      -r ROOT, --root ROOT  Root name for output field configuration files. If None,
                            based on the name of the provided target file. (default:
                            None)
      -o, --offset          Offset the default tiling by half the FOV (default:
                            False)
      --min_n MIN_N         Minimum number of objects in a tile to consider it
                            viable. (default: 0)
      --tile_only           Only show the tiling of the input coordinates. (default:
                            False)
      --max_nobs MAX_NOBS   Maximum number of times to revisit a FOBOS pointing to
                            acquire new targets. If None, pointings will be
                            revisited until there are no more objects to observe.
                            (default: None)
    