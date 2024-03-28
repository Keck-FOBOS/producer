.. code-block:: console

    $ fobos_allocsim -h
    usage: fobos_allocsim [-h] [-d DENSITY DENSITY DENSITY] [-s SIMS] [-m MODE] [-e]
    
    FOBOS Allocation Simulation
    
    options:
      -h, --help            show this help message and exit
      -d DENSITY DENSITY DENSITY, --density DENSITY DENSITY DENSITY
                            Target density sampling: minimum, maximum, number of
                            samples. Density is sampled geometrically. (default:
                            [2.5, 40, 5])
      -s SIMS, --sims SIMS  Number of simulations to use for mean and standard
                            deviation of trends. (default: 1)
      -m MODE, --mode MODE  The spectrograph mode. All spectrographs are put in the
                            same mode. Use 1 for MOS mode, 2 for multi-IFU mode.
                            (default: 1)
      -e, --embed           Embed using IPython before completing the script.
                            (default: False)
    