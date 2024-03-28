.. code-block:: console

    $ fobos_layout -h
    usage: fobos_layout [-h] [--mode [MODE ...]] [--design DESIGN]
    
    Show FOBOS aperture deployment for a given mode
    
    options:
      -h, --help         show this help message and exit
      --mode [MODE ...]  FOBOS spectrograph mode(s). Can be a single integer,
                         defining a single mode for all spectrographs, or a set of
                         three integers that set the mode for each spectrograph
                         individually. Spectrograph modes are: (1) single-fiber
                         apertures, (2) 37-fiber IFUs, (3) monolithic IFU. (default:
                         1)
      --design DESIGN    Focal-plane layout design. Currently testing two designs.
                         One where the spectrograph-to-module mapping distributes
                         all three spectrographs over the full FOBOS field-of-view
                         (design 1) and one where the modules for a given
                         spectrograph are confined to a coherent region in the focal
                         plane. (default: 1)
    