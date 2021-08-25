
from IPython import embed

from producer import data_file
from producer.design.fobosfocalplane import fobos_module_layout
from producer.deploy import FOBOSApertures


if __name__ == '__main__':
    fobos_module_layout(str(FOBOSApertures.module_src), str(FOBOSApertures.starbug_src),
                        version=FOBOSApertures.version)


