from pathlib import Path
import numpy
from producer import targets
from producer.tests.util import test_data_file

def main():
    p = test_data_file(filename='random_targets.db')
    # Radius is 1/3 degree (two adjacent FOBOS pointings)
    # Density is 5 per square arcmin
    x, y = targets.random_targets(1/3., density=5*60**2)
    numpy.savetxt(str(p), numpy.column_stack((x,y)), fmt=' %8.4f %8.4f')

if __name__ == '__main__':
    main()

