"""
Deduce the dependency graph for a simulation.

"""

# Ensure that all the componenets have been included in the
# set of possible providers and consumers
from faps.settings import config
from faps.calculators import Calculator
from faps.io import StructureReader, StructureWriter
from faps.managers import JobManager

def resolve(simulation):
    """
    Take a simulation and work out the set of simulations that need
    to be carried out to resolve all the dependencies.

    """
    calcualtors = Calculator()
    strcuture_readers = StructureReader()
    structure_writers = StructureWriter()
    for target in config.targets:
        # get dependencies
        # see what needs to be done
        # find a minimal set
        # return it as a....
