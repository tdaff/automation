"""
faps.backend.file interfaces the simulation with a simple file
reader and writer.

"""

import os.path
import json

from faps.backend import Backend
from faps.core.simulation import Simulation
from faps.settings import info, debug, error
from faps.util.system import terminate

class FileBackend(Backend):
    """
    Convert simulation to and from a JSON file on the disk.
    
    """

    @staticmethod
    def load(jobname):
        json_filename = "{}.json".format(simulation.jobname)
        try:
            with open(json_filename) as json_file:
                simulation_dict = json.load(json_file)
                return Simulation(jobname, drepr=simulation_dict)
        except ValueError:
            # Invalid json
            error("Invalid json encountered in {}.".format(json_filename))
            error("Cannot continue")
            terminate(44)
        except IOError:
            if os.path.exists(json_filename):
                error("Unable to load file {}".format(json_filename))
                error("Cannot continue")
                terminate(45)
            else:
                info("Starting new simulation")
                return Simulation(jobname)

    @staticmethod
    def save(simulation):
        json_filename = "{}.json".format(simulation.jobname)
        with open(json_filename) as json_file:
