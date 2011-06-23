"""
Job handlers

Machine specific job submission and tracking routines. Implements the
JobHandler class which should be initialized to the machine the calculations
are running on.
"""

class JobHandler(object):
    """
    Abstraction of the different submission routines management and scripting.

    """
    #TODO(tdaff): remember to re-run the script once job is finished!

    def __init__(self, options):
        """Initialize for machine specified by options."""
        pass

    def submit(self, options):
        pass

    def _wooki(self):
        pass
    
    def _orca(self):
        pass
