"""
Job handlers

Machine specific job submission and tracking routines. Implements the
JobHandler class which should be initialized to the machine the calculations
are running on.

"""

from subprocess import Popen, PIPE

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

    def jobcheck(self, jobid):
        """Get job status."""
        jobid = "%s" % jobid
        qstat = Popen(['qstat', '%s' % jobid], stdout=PIPE)
        for line in qstat.stdout.readlines():
            if "Unknown Job Id" in line:
                return False
            elif line.startswith(jobid):
                status = line[68:69]
                return status
        else:
            print("Failed to get job information.")
