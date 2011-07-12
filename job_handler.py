"""
Job handlers

Machine specific job submission and tracking routines. Implements the
JobHandler class which should be initialized to the machine the calculations
are running on.

"""

from subprocess import Popen, PIPE, STDOUT

class JobHandler(object):
    """
    Abstraction of the different submission routines management and scripting.

    """
    #TODO(tdaff): remember to re-run the script once job is finished!

    def __init__(self, options):
        """Initialize for machine specified by options."""
        self.queue = options.get('queue')
        if self.queue == 'wooki':
            self.submit = self._wooki_submit
            self.jobcheck = self.wooki_jobcheck
        else:
            self.submit = self._pbs_submit
            self.submit = self._pbs_jobcheck

    def _wooki_submit(self, job_type, nodes, **kwargs):
        """Submit a job to wooki; return the jobid."""
        pass

    def _pbs_submit(self, job_type, nodes, **kwargs):
        """Submit a generic pbs job; return the jobid."""
        pass

    def _wooki_jobcheck(self, jobid):
        """Get job status. Return status or False for non existent job."""
        jobid = "%s" % jobid
        qstat = Popen(['qstat', '%s' % jobid], stdout=PIPE, stderr=STDOUT)
        for line in qstat.stdout.readlines():
            if "Unknown Job Id" in line:
                return False
            elif line.startswith(jobid):
                status = line[68:69]
                return status
        else:
            print("Failed to get job information.")  # qstat parsing failed?
