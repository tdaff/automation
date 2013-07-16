"""
Serial Job handler

Runs external commands in-process without submitting to a queue system.

"""

import shlex
from logging import info
from subprocess import Popen, PIPE

from faps.managers import JobManager
from faps.settings import options

MAX_RETRY = 5

class SerialJobManager(JobManager):
    """
    Submit commands in-process and wait for them to finish.

    """

    def submit(job_type, input_file=None, input_args=None):
        """Run the exe in a subprocess. input_args must be a list if """

        # Bind some things locally, so we know what's going on
        job_name = options.get('job_name')
        exe = options.get('%s_exe' % job_type)
        serial_args = shlex.split(exe)

        if input_args is not None:
            serial_args.extend(input_args)

        # Some codes need the input file name
        if input_file is not None:
            input_file = open('%s' % input_file)

        # Output
        out_file = open('faps-%s.out' % job_name, 'wb')

        # run the job in process
        debug("Waiting for command: %s" % " ".join(serial_args))
        submit = Popen(serial_args, stdin=input_file, stdout=out_file)
        submit.wait()
        finished = submit.returncode
        info("%s job finished with return code %s" % (exe, finished))

        # always return True for a finished job
        return True

    def postrun(jobid = None, jobids = None):
        """Nothing to do here."""
        pass

    def jobcheck(jobid):
        """No job should be running."""
        return False

    def env(job_type):
        """Nothing to do here."""
        pass
