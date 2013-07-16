"""
Wooki Job handler

Run jobs through the queue for wooki.

"""

import subprocess
import sys
import time
from logging import debug, warn, error
from subprocess import Popen, PIPE, STDOUT

from faps.managers import JobManager
from faps.settings import options
from faps.util.path import argstrip


MAX_RETRY = 5


class WookiJobManager(JobManager):
    """
    Sumbit and monitor jobs on the wooki cluster running customised
    scripts on SGE.

    """

    def submit(job_type, *args, **kwargs):
        """
        Interface to the submission scripts on wooki. Let them deal with the
        node types, because it is too fiddly to care about.
        """
        submit_scripts = {
            'vasp': 'vasp-submit',
            'repeat': 'repeat-submit-faps',
            'siesta': 'siesta-submit',
            'fastmc': 'fastmc-submit',
            'gulp': 'gulp-submit-faps',
            'egulp': 'egulp-submit'
        }
        job_name = options.get('job_name')
        try:
            nodes = options.getint('%s_ncpu' % job_type)
        except AttributeError:
            nodes = 1

        submit_args = [submit_scripts[job_type], job_name, "%i" % nodes]

        debug("Submission command: %s" % " ".join(submit_args))
        submitted = False
        submit_count = 0
        while not submitted and submit_count <= MAX_RETRY:
            submit = Popen(submit_args, stdout=subprocess.PIPE)
            for line in submit.stdout.readlines():
                # Your job 123 ("vasp.testzif") has been submitted
                if "has been submitted" in line:
                    jobid = line.split()[2]
                    submitted = True
                    break
            else:
                submit_count += 1
                error("Job submission attempt %i failed." % submit_count)
                time.sleep(submit_count)

        return jobid


    def postrun(waitid):
        """
        Resubmit this script for the postrun on job completion. Will accept
        a single jobid or a list, as integers or strings.
        """
        # Magic makes everything into a set of strings
        if hasattr(waitid, '__iter__'):
            waitid = frozenset([("%s" % wid).strip() for wid in waitid])
        else:
            waitid = frozenset([("%s" % waitid).strip()])
        # No jobcheck here as we assume wooki works
        sge_script = ['#!/bin/bash\n',
                      '#$ -cwd\n',
                      '#$ -V\n',
                      '#$ -j y\n',
                      '#$ -N faps-post-%s\n' % '-'.join(sorted(waitid)),
                      '#$ -o faps-post-%s.out\n' % '-'.join(sorted(waitid)),
                      '#$ -hold_jid %s\n' % ','.join(waitid),
                      'python ',
                      ' '.join(argstrip(sys.argv))]

        sge_script = ''.join(sge_script)
        submit = Popen("qsub", shell=False, stdin=PIPE)
        submit.communicate(input=sge_script)


    def jobcheck(jobid):
        """Return true if job is still running or queued, or check fails."""
        # can deal with jobid as an int or a string
        jobid = ("%s" % jobid).strip()
        qstat = Popen(['qstat', '-j', jobid], stdout=PIPE, stderr=STDOUT)
        for line in qstat.stdout.readlines():
            if "Following jobs do not exist" in line:
                # Job finished and removed
                return False
            #TODO(tdaff): check actual job status
        else:
            #TODO(tdaff): check if scheduler failed
            return True


    def env(code, *args, **kwargs):
        """Hacks to get things working with wooki submission scripts"""
        pass

