"""
Wooki Job handler

Run jobs through the queue for wooki.

"""

import os
import re
import subprocess
import sys
import time
from logging import debug, error
from subprocess import Popen, PIPE, STDOUT

from faps.managers import JobManager
from faps.settings import options
from faps.util.path import argstrip, check_program

MAX_RETRY = 5


class SharcnetJobManager(JobManager):
    """
    Submit and monitor jobs on sharcnet clusters.

    """

    def submit(job_type, input_file=None, input_args=None):
        """Simple interface to the 'sqsub' submission on sharcnet"""
        # Threaded codes have different behaviour
        openmp_codes = options.gettuple('threaded_codes')

        # Bind some things locally, so we know what's going on
        job_name = options.get('job_name')
        exe = options.get('%s_exe' % job_type)
        try:
            nodes = options.getint('%s_ncpu' % job_type)
        except AttributeError:
            nodes = 1

        sqsub_args = ['sqsub']
        # Always use the dedicated queue; faster
        dedicated_queue = options.get('dedicated_queue')
        if dedicated_queue:
            sqsub_args.extend(['-q', dedicated_queue])
        # job_name
        sqsub_args.extend(['-j', 'faps-%s-%s' % (job_name, job_type)])
        # Is it a multiple CPU job?
        # Memory is mandatory; set depending on job type...
        if nodes > 1:
            # request nodes
            sqsub_args.extend(['-n', '%i' % nodes])
            if job_type in openmp_codes:
                # Some jobs are only openmp
                check_program(exe)
                sqsub_args.extend(['-f', 'threaded'])
                sqsub_args.extend(['--mpp=%fg' %
                                   options.getfloat('threaded_memory')])
            else:
                # Ensure mpi is enebaled
                check_program(exe, mpi=True)
                sqsub_args.extend(['-f', 'mpi'])
                # Pack up to 12 procs on a single node
                if nodes < 12:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['--pack'])
                elif nodes % 2 == 0 and nodes/2 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N2'])
                elif nodes % 4 == 0 and nodes/4 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N4'])
                elif nodes % 6 == 0 and nodes/6 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N6'])
                elif nodes % 8 == 0 and nodes/8 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N8'])
                elif nodes % 12 == 0 and nodes/12 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N12'])
                elif nodes % 24 == 0 and nodes/24 < 8:
                    sqsub_args.extend(['--mpp=1.66g'])
                    sqsub_args.extend(['-N24'])
                else:
                    sqsub_args.extend(['--mpp=2.66g'])

        else:
            check_program(exe)
            sqsub_args.extend(['--mpp=%fg' % options.getfloat('serial_memory')])
        # run-time estimate mandatory; 12 hours is plenty?
        sqsub_args.extend(['-r', '48h'])
        # Some codes need the input file name
        if input_file is not None:
            sqsub_args.extend(['-i', '%s' % input_file])
        # Output
        sqsub_args.extend(['-o', 'faps-%s.out' % job_name])
        # Which command?
        sqsub_args.extend([exe])
        if input_args is not None:
            sqsub_args.extend(input_args)

        debug("Submission command: %s" % " ".join(sqsub_args))
        submitted = False
        submit_count = 0
        while not submitted and submit_count <= MAX_RETRY:
            submit = Popen(sqsub_args, stdout=PIPE)
            for line in submit.stdout.readlines():
                if 'submitted as' in line:
                    jobid = int(line.split()[-1])
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
        # Check that job appears in sqjobs before submitting next one
        for loop_num in range(10):
            qstat = Popen(['qstat', '-u', '$USER'], stdout=PIPE, shell=True)
            if waitid.issubset(re.split('[\s.]', qstat.stdout.read())):
                # All jobs there
                break
            else:
                # Wait longer each time, in case the system is very slow
                time.sleep(loop_num)
        # Sumbit here, even if jobs never found in queue
        sqsub_args = [
            'sqsub',
            '-r', '20m',
            '-o', 'faps-post-%s.out' % '-'.join(sorted(waitid)),
            '--mpp=3g',
            '--waitfor=%s' % ','.join(waitid),
            ]
        #TODO(tdaff): options...
        if dedicated_queue:
            sqsub_args.extend(['-q', dedicated_queue])
        # Add the submitted program cleaned for instruction commands
        sqsub_args.extend(argstrip(sys.argv))
        # We can just call this as we don't care about the jobid
        debug("Postrun command: %s" % " ".join(sqsub_args))
        subprocess.call(sqsub_args)


    def jobcheck(jobid):
        """Return true if job is still running or queued, or check fails."""
        # can deal with jobid as an int or a string
        jobid = ("%s" % jobid).strip()
        running_status = ['Q', 'R', 'Z']
        qstat = Popen(['sqjobs', jobid], stdout=PIPE, stderr=STDOUT)
        for line in qstat.stdout.readlines():
            if "ERROR" in line:
                # Job finished and removed
                return False
            elif jobid in line:
                # use of 'in' should be fine as only this job will be shown
                # can't use positional slicing as columns resize
                status = line.split()[3]
                if status in running_status:
                    return True
                else:
                    return False
        else:
            error("Failed to get job information. Is sqjobs working?")
            # Act as if the job is still running, in case it hasn't finished
            return True


    def env(code, *args, **kwargs):
        """Update the running environment for specific codes."""
        if code == 'siesta':
            # Siesta uses pathscale libraries and is a module
            # the pathscale currently appends the 32 bit lib_dir
            # so we have to try and chop that out
            newenv = subprocess.Popen(
                ['/usr/bin/modulecmd', 'python', 'load', 'siesta'], stdout=PIPE)
            #newenv = newenv.stdout.read()
            exec(newenv.stdout)
            ld_lib = os.environ['LD_LIBRARY_PATH']
            ld_lib = re.sub('/opt/sharcnet/pathscale/(.*)/lib/(.*)/32',
                            '/opt/sharcnet/pathscale/\\1/lib/\\2', ld_lib)
            os.environ['LD_LIBRARY_PATH'] = ld_lib
        else:
            pass
