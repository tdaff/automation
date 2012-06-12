"""
Job handler

Machine specific job submission and tracking routines. Provides a unified
interface for all job handling. Implements the JobHandler class which will be
initialized to the machine the calculations are running on from the provided
options.

"""

import os
import re
import shlex
import subprocess
import sys
import time
from logging import debug, info, warn, error
from subprocess import Popen, PIPE, STDOUT


class JobHandler(object):
    """
    Abstraction of batch scheduler submission.

    """

    def __init__(self, options):
        """Initialize for machine specified by options."""
        self.queue = options.get('queue')
        if self.queue == 'wooki':
            self.submit = _wooki_submit
            self.postrun = _wooki_postrun
            self.jobcheck = _wooki_jobcheck
            self.env = _wooki_env
        elif self.queue == 'sharcnet':
            self.submit = _sharcnet_submit
            self.postrun = _sharcnet_postrun
            self.jobcheck = _sharcnet_jobcheck
            self.env = _sharcnet_env
        elif self.queue == 'local':
            self.submit = _local_run
            self.postrun = _pass
            self.jobcheck = _pass
            self.env = _pass
        else:
            print("ERROR unknown queue %s. Using null handler" % self.queue)
            self.submit = _pass
            self.postrun = _pass
            self.jobcheck = _pass
            self.env = _pass

    def _pbs_submit(self, job_type, nodes, **kwargs):
        """Submit a generic pbs job; return the jobid."""
        pass


def _sharcnet_submit(job_type, options, input_file=None, input_args=None):
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
    sqsub_args.extend(['-q', 'NRAP_20405'])
    # job_name
    sqsub_args.extend(['-j', 'faps-%s-%s' % (job_name, job_type)])
    # Is it a multiple CPU job?
    # Memory is mandatory; set depending on job type...
    if nodes > 1:
        # request nodes
        sqsub_args.extend(['-n', '%i' % nodes])
        if job_type in openmp_codes:
            # Some jobs are only openmp
            _check_program(exe)
            sqsub_args.extend(['-f', 'threaded'])
            sqsub_args.extend(['--mpp=%fg' %
                               options.getfloat('threaded_memory')])
        else:
            # Ensure mpi is enebaled
            _check_program(exe, mpi=True)
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
        _check_program(exe)
        sqsub_args.extend(['--mpp=%fg' % options.getfloat('threaded_memory')])
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
    submit = Popen(sqsub_args, stdout=PIPE)
    for line in submit.stdout.readlines():
        if 'submitted as' in line:
            jobid = int(line.split()[-1])
            break
    else:
        error("Job submission failed, no jobid received. Faps will crash now")

    return jobid


def _sharcnet_postrun(waitid):
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
        '-q', 'NRAP_20405',
        '-r', '20m',
        '-o', 'faps-post-%s.out' % '-'.join(sorted(waitid)),
        '--mpp=3g',
        '--waitfor=%s' % ','.join(waitid),
        ] + _argstrip(sys.argv)
    # We can just call this as we don't care about the jobid
    debug("Postrun command: %s" % " ".join(sqsub_args))
    subprocess.call(sqsub_args)


def _sharcnet_jobcheck(jobid):
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


def _sharcnet_env(code, *args, **kwargs):
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


def _wooki_submit(job_type, options, *args, **kwargs):
    """
    Interface to the submission scripts on wooki. Let them deal with the
    node types, because it is too fiddly to care about.
    """
    submit_scripts = {
        'vasp': 'vaspsubmit-faps',
        'repeat': 'repeatsubmit-faps',
        'siesta': 'siestasubmit-faps',
        'fastmc': 'fastmcsubmit-faps',
        'gulp': 'gulpsubmit-faps',
        'egulp': 'egulpsubmit-faps'
    }
    job_name = options.get('job_name')
    try:
        nodes = options.getint('%s_ncpu' % job_type)
    except AttributeError:
        nodes = 1

    submit_args = [submit_scripts[job_type], job_name, "%i" % nodes]

    if nodes != 1:
        submit_args.append("%i" % nodes)

    debug("Submission command: %s" % " ".join(submit_args))
    submit = Popen(submit_args, stdout=subprocess.PIPE)
    for line in submit.stdout.readlines():
        if "wooki" in line:
            jobid = line.split(".")[0]
            break
    else:
        print("Job submission failed?")

    return jobid


def _wooki_postrun(waitid):
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
    pbs_script = ['#PBS -N faps-post-%s\n' % '-'.join(sorted(waitid)),
                  '#PBS -m n\n',
                  '#PBS -o faps-post-%s.out\n' % '-'.join(sorted(waitid)),
                  '#PBS -j oe\n',
                  '#PBS -W depend=afterok:%s\n' % ':'.join(waitid),
                  'cd $PBS_O_WORKDIR\n',
                  'python ',
                  ' '.join(_argstrip(sys.argv))]

    pbs_script = ''.join(pbs_script)
    submit = Popen("/usr/local/bin/qsub", shell=False, stdin=PIPE)
    submit.communicate(input=pbs_script)


def _wooki_jobcheck(jobid):
    """Return true if job is still running or queued, or check fails."""
    # can deal with jobid as an int or a string
    jobid = ("%s" % jobid).strip()
    running_status = ['Q', 'R', 'Z']
    qstat = Popen(['/usr/local/bin/qstat', jobid], stdout=PIPE, stderr=STDOUT)
    for line in qstat.stdout.readlines():
        if "Unknown Job Id" in line:
            # Job finished and removed
            return False
        elif jobid in line:
            # use of 'in' should be fine as only this job will be shown
            status = line[68:69]
            if status in running_status:
                return True
            else:
                # Not running
                return False
    else:
        print("Failed to get job information.")  # qstat parsing failed?
        # Act as if the job is still running, in case it hasn't finished
        return True


def _wooki_env(code, *args, **kwargs):
    """Hacks to get things working with wooki submission scripts"""
    pass


def _local_run(job_type, options, input_file=None):
    """Run the exe in a subprocess."""
    # Threaded codes have different behaviour
#    openmp_codes = options.gettuple('threaded_codes')

    # Bind some things locally, so we know what's going on
    job_name = options.get('job_name')
    exe = options.get('%s_exe' % job_type)
    try:
        nodes = options.getint('%s_ncpu' % job_type)
    except AttributeError:
        nodes = 1

    # Some codes need the input file name
    if input_file is not None:
        sqsub_args.extend(['-i', '%s' % input_file])
    # Output
    out_file = open('faps-%s.out' % job_name, 'wb')

    run_args = shlex.split(exe)

    submit = Popen(sqsub_args, stdout=out_file)
    submit.wait()
    finished = submit.returncode

    return None


def _pass(*args, **kwargs):
    """Sometimes we want to do nothing."""
    pass


def _argstrip(arglist):
    """Some options might be best removed before resubmission."""
    to_remove = ['-i', '--import', 'dft', 'gcmc', 'charges']
    newargs = list(arglist)
    for item in to_remove:
        while item in newargs:
            newargs.remove(item)
    return newargs


def _check_program(program, mpi=False):
    """Test to see if the exe is in the path and might work."""
    exe = which(program)
    if exe is None:
        error("Could not find %s in path, job will probably fail." % program)
        return False
    elif mpi:
        try:
            binary = open(exe, 'rb')
            if re.search('mpi_init', binary.read(), re.IGNORECASE):
                return True
            else:
                warn("%s doesn't appear to be an mpi executable, job may fail" % program)
                return False
        except IOError:
            return False
    else:
        # TODO(tdaff): No easy way to check for OMP in general?
        return True


def which(program):
    """Return the equivalent of the 'which' command."""
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
