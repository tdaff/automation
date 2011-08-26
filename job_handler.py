"""
Job handler

Machine specific job submission and tracking routines. Implements the
JobHandler class which will be initialized to the machine the calculations
are running on from the provided options.

"""

import os
import getpass
import sys
from subprocess import Popen, PIPE, STDOUT


class JobHandler(object):
    """
    Abstraction of batch scheduler submission.

    """
    #TODO(tdaff): remember to re-run the script once job is finished!

    def __init__(self, options):
        """Initialize for machine specified by options."""
        self.queue = options.get('queue')
        if self.queue == 'wooki':
            self.submit = self._wooki_submit
            self.jobcheck = self.wooki_jobcheck
        elif self.queue == 'sharcnet':
            self.submit = _sharcnet_submit
            self.postrun = _sharcnet_postrun
            self.jobcheck = _sharcnet_jobcheck
        else:
            self.submit = self._pbs_submit
            self.jobcheck = self._pbs_jobcheck

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


def _sharcnet_submit(job_type, options):
    """Simple interface to the 'sqsub' submission on sharcnet"""
    # TODO(tdaff): self-resubmission?
    # sqsub -q DR_20293 -f mpi -n 48 -o std.out -j hmof-589 -r 6h --mpp=4g ~/bin/vasp-5.2.11-sequential
    # Bind some things locally, so we know what's going on
    job_name = options.get('job_name')
    exe = options.get('%s_exe' % job_type)
    try:
        nodes = options.getint('%s_ncpu' % job_type)
    except AttributeError:
        nodes = 1

    sqsub_args = ['sqsub']
    # Dedicated queue
    sqsub_args.extend(['-q', 'DR_20293'])
    # job_name
    sqsub_args.extend(['-j', 'faps-%s' % job_name])
    # Is it a multiple CPU job?
    if nodes > 1:
        # Ensure mpi is enebaled
        sqsub_args.extend(['-f', 'mpi'])
        # request nodes
        sqsub_args.extend(['-n', '%i' % nodes])
    # run-time estimate mandatory job type default?
    sqsub_args.extend(['-r', '6h'])
    # Memory might need increasing
    sqsub_args.extend(['--mpp=2.5g'])
    # Output
    sqsub_args.extend(['-o', 'faps-%s.out' % job_name])
    # Which command?
    sqsub_args.extend([exe])

    submit = Popen(sqsub_args, stdin=PIPE, stdout=PIPE)
    for line in submit.stdout.readlines():
        if 'submitted as' in line:
            jobid = int(line.split()[-1])
            break
    else:
        print("Job submission failed?")

    return jobid


def _sharcnet_postrun(waitfor_jobid):
    """Resubmit this script for the postrun on job completion."""
    sqsub_args = [
        'sqsub',
        '-q', 'DR_20293',
        '-r', '10m',
        '-o', 'faps-post-%s.out' % waitfor_jobid,
        '--waitfor=%s' % waitfor_jobid,
        ] + sys.argv
    print sqsub_args
    submit = Popen(sqsub_args, stdout=PIPE, stderr=STDOUT)
    # sqsub is slow sometimes, so we wait for it to complete
    submit.wait()
    print submit.stdout.readlines()

def _sharcnet_jobcheck(jobid):
    """Return true if job is still running or queued."""
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
        print("Failed to get job information.")  # qstat parsing failed?
        # Act as if the job is still running, in case it hasn't finished
        return True




def _wooki_generic(job_name, nodes=1, attributes=None):
    """Generic wooki submission to qsub."""

    max_cpus = 64  # Max that can be requested
    mem_default = None  # Passed to qsub
    default_queue = None  # leave blank ("") to use Wooki default (NOT attribute)
    shared_memory_only = False  # True if no infiniband version

    serial_scratch_dir = "/shared_scratch"
    parallel_scratch_dir = "/shared_scratch"

    vasp_versions = {
        "shm":
            {"mpi": "/home/system/LIBRARIES/mpich-intel10-shmem/bin/mpirun",
             "exe": "/home/program/CHEMISTRY/VASP/bin/vasp-5.2.11-ifort-shm"},
        "par":
            {"mpi": "/home/system/LIBRARIES/mvapich-intel/bin/mpirun",
             "exe": "/home/program/CHEMISTRY/VASP/bin/vasp-5.2.11-ifort-par"}}

    # Items in 'jobopts' are used to template the submit script
    job_opts = {}

    # Assume all files exist from calling script.
    job_opts["job_name"] = job_name
    job_opts["nodes"] = nodes

    # -----------------------------------------------------------------------------
    # As long as a specific attribute isn't selected as a command line option with
    # the '-a' option, then the script will attempt to run a parallel job in shared
    # memory mode first on the serial nodes, then on the parallel nodes. For a
    # single CPU job, if there are no serial nodes available, then it will attempt
    # to run on a parallel node
    # -----------------------------------------------------------------------------
    if not attributes is not None:
        freenodes_serial = int(Popen(['/home/program/bin/freenodes', '-m'],
                                     stdout=PIPE).communicate()[0])
        freenodes_parallel = int(Popen(['/home/program/bin/freenodes', '-p'],
                                       stdout=PIPE).communicate()[0])
        if nodes <= freenodes_serial:
#            print "Serial node is available to run in shared memory mode"
            attributes = "serial"
            scratch_dir = serial_scratch_dir
            shared_memory_job = True
        elif nodes <= freenodes_parallel:
#            print "Infiniband node available to run in shared memory mode"
            attributes = "infiniband"
            scratch_dir = parallel_scratch_dir
            shared_memory_job = False
        elif shared_memory_only == True:
        #------------------------------------------------------------------
        # If neither serial or parallel nodes are available and the job is a
        # shared memory only executable, then revert back to serial nodes
        #------------------------------------------------------------------
            attributes = "serial"
            scratch_dir = serial_scratch_dir
            shared_memory_job = True
        else:
            attributes = "infiniband"
            shared_memory_job = False
            print "Parallel execution on inifiband nodes "
            scratch_dir = parallel_scratch_dir
        if nodes == 1:
            scratch_dir = serial_scratch_dir
    else:
        # attributes specified
        if shared_memory_only == True or nodes == 1 or not "infini" in attributes:
            scratch_dir = serial_scratch_dir
            shared_memory_job = True
        else:
            scratch_dir = parallel_scratch_dir
            shared_memory_job = False

    job_opts["scratch"] = scratch_dir

    # ------------------------------------------------------
    # Make sure that user has a directory on /shared_scratch
    # ------------------------------------------------------
    job_dir = os.path.join("/shared_scratch/", getpass.getuser(), job_name)
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    #------------------------------------------------------
    # Construct PBS directives
    #------------------------------------------------------
    pbs_directives = ["#PBS -N fap-%s" % job_name,
                      "#PBS -m n",
                      "#PBS -o std.out",
                      "#PBS -j oe "]

    # ----------------------------------------------
    # -q --queue : specify a queue to run on
    # ----------------------------------------------
    #if options.queue_name:
        # optparse will set the default
    #    PBS_directives.append("#PBS -q %s" % options.queue_name)

    # ----------------------------------------------
    # -H --Hostname : specify the hostname to run on
    # ----------------------------------------------
    #if options.host_name:
    #    host_name = options.host_name
    #    print(">> Attempting to submit to host: %s" % host_name)
    #    print(" * -H/--Host option overrrides all default/specified attributes")

    #    PBS_directives.append("#PBS -l host=%s" % host_name)

        # -------------------------------------------------------------
        # Override other directives if host is specified
        # -------------------------------------------------------------
    #    options.min_mem = mem_default
    #    node_attribute = "serial"

    # -------------------------------------------------------------
    # -m --mem Memory options
    # -------------------------------------------------------------
    #if options.min_mem:
    #    print(">> Requesting node(s) with minimum memory of: %s" % options.min_mem)
    #    PBS_directives.append("#PBS -l mem=%s" % options.min_mem)

    # ----------------------------------------------
    # Specify Node numbers
    # ----------------------------------------------
    if nodes == 1:
        pbs_directives.append("#PBS -l nodes=1:%s" % attributes)
        try:
            job_opts.update(executable["ser"])
        except KeyError:
            try:
                job_opts.update(executable["shm"])
            except KeyError:
                print("no serial executable -- will now crash")
    elif shared_memory_job == True:
        pbs_directives.append("#PBS -l nodes=1:ppn=%i:%s" % (nodes, attributes))
        try:
            job_opts.update(executable["shm"])
        except KeyError:
            print("Shared memory version not available for %s")
    else:
        pbs_directives.append("#PBS -l nodes=%i:%s" % (nodes, attributes))
        try:
            job_opts.update(executable["par"])
        except KeyError:
            print("Parallel version not available for %s")



    #------------------------------------------------------
    # JOB type specific checks i.e. Gaussian specific
    #------------------------------------------------------

    # ------------------------------------------------------------------------
    # Define the bash script passed to qsub.  This can be cut out from the old
    # bash based submit scripts.
    #
    #  REMOVE the 'qsub << eof' line and all of the "#PBS " directive lines from the
    #  bash script.
    #    i.e.   remove th following from the script:
    #                qsub << eof
    #                #PBS -N %(jobname)s
    #                #PBS -l nodes=1:ppn=$NODES$ATTRIBUTE
    #                #PBS -m n
    #                 etc....
    #
    #  replace $QJOBNAME with %(jobname)s
    # ------------------------------------------------------------------------

    job_opts["startdir"] = os.getcwd()


    #******************************************************************************
    # QSUB bash script Begin
    #******************************************************************************
    line = """

    echo $HOSTNAME
    echo JOBID = $PBS_JOBID
    echo ==========================
    echo $PBS_NODEFILE
    cat $PBS_NODEFILE
    echo ==========================

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/program/LIBRARIES/CMKL/8.1/lib/em64t:/home/program/COMPILERS/Intel/fce/9.1.040/lib/:/home/program/COMPILERS/Intel/Compilers/11.0/069/mkl/lib/em64t/

    #---------------------------------------------------------------
    # make the run directory on /shared_scratch and a soft link to
    # the directory from the launch directory.
    #---------------------------------------------------------------
    if [ ! -e /shared_scratch/`whoami` ]; then
      mkdir /shared_scratch/`whoami`
    fi
    if [ ! -e /shared_scratch/`whoami`/%(job_name)s ]; then
      mkdir /shared_scratch/`whoami`/%(job_name)s
    fi
    if [ ! -e %(startdir)s/%(job_name)s.restart_DIR ]; then
      ln -s /shared_scratch/`whoami`/%(job_name)s %(startdir)s/%(job_name)s.restart_DIR
    fi
    cd /shared_scratch/`whoami`/%(job_name)s

    cp %(startdir)s/%(job_name)s.incar INCAR
    cp %(startdir)s/%(job_name)s.kpoints KPOINTS
    cp %(startdir)s/%(job_name)s.poscar POSCAR
    cp %(startdir)s/%(job_name)s.potcar POTCAR
    ln -sf %(startdir)s/%(job_name)s.outcar OUTCAR
    ln -sf %(startdir)s/%(job_name)s.contcar CONTCAR
    ln -sf %(startdir)s/%(job_name)s.log OSZICAR


    # ===============================
    # JOB EXECUTION
    # ===============================
    echo "Main Host:" $HOSTNAME >> VASP.stdout
    cat $PBS_NODEFILE | awk '{print $1".cluster";}' > machine_file
    echo "============================================" >> VASP.stdout
    echo "Contents of MPI machine_file" >> VASP.stdout
    echo " - This lists the nodes the job was run on" >> VASP.stdout
    echo "============================================" >> VASP.stdout
    cat machine_file >> VASP.stdout
    echo "============================================" >> VASP.stdout

    # vasp 5.2.2 will crash without this
    ulimit -s unlimited

    %(mpi)s -np %(nodes)i -machinefile machine_file %(exe)s >> VASP.stdout

    touch $PBS_JOBID

    rm machine_file

    """ % jobopts

    # TODO(tdaff): resubmit!
    #******************************************************************************
    # QSUB bash script END
    #******************************************************************************


    #-------------------------------------------------------------
    # Add the PBS directives to the qsub bash script and submit it
    #-------------------------------------------------------------
    line = "\n".join(pbs_directives) + "\n" + line
    submit = Popen("qsub", shell=False, stdin=PIPE)
    submit.communicate(input=line)


def _wooki_postrun(jobid):
    pbs_directives = ["#PBS -N fap-%s" % job_name,
                      "#PBS -m n",
                      "#PBS -o std.out",
                      "#PBS -j oe ",
                      "#PBS "
                      "cd $PBS_O_WORKDIR",
                      "python faps.py"]
