#!/usr/bin/env python

"""
faps -- Frontend for Automated Adsorption Analysis of Porous Solids.

aka

shpes -- Sorption analysis with a High throughput Python
         frontend to Examine binding in Structures

Strucutre adsorption property analysis for high throughput processing. Run
as a script, faps will automatically run complete analysis on a structure.
Sensible defaults are implemented, but calculations can be easily customised.
Faps also provides classes and methods for adapting the simulation or only
doing select parts.

"""

__version__ = "$Revision$"

import code
import ConfigParser
import glob
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
import textwrap
from copy import copy
from math import ceil
from logging import warning, debug, error, info, critical

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array, identity, dot, cross
from numpy.linalg import norm

from config import Options
from elements import WEIGHT, ATOMIC_NUMBER, UFF, VASP_PSEUDO_PREF
from job_handler import JobHandler
from logo import LOGO

# Global constants
DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859
EV2KCAL = 23.060542301389

# ID values for system state
NOT_RUN = 0
RUNNING = 1
FINISHED = 2
UPDATED = 3
SKIPPED = -1


class PyNiss(object):
    """
    PyNiss -- Negotiation of Intermediate System States

    A single property calculation for one structure. Instance with a set of
    options, then run the job_dispatcher() to begin the calculation. The
    calculation will pickle itself, or can be pickled at any time, by calling
    dump_state().

    """
    def __init__(self, options):
        """
        Instance an empty structure in the calculation; The dispatcher should
        be called to fill it up with data, as needed.

        """
        self.options = options
        self.structure = Structure(options.get('job_name'))
        self.state = {'init': (NOT_RUN, False),
                      'dft': (NOT_RUN, False),
                      'esp': (NOT_RUN, False),
                      'charges': (NOT_RUN, False),
                      'gcmc': {}}
        self.job_handler = JobHandler(options)

    def dump_state(self):
        """Write the .niss file holding the current system state."""
        job_name = self.options.get('job_name')
        info("Writing state file, %s.niss." % job_name)
        os.chdir(self.options.get('job_dir'))
        my_niss = open(job_name + ".niss", "wb")
        pickle.dump(self, my_niss)
        my_niss.close()

    def re_init(self, new_options):
        """Re initialize simulation (with updated options)."""
        if new_options.getbool('update_opts'):
            info("Using new options.")
            self.options = new_options
        else:
            # Just update command line stuff
            info("Using old options with new command line arguments.")
            self.options.args = new_options.args
            self.options.options = new_options.options
            self.options.cmdopts = new_options.cmdopts
        self.status(initial=True)

    def job_dispatcher(self):
        """
        Run parts explicity specified on the command line or do the next step
        in an automated run. Drop to interactive mode, if requested.

        """

        # In case options have changed, re-intitialize
        self.job_handler = JobHandler(self.options)

        if 'status' in self.options.args:
            self.status(initial=True)

        if self.options.getbool('interactive'):
            console = code.InteractiveConsole(locals())
            console.interact(
                banner="""See manual for instructions on interactive use.""")

        if self.options.getbool('import'):
            info("Importing results from a previous simulation")
            self.import_old()
            self.dump_state()

        if self.state['init'][0] == NOT_RUN:
            info("Reading in structure")
            # No structure, should get one
            self.structure.from_file(
                self.options.get('job_name'),
                self.options.get('initial_structure_format'),
                self.options)
            self.state['init'] = (UPDATED, False)
            self.dump_state()

        if self.state['dft'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_dft'):
                info("Skipping DFT step completely")
                warning("Job might fail later if you need the ESP")
                self.state['dft'] = (SKIPPED, False)
            elif self.state['dft'][0] == RUNNING:
                job_state = self.job_handler.jobcheck(self.state['dft'][1])
                if not job_state:
                    info("Queue reports DFT step has finished")
                    # Finished running update positions
                    self.structure.update_pos(self.options.get('dft_code'))
                    self.state['dft'] = (UPDATED, False)
                    self.dump_state()
                else:
                    # Still running
                    info("DFT still in progress")
                    terminate(0)

        if self.state['dft'][0] == NOT_RUN or 'dft' in self.options.args:
            self.run_dft()
            self.dump_state()
            terminate(0)

        if self.state['charges'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_charges'):
                info("Skipping charge calculation")
                self.state['charges'] = (SKIPPED, False)
            elif self.state['charges'][0] == RUNNING:
                job_state = self.job_handler.jobcheck(self.state['charges'][1])
                if not job_state:
                    info("Queue reports charge calculation has finished")
                    self.structure.update_charges(
                        self.options.get('charge_method'), self.options)
                    self.state['charges'] = (UPDATED, False)
                    self.dump_state()
                else:
                    info("Charge calculation still running")
                    terminate(0)

        if self.state['charges'][0] == NOT_RUN or 'charges' in self.options.args:
            self.run_charges()
            info("Running charge calculation")
            self.dump_state()
            terminate(0)

        if self.options.getbool('no_gcmc'):
            info("Skipping GCMC simulation")
        elif not self.state['gcmc'] or 'gcmc' in self.options.args:
            # The dictionary is empty before any runs
            info("Starting gcmc step")
            self.run_fastmc()
            self.dump_state()
            terminate(0)
        else:
            unfinished_gcmc = False
            for tp_point, tp_state in self.state['gcmc'].iteritems():
                if tp_state[0] == RUNNING:
                    new_state = self.job_handler.jobcheck(tp_state[1])
                    if not new_state:
                        info("Queue reports GCMC %s finished" % (tp_point,))
                        self.structure.update_gcmc(self.options.get('mc_code'), tp_point)
                        self.state['gcmc'][tp_point] = (UPDATED, False)
                        self.dump_state()
                    else:
                        info("GCMC %s still running" % (tp_point,))
                        unfinished_gcmc = True
            if not unfinished_gcmc:
                info("GCMC run has finished")
                self.post_summary()
                terminate(0)

    def status(self, initial=False):
        """Print the current status to the terminal."""
        valid_states = {NOT_RUN: 'Not run',
                        RUNNING: 'Running',
                        FINISHED: 'Finished',
                        UPDATED: 'Processed',
                        SKIPPED: 'Skipped'}

        if initial:
            info("Previous system state reported from .niss file "
                 "(running jobs may have already finished):")
        else:
            info("Current system status:")
        for step, state in self.state.iteritems():
            if step == 'gcmc':
                if not state:
                    info(" * State of GCMC: Not run")
                else:
                    for point, job in state.iteritems():
                        if job[0] is RUNNING:
                            info(" * GCMC %s: Running, jobid: %s" %
                                 (point, job[1]))
                        else:
                            info(" * GCMC %s: %s" %
                                 (point, valid_states[job[0]]))
            elif state[0] is RUNNING:
                info(" * State of %s: Running, jobid: %s" % (step, state[1]))
            else:
                info(" * State of %s: %s" % (step, valid_states[state[0]]))

    def post_summary(self):
        """Summarise any GCMC results."""
        info("Summary of GCMC results")
        nguests = len(self.structure.guests)
        for idx, guest in enumerate(self.structure.guests):
            csv = ["#T/K,p/bar,mol/uc,mmol/g,stdev,",
                   "v/v,stdev,hoa/kcal/mol,stdev,",
                   ",".join("p(g%i)" % gidx for gidx in range(nguests)), "\n"]
            info(guest.name)
            info(" mol/uc  mmol/g  vstp/v     hoa     T_P")
            for tp_point in sorted(guest.uptake):
                # <N>, sd, supercell
                uptake = guest.uptake[tp_point]
                hoa = guest.hoa[tp_point]
                vuptake = (guest.molar_volume*(uptake[0]/uptake[2])/
                           (6.023E-4*self.structure.volume))
                vuptake_stdev = (guest.molar_volume*(uptake[1]/uptake[2])/
                                 (6.023E-4*self.structure.volume))
                info("%7.2f %7.2f %7.2f %7.2f %s" % (
                    uptake[0]/uptake[2],
                    1000*uptake[0]/(uptake[2]*self.structure.weight),
                    vuptake,
                    hoa[0],
                    ("T=%s" % tp_point[0] +
                     ''.join(['P=%s' % x for x in tp_point[1]]))))
                csv.append("%f,%f,%f,%f,%f,%f,%f,%f,%f," % (
                    tp_point[0], tp_point[1][idx], uptake[0]/uptake[2],
                    1000*uptake[0]/(uptake[2]*self.structure.weight),
                    1000*uptake[1]/(uptake[2]*self.structure.weight),
                    vuptake, vuptake_stdev,
                    hoa[0], hoa[1]) +
                    ",".join("%f" % x for x in tp_point[1]) + "\n")
            csv_file = file('%s-%s.csv' %
                            (self.options.get('job_name'), guest.ident), 'wb')
            csv_file.writelines(csv)
            csv_file.close()

    def import_old(self):
        """Try and import any data from previous stopped simulation."""
        job_name = self.options.get('job_name')
        job_dir = self.options.get('job_dir')
        try:
            self.structure.from_file(
                job_name,
                self.options.get('initial_structure_format'))
            self.state['init'] = (UPDATED, False)
        except IOError:
            info("No initial structure found to import")
        try:
            self.structure.update_pos(self.options.get('dft_code'))
            self.state['dft'] = (UPDATED, False)
        except IOError:
            info("No optimized structure found to import")
        try:
            self.structure.update_charges(self.options.get('charge_method'))
            self.state['charges'] = (UPDATED, False)
        except IOError:
            info("No charges found to import")
        # Need to generate supercell here on import so that it is set, and
        # is based on the cell from dft, if changed
        self.structure.gen_supercell(self.options)

        guests = [Guest(x) for x in self.options.gettuple('guests')]
        if not same_guests(self.structure.guests, guests):
            info("Replacing old guests")
            self.structure.guests = guests
        else:
            # use existing guests that may have data
            debug("Retaining previous guests")
            guests = self.structure.guests
        temps = self.options.gettuple('mc_temperature', float)
        presses = self.options.gettuple('mc_pressure', float)
        indivs = self.options.gettuple('mc_state_points', float)
        for tp_point in state_points(temps, presses, indivs, len(guests)):
            try:
                self.structure.update_gcmc(self.options.get('mc_code'), tp_point)
                self.state['gcmc'][tp_point] = (UPDATED, False)
            except IOError:
                info("GCMC point %s not found" % str(tp_point))

        # Reset directory at end
        os.chdir(job_dir)

    def run_dft(self):
        """Select correct method for running the dft/optim."""
        dft_code = self.options.get('dft_code')
        info("Running a %s calculation" % dft_code)
        if dft_code == 'vasp':
            self.run_vasp()
        elif dft_code == 'siesta':
            self.run_siesta()
        else:
            critical("Unknown dft method: %s" % dft_code)
            terminate(92)

    def run_charges(self):
        """Select correct charge processing methods."""
        chg_method = self.options.get('charge_method')
        info("Calculating charges with %s" % chg_method)
        if chg_method == 'repeat':
            # Get ESP
            self.esp_to_cube()
            self.run_repeat()
        elif chg_method == 'gulp':
            self.run_qeq_gulp()
        else:
            critical("Unknown charge calculation method: %s" % chg_method)
            terminate(93)

    def run_vasp(self):
        """Make inputs and run vasp job."""
        job_name = self.options.get('job_name')
        nproc = self.options.getint('vasp_ncpu')
        # Keep things tidy in a subdirectory
        dft_code = self.options.get('dft_code')
        vasp_dir = os.path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, dft_code))
        mkdirs(vasp_dir)
        os.chdir(vasp_dir)
        debug("Running in %s" % vasp_dir)
        info("Running on %i nodes" % nproc)

        filetemp = open("POSCAR", "wb")
        filetemp.writelines(self.structure.to_vasp(self.options))
        filetemp.close()

        esp_grid = self.esp_grid

        filetemp = open("INCAR", "wb")
        filetemp.writelines(mk_incar(self.options, esp_grid=esp_grid))
        filetemp.close()

        filetemp = open("KPOINTS", "wb")
        filetemp.writelines(mk_kpoints(self.options.gettuple('kpoints', int)))
        filetemp.close()

        potcar_types = unique(self.structure.types)
        filetemp = open("POTCAR", "wb")
        potcar_dir = self.options.get('potcar_dir')
        for at_type in potcar_types:
            # Try and get the preferred POTCARS
            debug("Using %s pseudopotential for %s" %
                 (VASP_PSEUDO_PREF.get(at_type, at_type), at_type))
            potcar_src = os.path.join(potcar_dir,
                                      VASP_PSEUDO_PREF.get(at_type, at_type),
                                      "POTCAR")
            shutil.copyfileobj(open(potcar_src), filetemp)
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("Vasp input files generated; skipping job submission")
            self.state['dft'] = (SKIPPED, False)
        else:
            self.job_handler.env(dft_code, options=self.options)
            jobid = self.job_handler.submit(dft_code, self.options)
            info("Running VASP job in queue. Jobid: %s" % jobid)
            self.state['dft'] = (RUNNING, jobid)
            if self.options.getbool('run_all'):
                debug('Submitting postrun script')
                os.chdir(self.options.get('job_dir'))
                self.job_handler.postrun(jobid)
            else:
                debug('Postrun script not submitted')
        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))

    def run_siesta(self):
        """Make siesta input and run job."""
        job_name = self.options.get('job_name')
        nproc = self.options.getint('siesta_ncpu')
        # Keep things tidy in a subdirectory
        dft_code = self.options.get('dft_code')
        siesta_dir = os.path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, dft_code))
        mkdirs(siesta_dir)
        os.chdir(siesta_dir)
        debug("Running in %s" % siesta_dir)
        info("Running on %i nodes" % nproc)

        filetemp = open('%s.fdf' % job_name, 'wb')
        filetemp.writelines(self.structure.to_siesta(self.options))
        filetemp.close()

        psf_types = unique(self.structure.types)
        psf_dir = self.options.get('psf_dir')
        for at_type in psf_types:
            psf_atm = '%s.psf' % at_type
            psf_src = os.path.join(psf_dir, psf_atm)
            psf_dest = os.path.join(siesta_dir, psf_atm)
            try:
                if not os.path.exists(psf_atm):
                    os.symlink(psf_src, psf_dest)
            # symlinks not available pre 3.2 on windows
            except AttributeError:
                shutil.copy(psf_src, siesta_dir)
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("Siesta input files generated; skipping job submission")
            self.state['dft'] = (SKIPPED, False)
        else:
            # sharcnet does weird things for siesta
            self.job_handler.env(dft_code, options=self.options)
            jobid = self.job_handler.submit(dft_code, self.options,
                                            input_file='%s.fdf' % job_name)
            info("Running SIESTA job in queue. Jobid: %s" % jobid)
            self.state['dft'] = (RUNNING, jobid)
            if self.options.getbool('run_all'):
                debug('Submitting postrun script')
                os.chdir(self.options.get('job_dir'))
                self.job_handler.postrun(jobid)
            else:
                debug('Postrun script not submitted')
        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))

    def run_qeq_gulp(self):
        """Run GULP to calculate charge equilibration charges."""
        job_name = self.options.get('job_name')
        qeq_code = 'gulp'
        qeq_dir = os.path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, qeq_code))
        mkdirs(qeq_dir)
        os.chdir(qeq_dir)
        debug("Running in %s" % qeq_dir)

        filetemp = open('%s.gin' % job_name, 'wb')
        filetemp.writelines(self.structure.to_gulp())
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("GULP input files generated; skipping job submission")
            self.state['charges'] = (SKIPPED, False)
        else:
            jobid = self.job_handler.submit(qeq_code, self.options,
                                            input_file='%s.gin' % job_name)
            info("Running GULP job in queue. Jobid: %s" % jobid)
            self.state['charges'] = (RUNNING, jobid)
            if self.options.getbool('run_all'):
                debug('Submitting postrun script')
                os.chdir(self.options.get('job_dir'))
                self.job_handler.postrun(jobid)
            else:
                debug('Postrun script not submitted')
        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))

    def esp_to_cube(self):
        """Make the cube for repeat input."""
        job_name = self.options.get('job_name')
        esp_src = self.options.get('esp_src')
        repeat_dir = os.path.join(self.options.get('job_dir'),
                                  'faps_%s_repeat' % job_name)
        mkdirs(repeat_dir)
        src_dir = os.path.join(self.options.get('job_dir'),
                               'faps_%s_%s' % (job_name, esp_src))
        os.chdir(src_dir)
        if esp_src == 'vasp':
            esp_to_cube_args = shlex.split(self.options.get('vasp_to_cube'))
            info("Converting vasp esp to cube, this might take a minute...")
            submit = subprocess.Popen(esp_to_cube_args)
            submit.wait()
            # Cube should have job_name, but can get truncated;
            # therefore we try to look for it first
            cube_file = glob.glob('*.cube')
            if len(cube_file) == 1:
                cube_file = cube_file[0]
            elif len(cube_file) > 1:
                cube_file = cube_file[0]
                warning("More or than one .cube found; using %s" % cube_file)
            else:
                critical("No cube files found; check vasp_to_cube output")
                terminate(101)
            # Move it to the repeat directory and give a proper name
            move_and_overwrite(cube_file,
                               os.path.join(repeat_dir, job_name + '.cube'))
            unneeded_files = self.options.gettuple('vasp_delete_files')
            remove_files(unneeded_files)
            keep_files = self.options.gettuple('vasp_compress_files')
            compress_files(keep_files)
        elif esp_src == 'siesta':
            esp_to_cube_args = shlex.split(self.options.get('siesta_to_cube'))
            esp_grid = self.esp_grid
            info("Generating ESP grid of %ix%ix%i" % esp_grid)
            siesta_to_cube_input = [
                "%s\n" % job_name,
                "%f %f %f\n" % (0.0, 0.0, 0.0),
                "%i %i %i\n" % esp_grid]
            info("Converting siesta esp to cube, this might take a minute...")
            submit = subprocess.Popen(esp_to_cube_args, stdin=subprocess.PIPE)
            submit.communicate(input=''.join(siesta_to_cube_input))
            move_and_overwrite(job_name + '.cube', repeat_dir)
            unneeded_files = self.options.gettuple('siesta_delete_files')
            remove_files(unneeded_files)
            keep_files = self.options.gettuple('siesta_compress_files')
            compress_files(keep_files)

        os.chdir(self.options.get('job_dir'))

    def run_repeat(self):
        """Submit the repeat calc to the queue."""
        job_name = self.options.get('job_name')
        charge_code = self.options.get('charge_method')
        repeat_dir = os.path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, charge_code))
        mkdirs(repeat_dir)
        os.chdir(repeat_dir)

        mk_repeat(cube_name=job_name + '.cube')
        if self.options.getbool('no_submit'):
            info("REPEAT input files generated; skipping job submission")
            self.state['charges'] = (SKIPPED, False)
        else:
            jobid = self.job_handler.submit(charge_code, self.options)
            info("Running REPEAT calculation in queue: Jobid %s" % jobid)
            self.state['charges'] = (RUNNING, jobid)
            if self.options.getbool('run_all'):
                debug('Submitting postrun script')
                os.chdir(self.options.get('job_dir'))
                self.job_handler.postrun(jobid)
            else:
                debug('Postrun script not submitted')

        os.chdir(self.options.get('job_dir'))

    def run_fastmc(self):
        """Submit a fastmc job to the queue."""
        job_name = self.options.get('job_name')
        mc_code = self.options.get('mc_code')

        # Set the guests before generating the files
        # Load here as options may change in each run
        # and before changing directory, or it will not find guests.lib
        guests = [Guest(x) for x in self.options.gettuple('guests')]
        if not same_guests(self.structure.guests, guests):
            info("Replacing old guests")
            self.structure.guests = guests
        else:
            # use existing guests that may have data
            debug("Retaining previous guests")
            guests = self.structure.guests

        gcmc_dir = os.path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, mc_code))
        mkdirs(gcmc_dir)
        os.chdir(gcmc_dir)

        config, field = self.structure.to_fastmc(self.options)

        filetemp = open("CONFIG", "wb")
        filetemp.writelines(config)
        filetemp.close()

        filetemp = open("FIELD", "wb")
        filetemp.writelines(field)
        filetemp.close()

        temps = self.options.gettuple('mc_temperature', float)
        presses = self.options.gettuple('mc_pressure', float)
        indivs = self.options.gettuple('mc_state_points', float)
        jobids = []
        for tp_point in state_points(temps, presses, indivs, len(guests)):
            temp = tp_point[0]
            press = tp_point[1]
            info("Running GCMC: T=%.1f " % temp +
                 " ".join(["P=%.2f" % x for x in press]))
            tp_path = ('T%s' % temp +
                       ''.join(['P%.2f' % x for x in press]))
            mkdirs(tp_path)
            try_symlink(os.path.join(gcmc_dir, 'CONFIG'),
                        os.path.join(tp_path, 'CONFIG'))
            try_symlink(os.path.join(gcmc_dir, 'FIELD'),
                        os.path.join(tp_path, 'FIELD'))
            os.chdir(tp_path)
            filetemp = open("CONTROL", "wb")
            filetemp.writelines(mk_gcmc_control(temp, press, self.options,
                                                guests))
            filetemp.close()

            if self.options.getbool('no_submit'):
                info("FastMC input files generated; "
                     "skipping job submission")
                self.state['gcmc'][(temp, press)] = (SKIPPED, False)
            else:
                jobid = self.job_handler.submit(mc_code, self.options)
                info("Running FastMC in queue: Jobid %s" % jobid)
                self.state['gcmc'][(temp, press)] = (RUNNING, jobid)
                jobids.append(jobid)
            os.chdir('..')

        # Postrun after all submitted so don't have to deal with messy
        # directory switching
        os.chdir(self.options.get('job_dir'))
        # jobids will be empty if nothing has been submitted
        if self.options.getbool('run_all') and jobids:
            debug('Submitting postrun script')
            self.job_handler.postrun(jobids)
        else:
            debug('Postrun script not submitted')

    @property
    def esp_grid(self):
        """Estimate the esp grid based on resolution and memory."""
        # If repeat is unsing double precision, use 4 for single
        repeat_prec = 8
        # User defined resolution, try to use this
        resolution = self.options.getfloat('esp_resolution')
        repeat_ncpu = self.options.getint('repeat_ncpu')
        if repeat_ncpu == 1:
            vmem = self.options.getfloat('serial_memory')
        else:
            vmem = self.options.getfloat('threaded_memory')
        # Nice even grids might scale better in parallel repeat
        esp_grid = tuple([int(4*np.ceil(x/(4*resolution)))
                          for x in self.structure.cell.params[:3]])
        memory_guess = prod(esp_grid)*self.structure.natoms*repeat_prec/1e9
        if memory_guess > vmem:
            warning("ESP at this resolution might need up to %.1f GB of "
                    "memory but calculation will only request %.1f" %
                    (memory_guess, vmem))
            resolution = resolution/pow(vmem/memory_guess, 1.0/3)
            esp_grid = tuple([int(4*np.ceil(x/(4*resolution)))
                              for x in self.structure.cell.params[:3]])
            warning("Reduced grid to %.2f A resolution to fit" % resolution)

        return esp_grid


class Structure(object):
    """
    The current state of the structure; update as the calculations proceed.

    Structure provides methods to produce input files for and take output from
    various computational chemistry packages but needs to be told what to do.
    Internal energy units are kcal/mol.

    Methods are grouped:
    * Initial structure parsers
    * Output file parsers to update structure
    * Input file generation
    * Internal manipulation methods

    """
    # FIXME: symmetry not considered
    # TODO: dft energy?
    def __init__(self, name):
        """Just instance an empty structure initially."""
        self.name = name
        self.cell = Cell()
        self.atoms = []
        self.esp = None
        self.dft_energy = 0.0
        self.guests = []
        self.properties = {}

    def from_file(self, basename, filetype, defaults):
        """Select the correct file parser."""
        filetype = filetype.lstrip('.')
        if filetype.lower() in ['pdb']:
            self.from_pdb(basename + '.' + filetype)
        elif filetype.lower() in ['vasp', 'poscar', 'contcar']:
            listdir = os.listdir('.')
            test_files = [
                basename + '.contcar', basename + '.CONTCAR', 'CONTCAR',
                basename + '.poscar', basename + '.POSCAR', 'POSCAR']
            for filename in test_files:
                if filename in listdir:
                    self.from_vasp(filename)
                    break
        elif filetype.lower() in ['cif']:
            self.from_cif(basename + '.' + filetype)
        elif filetype.lower() in ['xyz']:
            cell = defaults.gettuple('default_cell', float)
            self.from_xyz(basename + '.' + filetype, cell=cell)
        else:
            error("Unknown filetype %s" % filetype)

    def update_pos(self, opt_code):
        """Select the method for updating atomic positions."""
        opt_path = os.path.join('faps_%s_%s' % (self.name, opt_code))
        info("Updating positions from %s" % opt_code)
        if opt_code == 'vasp':
            self.from_vasp(os.path.join(opt_path, 'CONTCAR'), update=True)
        elif opt_code == 'siesta':
            self.from_siesta(os.path.join(opt_path, '%s.STRUCT_OUT' % self.name))
        else:
            error("Unknown positions to import %s" % opt_code)

    def update_charges(self, charge_method, options=None):
        """Select the method for updating charges."""
        charge_path = os.path.join('faps_%s_%s' % (self.name, charge_method))
        if charge_method == 'repeat':
            info("Updating charges from repeat")
            self.charges_from_repeat(
                os.path.join(charge_path, 'faps-%s.out' % self.name))
            # Cleanup of REPEAT files
            unneeded_files = options.gettuple('repeat_delete_files')
            remove_files(unneeded_files, charge_path)
            keep_files = options.gettuple('repeat_compress_files')
            compress_files(keep_files, charge_path)
        elif charge_method == 'gulp':
            info("Updating charges from GULP QEq")
            self.charges_from_gulp(
                os.path.join(charge_path, 'faps-%s.out' % self.name))
        else:
            error("Unknown charge method to import %s" % charge_method)

    def update_gcmc(self, gcmc_code, tp_point):
        """Select the source for GCMC results and import."""
        gcmc_path = os.path.join('faps_%s_%s' % (self.name, gcmc_code))
        # Runs in subdirectories
        tp_path = ('T%s' % tp_point[0] +
                   ''.join(['P%.2f' % x for x in tp_point[1]]))
        if gcmc_code == 'fastmc':
            info("Importing results from FastGCMC")
            self.fastmc_postproc(
                os.path.join(gcmc_path, tp_path, 'OUTPUT'), tp_point)
        else:
            error("Unknown gcmc method to import %s" % gcmc_code)

    def from_pdb(self, filename):
        """Read an initial structure from a pdb file."""
        info("Reading positions from pdb file: %s" % filename)
        filetemp = open(filename)
        pdb_file = filetemp.readlines()
        filetemp.close()
        # Build a local list before setting attribute
        newatoms = []
        for line in pdb_file:
            if line.lower().startswith('cryst1'):
                self.cell.from_pdb(line)
            elif line.lower().startswith('atom'):
                newatom = Atom()
                newatom.from_pdb(line)
                newatoms.append(newatom)
            elif line.lower().startswith('hetatm'):
                newatom = Atom()
                newatom.from_pdb(line)
                newatoms.append(newatom)

        self.atoms = newatoms
        self.order_by_types()

    def from_cif(self, filename):
        """Genereate structure from a .cif file."""
        info("Reading positions from cif file: %s" % filename)
        filetemp = open(filename)
        cif_file = filetemp.readlines()
        filetemp.close()
        cif_file = strip_blanks(cif_file)
        params = [None, None, None, None, None, None]
        atoms = []
        symmetry = []
        loops = []
        idx = 0
        while idx < len(cif_file):
            line = cif_file[idx].lower()
            if '_cell_length_a' in line:
                params[0] = ufloat(line.split()[1])
            elif '_cell_length_b' in line:
                params[1] = ufloat(line.split()[1])
            elif '_cell_length_c' in line:
                params[2] = ufloat(line.split()[1])
            elif '_cell_angle_alpha' in line:
                params[3] = ufloat(line.split()[1])
            elif '_cell_angle_beta' in line:
                params[4] = ufloat(line.split()[1])
            elif '_cell_angle_gamma' in line:
                params[5] = ufloat(line.split()[1])
            elif 'loop_' in line:
                # loops for _atom_site and _symmetry
                heads = []
                body = []
                while '_' in line:
                    heads.extend(line.split())
                    idx += 1
                    line = cif_file[idx]
                while idx < len(cif_file) and '_' not in line:
                    # shlex keeps 'quoted items' as one
                    body.extend(shlex.split(line))
                    idx += 1
                    try:
                        line = cif_file[idx]
                    except IndexError:
                        line = ''
                if 'loop_' in heads:
                    heads.remove('loop_')
                loops.append((heads, body))
                continue
            idx += 1

        # cell first
        if np.all(params):
            self.cell.params = params
        else:
            error("No cell or incomplete cell found in cif file")

        # parse loop contents
        for heads, body in loops:
            if '_atom_site_fract_x' in heads:
                while body:
                    atoms.append(dict(zip(heads, body)))
                    body = body[len(heads):]
            if '_symmetry_equiv_pos_as_xyz' in heads:
                while body:
                    sym_dict = dict(zip(heads, body))
                    symmetry.append(
                        Symmetry(sym_dict['_symmetry_equiv_pos_as_xyz']))
                    body = body[len(heads):]

        if not symmetry:
            debug('No symmetry found; assuming identity')
            symmetry = [Symmetry('x,y,z')]

        newatoms = []
        for atom in atoms:
            for sym_op in symmetry:
                newatom = Atom()
                newatom.from_cif(atom, self.cell.cell, sym_op)
                newatoms.append(newatom)

        self.atoms = newatoms
        if len(symmetry) > 1:
            # can skip if just identity operation as it's slow for big systems
            self.remove_duplicates()
        self.order_by_types()

    def from_vasp(self, filename='CONTCAR', update=False):
        """Read a structure from a vasp [POS,CONT]CAR file."""
        #TODO(tdaff): difference between initial and update?
        info("Reading positions from vasp file: %s" % filename)
        filetemp = open(filename)
        contcar = filetemp.readlines()
        filetemp.close()
        atom_list = []
        atom_types = []
        scale = float(contcar[1])
        self.cell.from_lines(contcar[2:5], scale)
        if contcar[5].split()[0].isalpha():
            # vasp 5 with atom names
            atom_types = contcar[5].split()
            del contcar[5]
        poscar_counts = [int(x) for x in contcar[5].split()]
        natoms = sum(poscar_counts)
        if contcar[6].strip()[0].lower() in "s":
            # 's'elective dynamics line; we don't care
            del contcar[6]

        # mcell converts frac -> cart if necessary and scales
        if contcar[6].strip()[0].lower() in "ck":
            mcell = identity(3) * scale
        else:
            mcell = self.cell.cell

        # parsing positions
        if update:
            for atom, at_line in zip(self.atoms, contcar[7:7+natoms]):
                atom.from_vasp(at_line, cell=mcell)
        elif not atom_types:
            critical("Will not extract structure from older vasp files")
        else:
            line_idx = 6
            for at_type, at_count in zip(atom_types, poscar_counts):
                for _atom_idx in range(at_count):
                    line_idx += 1
                    this_atom = Atom()
                    this_atom.from_vasp(contcar[line_idx], at_type, mcell)
                    atom_list.append(this_atom)
            self.atoms = atom_list
            self.order_by_types()

    def from_siesta(self, filename):
        """Update the structure from the siesta output."""
        info("Updating positions from file: %s" % filename)
        filetemp = open(filename)
        struct_out = filetemp.readlines()
        filetemp.close()
        self.cell.from_lines(struct_out[:3])
        for atom, line in zip(self.atoms, struct_out[4:]):
            atom.from_siesta(line, self.cell.cell)

    def from_xyz(self, filename, update=False, cell=None):
        """Read a structure from an file."""
        info("Reading positions from xyz file: %s" % filename)
        filetemp = open(filename)
        xyz_file = filetemp.readlines()
        filetemp.close()
        # Setting the cell
        if len(cell) == 6:
            self.cell.params = cell
        elif len(cell) == 9:
            self.cell.cell = array(cell).reshape((3, 3))
        elif cell is not None:
            error("Invalid cell specification %s" % str(cell))
        # Build a local list before setting attribute
        newatoms = []
        natoms = int(xyz_file[0])
        self.properties['header'] = xyz_file[1].strip()
        for line in xyz_file[2:2+natoms]:
            newatom = Atom()
            newatom.from_xyz(line)
            newatoms.append(newatom)
        if update:
            if natoms != self.natoms:
                critical("Incorrect number of atoms to update")
                terminate(96)
            for atom, newatom in zip(self.atoms, newatoms):
                if atom.type != newatom.type:
                    error("Atom order may have changed")
                atom.pos = newatom.pos
        else:
            self.atoms = newatoms
            self.order_by_types()

    def charges_from_repeat(self, filename):
        """Parse charges and update structure."""
        info("Getting charges from file: %s" % filename)
        charges = []
        filetemp = open(filename)
        for line in filetemp:
            if line.startswith(" Charge"):
                line = line.split()
                # index, type, charge
                charges.append((int(line[1]), int(line[4]), float(line[6])))
            if "Error" in line:
                if float(line.split()[-1]) > 0.6:
                    warning("Error in repeat charges is very high - check cube!")
        filetemp.close()
        # TODO(tdaff): no symmetry yet
        if len(charges) != len(self.atoms):
            error("Incorrect number of charges; check REPEAT output")
            terminate(90)
        for atom, charge in zip(self.atoms, charges):
            atom.charge = charge[2]

    def charges_from_gulp(self, filename):
        """Parse QEq charges from GULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        start_line = gout.index('  Final charges from QEq :\n') + 7
        for atom, chg_line in zip(self.atoms, gout[start_line:]):
            atom.charge = float(chg_line.split()[2])

    def to_vasp(self, options):
        """Return a vasp5 poscar as a list of lines."""
        optim_h = options.getbool('optim_h')
        optim_all = options.getbool('optim_all')
        types = [atom.type for atom in self.atoms]
        ordered_types = unique(types)
        poscar = ["%s\n" % self.name[:80],
                  " 1.0\n"]
        # Vasp does 16 dp but we get rounding errors (eg cubic) use 14
        poscar.extend(self.cell.to_vector_strings(fmt="%23.14f"))
        poscar.append("".join("%5s" % x for x in ordered_types) + "\n")
        poscar.append("".join("%6i" % types.count(x)
                              for x in ordered_types) + "\n")
        # We always have the T or F so turn on selective dynamics for
        # fixed pos variable cell
        poscar.extend(["Selective dynamics\n", "Cartesian\n"])
        if optim_all:
            info("Optimizing all atom positions")
            fix_h = 'T'
            fix_all = 'T'
        elif optim_h:
            info("Optimizing hydrogen positions")
            fix_h = 'T'
            fix_all = 'F'
        else:
            info("All atom positions fixed")
            fix_h = 'F'
            fix_all = 'F'

        # assume atoms are ordered
        # Atom lines differ from vasp to ensure spaces between numbers
        # with <-10 values
        for atom in self.atoms:
            if atom.type == "H":
                poscar.append("%20.15f %19.15f %19.15f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_h, fix_h, fix_h))
            else:
                poscar.append("%20.15f %19.15f %19.15f" % tuple(atom.pos) +
                              "%4s%4s%4s\n" % (fix_all, fix_all, fix_all))
        return poscar

    def to_siesta(self, options):
        """Return a siesta input file as a list of lines."""
        job_name = options.get('job_name')
        siesta_accuracy = options.get("siesta_accuracy").lower()
        if siesta_accuracy in ['high']:
            info("Using 'high' accuracy siesta settings")
            basis = ('DZP', 100, 200)
        elif siesta_accuracy in ['low']:
            info("Using 'low' accuracy siesta settings")
            basis = ('SZ', 200, 100)
        else:
            info("Using default siesta accuracy settings")
            basis = ('DZ', 150, 150)

        u_atoms = unique(self.atoms, key=lambda x: x.type)
        u_types = unique(self.types)
        fdf = ([
            "SystemName %s\n" % job_name,
            "SystemLabel %s\n" % job_name,
            "\n",
            "NumberOfAtoms %i\n" % len(self.atoms),
            "NumberOfSpecies %i\n" % len(u_atoms),
            "\n",
            "SaveElectrostaticPotential .true.\n",
            "WriteMullikenPop 1\n"
            "WriteXML F\n",
            "\n",
            "# Accuracy bits\n",
            "\n",
            "PAO.BasisSize %s\n" % basis[0],
            "PAO.EnergyShift %i meV\n" % basis[1],
            "MeshCutoff %i Ry\n" % basis[2],
            "\n",
            "MaxSCFIterations 100\n",
            "XC.Functional GGA\n",
            "XC.Authors PBE\n",
            "SolutionMethod diagon\n",
            "ElectronicTemperature 25 K\n",
            "DM.NumberPulay 5\n",
            "DM.MixingWeight 0.05\n",
            "\n",
            "%block ChemicalSpeciesLabel\n"] + [
            "%6i %6i %6s\n" % ((idx + 1), atom.atomic_number, atom.type)
                for idx, atom in enumerate(u_atoms)] + [
            "%endblock ChemicalSpeciesLabel\n",
            "\n",
            "LatticeConstant 1 Ang\n"
            "%block LatticeVectors\n"] +
            self.cell.to_vector_strings() + [
            "%endblock LatticeVectors\n",
            "\n",
            "AtomicCoordinatesFormat   Ang\n",
            "\n",
            "%block AtomicCoordinatesAndAtomicSpecies\n"] + [
            "%12.8f %12.8f %12.8f" % tuple(atom.pos) + "%6i\n" %
                (u_types.index(atom.type) + 1) for atom in self.atoms] + [
            "%endblock AtomicCoordinatesAndAtomicSpecies\n"])
        if "Br" in u_types:
            fdf.extend(["\n%block PS.lmax\n",
                        "   Br 2\n",
                        "%endblock PS.lmax\n"])

        optim_h = options.getbool('optim_h')
        optim_all = options.getbool('optim_all')
        optim_cell = options.getbool('optim_cell')

        constraint = []

        if optim_all:
            info("Optimizing all atom positions")
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 800)
        elif optim_h and "H" in self.types:
            info("Optimizing hydrogen positions")
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 300)
            constraint = ["%i" % (idx+1)
                          for idx, species in enumerate(self.types)
                          if species not in ["H"]]
        elif optim_cell:
            fdf.append("\nMD.TypeOfRun  CG\n")
            fdf.append("MD.NumCGSteps  %i\n" % 300)
            constraint = ["%i" % (idx+1) for idx in range(self.natoms)]

        if optim_cell:
            info("Cell vectors will be optimized")
            fdf.append("MD.VariableCell .true.\n")
        if constraint:
            constraint = textwrap.fill(" ".join(constraint),
                                       initial_indent='position ',
                                       subsequent_indent='position ')
            fdf.extend(["\n%block GeometryConstraints\n",
                        constraint, "\n",
                        "%endblock GeometryConstraints\n"])

        return fdf

    def to_gulp(self):
        """Return a GULP file to use for the QEq charges."""
        gin_file = [
            "single conp qeq\n\n",
            "name %s\n\n" % self.name,
            "vectors\n"] + self.cell.to_vector_strings() + [
            "cartesian\n"]
        for atom in self.atoms:
            gin_file.extend(["%s core " % atom.type,
                             "%f %f %f\n" % tuple(atom.pos)])
        gin_file.append("\n\nprint 1\n\n")
        return gin_file

    def to_fastmc(self, options):
        """Return the FIELD and CONFIG needed for a fastmc run."""
        # CONFIG
        self.gen_supercell(options)
        supercell = self.gcmc_supercell
        levcfg = 0  # always
        imcon = self.cell.imcon
        natoms = len(self.atoms) * prod(supercell)
        config = ["%s\n" % self.name[:80],
                  "%10i%10i%10i\n" % (levcfg, imcon, natoms)]
        config.extend(self.cell.to_vector_strings(scale=supercell))
        for idx, atom in enumerate(self.supercell(supercell)):
            # idx+1 for 1 based indexes in CONFIG
            config.extend(["%-6s%10i\n" % (atom.type, idx + 1),
                           "%20.12f%20.12f%20.12f\n" % tuple(atom.pos)])

        # FIELD
        ntypes = len(self.guests) + 1
        field = ["%s\n" % self.name[:80],
                 "UNITS   kcal\n",
                 "molecular types %i\n" % ntypes]
        # Guests
        for guest in self.guests:
            field.extend(["&guest %s: %s\n" % (guest.name, guest.source),
                          "NUMMOLS %i\n" % 0,
                          "ATOMS %i\n" % len(guest.atoms)])
            for atom in guest.atoms:
                field.append(("%-6s %12.6f %12.6f" %
                              tuple([atom.type, atom.mass, atom.charge])) +
                             ("%12.6f %12.6f %12.6f\n" % tuple(atom.pos)))
            field.append("finish\n")
        # Framework
        field.extend(["Framework\n",
                      "NUMMOLS %i\n" % prod(supercell),
                      "ATOMS %i\n" % len(self.atoms)])
        for atom in self.atoms:
            field.append("%-6s %12.6f %20.14f %6i %6i\n" %
                         (atom.type, atom.mass, atom.charge, 1, 1))
        field.append("finish\n")
        # VDW potentials
        atom_set = [atom.type for atom in self.atoms]
        for guest in self.guests:
            atom_set.extend(atom.type for atom in guest.atoms)
        atom_set = unique(atom_set)
        field.append("VDW %i\n" % ((len(atom_set) * (len(atom_set) + 1)) / 2))

        # modify local ff to deal with guests
        force_field = copy(UFF)
        for guest in self.guests:
            force_field.update(guest.potentials)

        for idxl in range(len(atom_set)):
            for idxr in range(idxl, len(atom_set)):
                left = atom_set[idxl]
                right = atom_set[idxr]
                try:
                    sigma, epsilon = lorentz_berthelot(force_field[left],
                                                       force_field[right])
                except KeyError:
                    # catch this if not in the UFF -> zero
                    warning("No potential defined for %s %s; defaulting to 0" %
                         (left, right))
                    sigma, epsilon = 0.0, 0.0
                field.append("%-6s %-6s lj %f %f\n" %
                             (left, right, epsilon, sigma))
        # EOF
        field.append("close\n")

        return config, field

    def fastmc_postproc(self, filename, tp_point):
        """Update structure properties from gcmc OUTPUT."""
        filetemp = open(filename)
        output = filetemp.readlines()
        filetemp.close()

        supercell_mult = prod(self.gcmc_supercell)
        for idx, line in enumerate(output):
            if line.startswith('   final stats'):
                guest_id = int(line.split()[4]) - 1
                self.guests[guest_id].uptake[tp_point] = (
                    float(output[idx + 3].split()[-1]),
                    float(output[idx + 4].split()[-1]),
                    supercell_mult)
                # This will sometimes be NaN
                self.guests[guest_id].hoa[tp_point] = (
                    float(output[idx + 7].split()[-1]),
                    float(output[idx + 8].split()[-1]))

    def remove_duplicates(self, epsilon=0.0002):
        """Find overlapping atoms and remove them."""
        uniq_atoms = []
        found_atoms = []
        cell = self.cell.cell
        for atom in self.atoms:
            ifpos = atom.ifpos(cell)
            for found_atom in found_atoms:
                if frac_near(ifpos, found_atom, epsilon=epsilon):
                    break
            # else excutes when not found here
            else:
                uniq_atoms.append(atom)
                found_atoms.append(ifpos)
        debug("Found %i unique atoms in %i" % (len(uniq_atoms), self.natoms))
        self.atoms = uniq_atoms

    def gen_supercell(self, options):
        """Cacluate the smallest satisfactory supercell and set attribute."""
        config_supercell = options.gettuple('mc_supercell', int)
        config_cutoff = options.getfloat('mc_cutoff')
        if config_cutoff < 12:
            warning("Simulation is using a very small cutoff! I hope you "
                 "know, what you are doing!")
        minimum_supercell = self.cell.minimum_supercell(config_cutoff)
        supercell = tuple(max(i, j)
                          for i, j in zip(config_supercell, minimum_supercell))
        self.gcmc_supercell = supercell
        info("%s supercell requested in config" % str(config_supercell))
        info("%s minimum supercell for a %.1f cutoff" %
             (str(minimum_supercell), config_cutoff))
        info("Constructing %s supercell for gcmc." % str(supercell))

    def supercell(self, scale):
        """Iterate over all the atoms of supercell."""
        if isinstance(scale, (int, long)):
            scale = (scale, scale, scale)
        for x_super in range(scale[0]):
            for y_super in range(scale[1]):
                for z_super in range(scale[2]):
                    offset = dot((x_super, y_super, z_super), self.cell.cell)
                    for atom in self.atoms:
                        newatom = copy(atom)
                        newatom.translate(offset)
                        yield newatom

    def order_by_types(self):
        """Sort the atoms alphabetically and group them."""
        self.atoms.sort(key=lambda x: (x.type, x.site))

    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def atominc_numbers(self):
        """Ordered list of atomic numbers."""
        return [atom.atomic_number for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def volume(self):
        """Unit cell volume."""
        return self.cell.volume

    @property
    def natoms(self):
        """Number of atoms in the unit cell."""
        return len(self.atoms)

    def get_gcmc_supercell(self):
        """Supercell used for gcmc."""
        return self.properties.get('supercell', (1, 1, 1))

    def set_gcmc_supercell(self, value):
        """Set the supercell property for the structure."""
        self.properties['supercell'] = value

    gcmc_supercell = property(get_gcmc_supercell, set_gcmc_supercell)
    # TODO(tdaff): properties: density, surface area, dft_energy, absorbance


class Cell(object):
    """
    Crystollagraphic cell representations and interconversion methods.

    Setter methods can be defined for different file types, however
    .cell and .params will be self-consistent if set directly.

    """

    def __init__(self):
        """Default to a 1A cubic box."""
        self._cell = array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])
        self._params = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    def from_pdb(self, line):
        """Extract cell from CRYST1 line in a pdb."""
        # Must use fixed widths as -ve numbers do not leave gaps to .split()
        self.params = (float(line[6:15]),
                       float(line[15:24]),
                       float(line[24:33]),
                       float(line[33:40]),
                       float(line[40:47]),
                       float(line[47:54]))

    def from_lines(self, lines, scale=1.0):
        """Extract cell from a 3-line POSCAR cell representation."""
        self.cell = array([[float(x) * scale for x in lines[0].split()],
                           [float(x) * scale for x in lines[1].split()],
                           [float(x) * scale for x in lines[2].split()]])

    def to_vector_strings(self, scale=1, bohr=False, fmt="%20.12f"):
        """Generic [Super]cell vectors in Angstrom as a list of strings."""
        out_format = 3 * fmt + "\n"
        if isinstance(scale, (int, long)):
            scale = [scale, scale, scale]
            # else assume an iterable 3-vector
        if bohr:
            scale = [x / BOHR2ANG for x in scale]
        return [out_format % tuple(scale[0] * self.cell[0]),
                out_format % tuple(scale[1] * self.cell[1]),
                out_format % tuple(scale[2] * self.cell[2])]

    def minimum_supercell(self, cutoff):
        """Calculate the smallest supercell with a half-cell width cutoff."""
        a_cross_b = cross(self.cell[0], self.cell[1])
        b_cross_c = cross(self.cell[1], self.cell[2])
        c_cross_a = cross(self.cell[2], self.cell[0])

        volume = dot(self.cell[0], b_cross_c)

        widths = [volume / norm(b_cross_c),
                  volume / norm(c_cross_a),
                  volume / norm(a_cross_b)]

        return tuple(int(ceil(2*cutoff/x)) for x in widths)

    @property
    def imcon(self):
        """Guess cell shape and return DL_POLY imcon key."""
        keys = {'none': 0,
                'cubic': 1,
                'orthorhombic': 2,
                'parallelepiped': 3}
        if np.all(self.cell == 0):
            return keys['none']
        elif np.all(self.params[3:] == 90):
            if self.params[0] == self.params[1] == self.params[2]:
                return keys['cubic']
            else:
                return keys['orthorhombic']
        else:
            return keys['parallelepiped']

    @property
    def volume(self):
        """Calculate cell volume a.bxc."""
        b_cross_c = cross(self.cell[1], self.cell[2])
        return dot(self.cell[0], b_cross_c)

    def get_cell(self):
        """Get the 3x3 vector cell representation."""
        return self._cell

    def set_cell(self, value):
        """Set cell and params from the cell representation."""
        self._cell = value
        self.__mkparam()

    # Property so that params are updated when cell is set
    cell = property(get_cell, set_cell)

    def get_params(self):
        """Get the six parameter cell representation."""
        return self._params

    def set_params(self, value):
        """Set cell and params from the cell parameters."""
        self._params = value
        self.__mkcell()

    # Property so that cell is updated when params are set
    params = property(get_params, set_params)

    # Implementation details -- directly access the private _{cell|param}
    # attributes; please don't break.
    def __mkcell(self):
        """Update the cell representation to match the parameters."""
        a_mag, b_mag, c_mag = self.params[:3]
        alpha, beta, gamma = [x * DEG2RAD for x in self.params[3:]]
        a_vec = array([a_mag, 0.0, 0.0])
        b_vec = array([b_mag * cos(gamma), b_mag * sin(gamma), 0.0])
        c_x = c_mag * cos(beta)
        c_y = c_mag * (cos(alpha) - cos(gamma) * cos(beta)) / sin(gamma)
        c_vec = array([c_x, c_y, (c_mag**2 - c_x**2 - c_y**2)**0.5])
        self._cell = array([a_vec, b_vec, c_vec])

    def __mkparam(self):
        """Update the parameters to match the cell."""
        cell_a = sqrt(sum(x**2 for x in self.cell[0]))
        cell_b = sqrt(sum(x**2 for x in self.cell[1]))
        cell_c = sqrt(sum(x**2 for x in self.cell[2]))
        alpha = arccos(sum(self.cell[1, :] * self.cell[2, :]) /
                       (cell_b * cell_c)) * 180 / pi
        beta = arccos(sum(self.cell[0, :] * self.cell[2, :]) /
                      (cell_a * cell_c)) * 180 / pi
        gamma = arccos(sum(self.cell[0, :] * self.cell[1, :]) /
                       (cell_a * cell_b)) * 180 / pi
        self._params = (cell_a, cell_b, cell_c, alpha, beta, gamma)


class Atom(object):
    """Base atom object."""

    def __init__(self, at_type=False, pos=False, **kwargs):
        """Accept arbritary kwargs as attributes."""
        self.type = at_type
        self.pos = pos
        self.charge = 0.0
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None
        # Sets anything else specified as an attribute
        for key, val in kwargs.iteritems():
            setattr(self, key, val)

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def fpos(self, cell):
        """Fractional position within a given cell."""
        return np.linalg.solve(cell.T, self.pos)

    def ifpos(self, cell):
        """In cell fractional position."""
        return [i % 1 for i in self.fpos(cell)]

    def ipos(self, cell):
        """In cell cartesian position."""
        return np.dot(self.ifpos(cell), cell)

    def from_cif(self, at_dict, cell, symmetry=None):
        """Extract an atom description from dictionary of cif items."""
        self.site = at_dict['_atom_site_label']
        # type_symbol takes precedence but need not be specified
        self.type = at_dict.get('_atom_site_type_symbol', self.site)
        self.mass = WEIGHT[self.type]
        frac_pos = [ufloat(at_dict['_atom_site_fract_x']),
                    ufloat(at_dict['_atom_site_fract_y']),
                    ufloat(at_dict['_atom_site_fract_z'])]
        if symmetry is not None:
            frac_pos = symmetry.trans_frac(frac_pos)
        self.pos = dot(frac_pos, cell)

    def from_pdb(self, line):
        """Parse the ATOM line from a pdb file."""
        # pdb is defined with fixed width fields rather than splitting
        self.idx = try_int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = try_int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()
        self.mass = WEIGHT[self.type]

    def from_vasp(self, line, at_type=None, cell=identity(3)):
        """Set the atom data from vasp input. Only pass cell if fractional."""
        self.pos = dot([float(x) for x in line.split()[:3]], cell)
        if at_type is not None:
            self.type = at_type
            self.mass = WEIGHT[at_type]

    def from_siesta(self, line, cell):
        """Parse line from SIESTA.STRUCT_OUT file."""
        self.pos = dot([float(x) for x in line.split()[2:5]], cell)
        self.atomic_number = int(line.split()[1])
        self.mass = WEIGHT[self.type]

    def from_xyz(self, line):
        """Parse line from generic xyz file."""
        self.pos = [float(x) for x in line.split()[1:4]]
        self.type = line.split()[0]
        self.mass = WEIGHT[self.type]

    def translate(self, vec):
        """Move the atom by the given vector."""
        self.pos = [x + y for x, y in zip(self.pos, vec)]

    def get_atomic_number(self):
        """The atomic number for the element, or closest match."""
        name = self.type
        atomic_number = None
        while name:
            try:
                atomic_number = ATOMIC_NUMBER.index(name)
                break
            except ValueError:
                name = name[-1]
        return atomic_number

    def set_atomic_number(self, value):
        """Set the atom type based on the atomic number."""
        self.type = ATOMIC_NUMBER[value]

    atomic_number = property(get_atomic_number, set_atomic_number)


class Guest(object):
    """Guest molecule and properties."""
    def __init__(self, ident=None, path=None):
        """Populate an empty guest then load from library if required."""
        self.ident = ''
        self.name = "Unknown guest"
        self.potentials = {}
        self.probability = []
        self.atoms = []
        self.source = "Unknown source"
        self.uptake = {}
        self.hoa = {}
        # only load if asked, set the ident in the loader
        if ident:
            self.load_guest(ident, path=None)

    def load_guest(self, ident, path=None):
        """Look in guests.lib in submit directory and default."""
        # Ident set here to keep consistent
        self.ident = ident
        # Need the different directories
        if path is None:
            job_dir = os.getcwd()
        else:
            job_dir = path
        if __name__ != '__main__':
            script_dir = os.path.dirname(__file__)
        else:
            script_dir = os.path.abspath(sys.path[0])
        job_guests = ConfigParser.SafeConfigParser()
        lib_guests = ConfigParser.SafeConfigParser()
        # Try and find guest in guests.lib
        job_guests.read(os.path.join(job_dir, 'guests.lib'))
        lib_guests.read(os.path.join(script_dir, 'guests.lib'))
        if job_guests.has_section(ident):
            debug("%s found in job dir" % ident)
            self._parse_guest(job_guests.items(ident))
        elif lib_guests.has_section(ident):
            debug("%s found in library" % ident)
            self._parse_guest(lib_guests.items(ident))
        else:
            error("Guest not found: %s" % ident)

    def _parse_guest(self, raw_text):
        """Set attributes according to the raw input."""
        for key, val in raw_text:
            if key == 'atoms':
                # Build a local list to replace what is there
                new_atoms = []
                # Only use non blank lines
                atoms = strip_blanks(val.splitlines())
                for atom in atoms:
                    atom = atom.split()
                    new_atoms.append(Atom(
                        at_type=atom[0],
                        mass=float(atom[1]),
                        charge=float(atom[2]),
                        pos=tuple(float(x) for x in atom[3:6])))
                self.atoms = new_atoms
            elif key == 'potentials':
                potens = strip_blanks(val.splitlines())
                for poten in potens:
                    poten = poten.split()
                    self.potentials[poten[0]] = tuple(float(x)
                                                      for x in poten[1:])
            elif key == 'probability':
                prob = strip_blanks(val.splitlines())
                self.probability = [tuple(int(y) for y in x.split())
                                    for x in prob]
            elif key == 'molar volume':
                self._molar_volume = val
            else:
                # Arbitrary attributes can be set
                # Might also enable breakage
                setattr(self, key, val)

    # Simple attribute-like calculations
    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def weight(self):
        """Unit cell weight."""
        return sum([atom.mass for atom in self.atoms])

    @property
    def molar_volume(self):
        """Molar volume at STP in dm3/mol"""
        if hasattr(self, '_molar_volume'):
            return self._molar_volume
        else:
            return 22.414


class Symmetry(object):
    """Apply symmetry operations to atomic coordinates."""
    def __init__(self, blob):
        """Read the operation from the argument."""
        self.sym_ops = []
        self.blob = blob
        self.parse_blob()
        self.cell = None

    def parse_blob(self):
        """Interpret a symmetry line from a cif."""
        # convert integers to floats to avoid integer division
        self.sym_ops = [re.sub(r'([\d]+)', r'\1.0', x.strip())
                        for x in re.split(',', self.blob) if x.strip()]

    def trans_frac(self, pos):
        """Apply symmetry operation to the supplied position."""
        new_pos = [eval(sym_op.replace('x', str(pos[0]))
                        .replace('y', str(pos[1]))
                        .replace('z', str(pos[2]))) for sym_op in self.sym_ops]
        return new_pos


def mk_repeat(cube_name='REPEAT_ESP.cube', symmetry=False):
    """Standard REPEAT input file."""
    if symmetry:
        symmetry_flag = 1
        # TODO(tdaff): connectivity.ff
    else:
        symmetry_flag = 0
    repeat_input = [
        "Input ESP file name in cube format\n",
        "%s\n" % cube_name,
        "Fit molecular(0) or periodic(1:default) system?\n",
        "1\n",
        "van der Waals scaling factor (default = 1.0)\n",
        "1.00000\n",
        "Apply RESP penalties?, no(0:default), yes(1)\n",
        "0\n",
        "Read cutoff radius? no(0), yes(1:default)\n",
        "1\n",
        "If flag above=1 provide R_cutoff next (in Bohrs)\n",
        "20.00000\n",
        "Apply symmetry restrain? no(0:default), yes(1)\n",
        "%i\n" % symmetry_flag,
        "Use Goddard-type restrains? no(0:default), yes(1)\n",
        "0\n",
        "If flag above=1 then provide weight next\n",
        "0.00000\n",
        "Enter total charge of the system\n",
        "0.00000\n"]

    filetemp = open('REPEAT_param.inp', 'w')
    filetemp.writelines(repeat_input)
    filetemp.close()


def mk_incar(options, esp_grid=None):
    """Basic vasp INCAR; use defaults as much as possible."""
    # We need these options
    job_name = options.get('job_name')
    spin = options.getbool('spin')
    optim_h = options.getbool('optim_h')
    optim_all = options.getbool('optim_all')
    optim_cell = options.getbool('optim_cell')
    dispersion = options.getbool('dispersion')

    incar = ["SYSTEM  = %s\n" % job_name,
             "ALGO    = Fast\n",
             "EDIFF   = 1E-5\n",
             "EDIFFG  = -0.02\n",
             "POTIM   = 0.4\n",
             "LREAL   = Auto\n",
             "LVTOT   = .TRUE.\n",
             "LVHAR   = .TRUE.\n",
             "ISMEAR  = 0\n",
             "SIGMA   = 0.05\n",
             "NWRITE  = 0\n"]
    if optim_cell:
        # Positions will be fixed by selective dynamics
        info("Cell vectors will be optimized")
        incar.extend(["ENCUT = 520\n",
                      "IBRION  = 2\n",
                      "NSW     = 800\n",
                      "ISIF    = 3\n"])
    elif optim_all or optim_h:
        # Just move positions
        incar.extend(["IBRION  = 2\n",
                      "NSW     = 300\n",
                      "ISIF    = 2\n"])
    else:
        # Single point energy
        info("Single point calculation")
        incar.extend(["IBRION  = 0\n",
                      "NSW     = 0\n",
                      "ISIF    = 0\n"])
    if spin:
        info("Spin polarised calculation")
        incar.append("ISPIN   = 2\n")

    if dispersion:
        info("Dispersion correction will be used")
        incar.append("LVDW    = .TRUE.\n")

    if esp_grid is not None:
        info("Changing FFT grid to %ix%ix%i" % esp_grid)
        incar.append("NGXF = %i ; NGYF = %i ; NGZF = %i\n" % esp_grid)

    return incar


def mk_kpoints(kpoints):
    """Defaults to gamma point only, or specified number."""
    if len(kpoints) != 3:
        error("kpoints specified incorectly; should be (i, i, i)")
    kpoints = [
        "Auto\n",
        "0\n",
        "Gamma\n",
        "%i %i %i\n" % tuple(kpoints),
        "0 0 0\n"]
    return kpoints


def mk_gcmc_control(temperature, pressures, options, guests):
    """Standard GCMC CONTROL file."""
    control = [
        "GCMC Run\n"
        "temperature  %f\n" % temperature,
        "steps    %i\n" % options.getint('mc_prod_steps'),
        "equilibration    %i\n" % options.getint('mc_eq_steps'),
        "cutoff          %f angstrom\n" % options.getfloat('mc_cutoff'),
        "delr            1.0 angstrom\n",
        "ewald precision  1d-6\n",
        "numguests %i\n" % options.getint('mc_numguests_freq'),
        "history %i\n" % options.getint('mc_history_freq')]
    if options.getbool('mc_jobcontrol'):
        control.append("jobcontrol\n")
    else:
        control.append("# jobcontrol\n")
    # Guest stuff
    if options.getbool('mc_probability_plot'):
        # We just comment out the probabilty plot lines if unneeded
        prob_on = ""
    else:
        prob_on = "# "
    if len(guests) > 1:
        gp_zip = zip(guests, pressures)
    else:
        gp_zip = zip(guests, [pressures])
    guest_count = 0
    for guest, press in gp_zip:
        guest_count += 1
        control.append("&guest %i\n" % guest_count)
        control.append("  pressure (bar) %f\n" % press)
        control.append("%s  probability %i\n" %
                       (prob_on, len(guest.probability)))
        for prob in guest.probability:
            control.append("%s  %i  " % (prob_on, len(prob)) +
                           "  ".join(["%i" % x for x in prob]) + "\n")
        control.append("&end\n")

    control.append("finish\n")
    return control


def unique(in_list, key=None):
    """Unique values in list ordered by first occurance"""
    uniq = []
    if key is not None:
        keys = []
        for item in in_list:
            item_key = key(item)
            if item_key not in keys:
                uniq.append(item)
                keys.append(item_key)
    else:
        for item in in_list:
            if item not in uniq:
                uniq.append(item)
    return uniq


def lorentz_berthelot(left, right):
    """Lorentz-Berthelot mixing rules for (sigma, epsilon) tuples."""
    sigma = (left[0] + right[0]) / 2.0
    epsilon = (left[1] * right[1])**0.5
    if epsilon == 0:
        sigma = 0
    return sigma, epsilon


# General utility functions
def mkdirs(directory):
    """Create a directory if it does not exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def terminate(exit_code=0):
    """Exit and announce if faps is terminating normally (default)."""
    if exit_code == 0:
        info("Faps terminated normally")
        raise SystemExit
    else:
        warning("Abnormal termination of faps; exit code %i" % exit_code)
        raise SystemExit(exit_code)


def move_and_overwrite(src, dest):
    """Move src to dest and overwrite if it is an existing file."""
    # As src and dest can be files or directories, do some checks.
    if os.path.exists(dest):
        if os.path.isdir(dest):
            dest_full = os.path.join(dest, os.path.basename(src))
            if os.path.exists(dest_full):
                if os.path.isfile(dest_full):
                    os.remove(dest_full)
                    shutil.move(src, dest)
                else:
                    raise OSError("Directory %s already exists" % dest_full)
            else:
                shutil.move(src, dest)
        elif os.path.isfile(dest):
            os.remove(dest)
            shutil.move(src, dest)
        else:
            raise OSError("%s is not a folder or file" % dest)
    else:
        shutil.move(src, dest)


def try_symlink(src, dest):
    """Delete an existing dest file, symlink a new one if possible."""
    if os.path.lexists(dest):
        os.remove(dest)
    try:
        os.symlink(src, dest)
    except AttributeError:
        shutil.copy(src, dest)


def same_guests(base, other):
    """Test if the guests are the same index and order."""
    return [item.ident for item in base] == [item.ident for item in other]


def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub('\(.*\)', '', text))


def try_int(text, default=0):
    """Try to parse an integer but return a default if it fails."""
    try:
        return int(text)
    except ValueError:
        return default


def prod(seq):
    """Calculate the product of all members of a sequence."""
    # numpy.prod will silently overflow 32 bit integer values
    # so we can use python bignums natively
    product = 1
    for item in seq:
        product *= item
    return product


def frac_near(pos_a, pos_b, epsilon=0.0002):
    """Return true if fractional points are close."""
    if epsilon < abs(pos_a[0] - pos_b[0]) < (1 - epsilon):
        return False
    elif epsilon < abs(pos_a[1] - pos_b[1]) < (1 - epsilon):
        return False
    elif epsilon < abs(pos_a[2] - pos_b[2]) < (1 - epsilon):
        return False
    else:
        return True


def remove_files(files, directory='.'):
    """Delete any of the files if they exist, or ignore if not found."""
    del_list = []
    for file_name in files:
        del_list.extend(glob.glob(os.path.join(directory, file_name)))
    for del_name in del_list:
        debug("deleting %s" % del_name)
        os.remove(del_name)


def compress_files(files, directory='.'):
    """Gzip any big files to keep."""
    zip_list = []
    for file_name in files:
        zip_list.extend(glob.glob(os.path.join(directory, file_name)))
    for zip_name in zip_list:
        debug("compressing %s" % zip_name)
        gzip_command = ['gzip', '-f', zip_name]
        subprocess.call(gzip_command)


def strip_blanks(lines):
    """Strip lines and remove blank lines."""
    return [line.strip() for line in lines if line.strip() != '']


def subgroup(iterable, width, itype=None):
    """Split an iterable into nested sub-itypes of width members."""
    # Return the same type as iterable
    if itype is None:
        if isinstance(iterable, list):
            itype = list
        else:
            itype = tuple
    # Will leave short groups if not enough members
    return itype([itype(iterable[x:x+width])
                  for x in range(0, len(iterable), width)])


def state_points(temperatures, pressures, individual, nguests):
    """Group temperatures and pressures and append given state point tuples."""
    # No error checking done, assume know what they are doing
    for index in range(0, len(individual), nguests+1):
        yield (individual[index], tuple(individual[index+1:index+nguests+1]))
    for temp in temperatures:
        for pressure in subgroup(pressures, nguests):
            yield (temp, pressure)


def welcome():
    """Print any important messages."""
    print(LOGO)
    print(("faps 0.999-r%s" % __version__.strip('$Revision: ')).rjust(79))

def main():
    """Do a standalone calculation when run as a script."""
    main_options = Options()
    info("Starting faps version 0.999-r%s" % __version__.strip('$Revision: '))
    # try to unpickle the job or
    # fall back to starting a new simulation
    niss_name = main_options.get('job_name') + ".niss"
    if os.path.exists(niss_name):
        info("Existing simulation found: %s; loading..." % niss_name)
        load_niss = open(niss_name)
        my_simulation = pickle.load(load_niss)
        load_niss.close()
        my_simulation.re_init(main_options)
    else:
        if not main_options.getbool('quiet'):
            welcome()
        info("Starting a new simulation...")
        my_simulation = PyNiss(main_options)

    # run requested jobs
    my_simulation.job_dispatcher()
    my_simulation.dump_state()
    terminate(0)


if __name__ == '__main__':
    main()
