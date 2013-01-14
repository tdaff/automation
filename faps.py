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
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import glob
import os
import pickle
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import textwrap
import time
from copy import copy
from itertools import count
from logging import warning, debug, error, info, critical
from math import ceil, log
from os import path

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos
from numpy import array, identity, dot, cross
from numpy.linalg import norm

from config import Options
from elements import WEIGHT, ATOMIC_NUMBER, UFF, VASP_PSEUDO_PREF
from elements import CCDC_BOND_ORDERS, GULP_BOND_ORDERS, METALS
from elements import COVALENT_RADII, UFF_FULL
from job_handler import JobHandler
from logo import LOGO

# Global constants
DEG2RAD = pi / 180.0
BOHR2ANG = 0.52917720859
EV2KCAL = 23.060542301389
NAVOGADRO = 6.02214129E23
INFINITY = float('inf')
KCAL_TO_KJ = 4.1868  # Steam tables from openbabel

FASTMC_DEFAULT_GRID_SPACING = 0.1

# ID values for system state
NOT_RUN = 0
RUNNING = 1
FINISHED = 2
UPDATED = 3
SKIPPED = -1
NOT_SUBMITTED = -2

# Possible folder names; need these so that similar_ones_with_underscores are
# not globbed
FOLDER_SUFFIXES = ['gulp', 'gulp_opt', 'gulp_fit', 'siesta', 'vasp', 'egulp',
                   'repeat', 'fastmc', 'properties']


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
                      'ff_opt': (NOT_RUN, False),
                      'dft': (NOT_RUN, False),
                      'esp': (NOT_RUN, False),
                      'charges': (NOT_RUN, False),
                      'properties': (NOT_RUN, False),
                      'gcmc': {}}
        self.job_handler = JobHandler(options)

    def dump_state(self):
        """Write the .niss file holding the current system state."""
        job_name = self.options.get('job_name')
        info("Writing state file, %s.niss." % job_name)
        os.chdir(self.options.get('job_dir'))
        # Don't save the job handler in case it changes
        save_handler = self.job_handler
        self.job_handler = None
        my_niss = open(job_name + ".niss", "wb")
        pickle.dump(self, my_niss)
        my_niss.close()
        # put the job handler back and continue
        self.job_handler = save_handler

    def re_init(self, new_options):
        """Re initialize simulation (with updated options)."""
        if new_options.getbool('update_opts'):
            info("Using new options.")
            self.options = new_options
            self.structure.name = new_options.get('job_name')
        else:
            # Just update command line stuff
            info("Using old options with new command line arguments.")
            self.options.args = new_options.args
            self.options.options = new_options.options
            self.options.cmdopts = new_options.cmdopts
            self.structure.name = new_options.get('job_name')
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
            code_locals = locals()
            code_locals.update(globals())
            console = code.InteractiveConsole(code_locals)
            console.push('import rlcompleter, readline')
            console.push('readline.parse_and_bind("tab: complete")')
            banner = ("((-----------------------------------------------))\n"
                      "((             Interactive faps mode             ))\n"
                      "((             =====================             ))\n"
                      "((                                               ))\n"
                      "(( WARNING: mode is designed for devs and        ))\n"
                      "(( experts only!                                 ))\n"
                      "(( Current simulation is accessed as 'self' and  ))\n"
                      "(( associated methods. Type 'dir()' to see the   ))\n"
                      "(( methods in the local namespace and 'help(x)'  ))\n"
                      "(( for help on any object.                       ))\n"
                      "(( Use 'self.dump_state()' to save any changes.  ))\n"
                      "((-----------------------------------------------))")
            console.interact(banner=banner)

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

        self.step_force_field()

        self.step_dft()

        self.step_charges()

        if self.options.getbool('qeq_fit'):
            if not 'qeq_fit' in self.state and self.state['charges'][0] == UPDATED:
                info("QEq parameter fit requested")
                self.run_qeq_gulp(fitting=True)
                self.dump_state()

        self.step_gcmc()

        self.step_properties()

        self.send_to_database()

        self.post_summary()


    def status(self, initial=False):
        """Print the current status to the terminal."""
        valid_states = {NOT_RUN: 'Not run',
                        RUNNING: 'Running',
                        FINISHED: 'Finished',
                        UPDATED: 'Processed',
                        SKIPPED: 'Skipped',
                        NOT_SUBMITTED: 'Not submitted'}

        if initial:
            info("Previous system state reported from .niss file "
                 "(running jobs may have already finished):")
        else:
            info("Current system status:")
        for step, state in self.state.items():
            if step == 'gcmc':
                if not state:
                    info(" * State of GCMC: Not run")
                else:
                    for point, job in state.items():
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

    def send_to_database(self):
        """If using a database, store the results"""

        # we can skip if not using a database
        if not 'sql' in self.options.get('initial_structure_format'):
            return

        # extract the database and structure names
        db_params = self.options.get('job_name').split('.')

        # import this here so sqlalchemy is not required generally
        from backend.sql import AlchemyBackend

        database = AlchemyBackend(db_params[0])
        info("Storing results in database")
        database.store_results(db_params[1], int(db_params[2]), self.structure)
        debug("Database finished")

    def post_summary(self):
        """Summarise any results for GCMC, properties..."""
        # Also includes excess calculation if void volume calculated
        # N = pV/RT
        R_GAS = 8.3144621E25 / NAVOGADRO # A^3 bar K-1 molecule
        job_name = self.options.get('job_name')
        info("Summary of GCMC results")
        info("======= ======= ======= ======= =======")
        nguests = len(self.structure.guests)
        for idx, guest in enumerate(self.structure.guests):
            # Determine whether we can calulate the excess for
            # any different probes
            void_volume = self.structure.sub_property('void_volume')
            he_excess, guest_excess = "", ""
            if 1.0 in void_volume:
                he_excess = 'He-xs-molc/uc,He-xs-mmol/g,He-xs-v/v,He-xs-wt%,'
            if hasattr(guest, 'probe_radius'):
                if guest.probe_radius != 1.0 and guest.probe_radius in void_volume:
                    guest_excess = 'xs-molc/uc,xs-mmol/g,xs-v/v,xs-wt%,'
            # Generate headers separately
            csv = ["#T/K,p/bar,molc/uc,mmol/g,stdev,",
                   "v/v,stdev,wt%,stdev,hoa/kcal/mol,stdev,",
                   guest_excess, he_excess,
                   ",".join("p(g%i)" % gidx for gidx in range(nguests)), "\n"]
            info(guest.name)
            info("---------------------------------------")
            info("molc/uc  mmol/g  vstp/v     hoa     T_P")
            info("======= ======= ======= ======= =======")
            for tp_point in sorted(guest.uptake):
                # <N>, sd, supercell
                uptake = guest.uptake[tp_point]
                uptake = [uptake[0]/uptake[2], uptake[1]/uptake[2]]
                hoa = guest.hoa[tp_point]
                # uptake in mmol/g
                muptake = 1000*uptake[0]/self.structure.weight
                muptake_stdev = 1000*uptake[1]/self.structure.weight
                # volumetric uptake
                vuptake = (guest.molar_volume*uptake[0]/
                           (6.023E-4*self.structure.volume))
                vuptake_stdev = (guest.molar_volume*uptake[1]/
                                 (6.023E-4*self.structure.volume))
                # weight percent uptake
                wtpc = 100*(1 - self.structure.weight/
                            (self.structure.weight + uptake[0]*guest.weight))
                wtpc_stdev = 100*(1 - self.structure.weight/
                                  (self.structure.weight + uptake[1]*guest.weight))
                info("%7.2f %7.2f %7.2f %7.2f %s" % (
                    uptake[0], muptake, vuptake, hoa[0],
                    ("T=%s" % tp_point[0] +
                     ''.join(['P=%s' % x for x in tp_point[1]]))))
                csv.append("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f," % (
                    tp_point[0], tp_point[1][idx], uptake[0],
                    muptake, muptake_stdev,
                    vuptake, vuptake_stdev,
                    wtpc, wtpc_stdev,
                    hoa[0], hoa[1]))
                if guest_excess:
                    guest_void = void_volume[guest.probe_radius]
                    n_bulk = (tp_point[1][idx]*guest_void)/(tp_point[0]*R_GAS)
                    xs_uptake = uptake[0]-n_bulk
                    # uptake in mmol/g
                    muptake = 1000*xs_uptake/self.structure.weight
                    # volumetric uptake
                    vuptake = (guest.molar_volume*xs_uptake/
                               (6.023E-4*self.structure.volume))
                    # weight percent uptake
                    wtpc = 100*(1 - self.structure.weight/
                                (self.structure.weight + xs_uptake*guest.weight))
                    csv.append("%f,%f,%f,%f," % (
                        xs_uptake, muptake, vuptake, wtpc,))
                if he_excess:
                    guest_void = void_volume[1.0]
                    n_bulk = (tp_point[1][idx]*guest_void)/(tp_point[0]*R_GAS)
                    xs_uptake = uptake[0]-n_bulk
                    # uptake in mmol/g
                    muptake = 1000*xs_uptake/self.structure.weight
                    # volumetric uptake
                    vuptake = (guest.molar_volume*xs_uptake/
                               (6.023E-4*self.structure.volume))
                    # weight percent uptake
                    wtpc = 100*(1 - self.structure.weight/
                                (self.structure.weight + xs_uptake*guest.weight))
                    csv.append("%f,%f,%f,%f," % (
                        xs_uptake, muptake, vuptake, wtpc,))
                # list all the other guest pressures and start a new line
                csv.append(",".join("%f" % x for x in tp_point[1]) + "\n")

            csv_file = open('%s-%s.csv' % (job_name, guest.ident), 'w')
            csv_file.writelines(csv)
            csv_file.close()
        info("======= ======= ======= ======= =======")

        info("Structure properties")

        # Internally calculated surface area
        surf_area_results = self.structure.surface_area()
        if surf_area_results:
            info("Summary of faps surface areas")
            info("========= ========= ========= =========")
            info(" radius/A total/A^2  m^2/cm^3     m^2/g")
            info("========= ========= ========= =========")
            for probe, area in surf_area_results.items():
                vol_area = 1E4*area/self.structure.volume
                specific_area = NAVOGADRO*area/(1E20*self.structure.weight)
                info("%9.3f %9.2f %9.2f %9.2f" %
                     (probe, area, vol_area, specific_area))
            info("========= ========= ========= =========")

        # Messy, but check individual properties that might not be there
        # and dump them to the screen
        info("weight (u): %f" % self.structure.weight)
        if hasattr(self.structure, 'pore_diameter'):
            info("pores (A): %f %f %f" % self.structure.pore_diameter)
        channel_results = self.structure.sub_property('dimensionality')
        if channel_results:
            for probe, channels in channel_results.items():
                info(("channels %.2f probe: " % probe) +
                     " ".join("%i" % x for x in channels))
        # The table is copied from above as it does some calculating
        surf_area_results = self.structure.sub_property('zeo_surface_area')
        if surf_area_results:
            info("Summary of zeo++ surface areas")
            info("========= ========= ========= =========")
            info(" radius/A total/A^2  m^2/cm^3     m^2/g")
            info("========= ========= ========= =========")
            for probe, area in surf_area_results.items():
                vol_area = 1E4*area/self.structure.volume
                specific_area = NAVOGADRO*area/(1E20*self.structure.weight)
                info("%9.3f %9.2f %9.2f %9.2f" %
                     (probe, area, vol_area, specific_area))
            info("========= ========= ========= =========")

        info("volume (A^3): %f" % self.structure.volume)

        void_volume_results = self.structure.sub_property('void_volume')
        if surf_area_results:
            info("Summary of zeo++ void volumes")
            info("========= ========= ========= =========")
            info(" radius/A total/A^3  fraction    cm^3/g")
            info("========= ========= ========= =========")
            for probe, void in void_volume_results.items():
                void_fraction = void/self.structure.volume
                specific_area = NAVOGADRO*void/(1E24*self.structure.weight)
                info("%9.3f %9.2f %9.5f %9.4f" %
                     (probe, void, void_fraction, specific_area))
            info("========= ========= ========= =========")


    def step_force_field(self):
        """Check the force field step of the calculation."""
        end_after = False

        if 'ff_opt' not in self.state:
            self.state['ff_opt'] = (NOT_RUN, False)

        if self.state['ff_opt'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_force_field_opt'):
                info("Skipping force field optimisation")
                self.state['ff_opt'] = (SKIPPED, False)
            elif self.state['ff_opt'][0] == RUNNING:
                job_state = self.job_handler.jobcheck(self.state['ff_opt'][1])
                if not job_state:
                    info("Queue reports force field optimisation has finished")
                    self.state['ff_opt'] = (FINISHED, False)
                else:
                    # Still running
                    info("Force field optimisation still in progress")
                    end_after = True

        if self.state['ff_opt'][0] == NOT_RUN or 'ff_opt' in self.options.args:
            jobid = self.run_ff_opt()
            sys_argv_strip('ff_opt')
            end_after = self.postrun(jobid)
            self.dump_state()

        if self.state['ff_opt'][0] == FINISHED:
            self.structure.update_pos(self.options.get('ff_opt_code'))
            self.state['ff_opt'] = (UPDATED, False)
            self.dump_state()

        # If postrun is submitted then this script is done!
        if end_after:
            terminate(0)

    def step_dft(self):
        """Check the DFT step of the calculation."""
        end_after = False

        if self.state['dft'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_dft'):
                info("Skipping DFT step completely")
                info("Job might fail later if you need the ESP")
                self.state['dft'] = (SKIPPED, False)
            elif self.state['dft'][0] == RUNNING:
                job_state = self.job_handler.jobcheck(self.state['dft'][1])
                if not job_state:
                    info("Queue reports DFT step has finished")
                    self.state['dft'] = (FINISHED, False)
                else:
                    # Still running
                    info("DFT still in progress")
                    end_after = True

        if self.state['dft'][0] == NOT_RUN or 'dft' in self.options.args:
            jobid = self.run_dft()
            sys_argv_strip('dft')
            end_after = self.postrun(jobid)
            self.dump_state()

        if self.state['dft'][0] == FINISHED:
            self.structure.update_pos(self.options.get('dft_code'))
            self.state['dft'] = (UPDATED, False)
            self.dump_state()

        # If postrun is submitted then this script is done!
        if end_after:
            terminate(0)

    def step_charges(self):
        """Check the charge step of the calculation."""
        end_after = False

        if self.state['charges'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_charges'):
                info("Skipping charge calculation")
                self.state['charges'] = (SKIPPED, False)
            elif self.state['charges'][0] == RUNNING:
                job_state = self.job_handler.jobcheck(self.state['charges'][1])
                if not job_state:
                    info("Queue reports charge calculation has finished")
                    self.state['charges'] = (FINISHED, False)
                else:
                    info("Charge calculation still running")
                    end_after = True

        if self.state['charges'][0] == NOT_RUN or 'charges' in self.options.args:
            jobid = self.run_charges()
            sys_argv_strip('charges')
            end_after = self.postrun(jobid)
            self.dump_state()

        if self.state['charges'][0] == FINISHED:
            self.structure.update_charges(self.options.get('charge_method'),
                                          self.options)
            self.state['charges'] = (UPDATED, False)
            self.dump_state()

        # If postrun is submitted then this script is done!
        if end_after:
            terminate(0)

    def step_gcmc(self):
        """Check the GCMC step of the calculation."""
        end_after = False
        jobids = {}
        postrun_ids = []

        if self.options.getbool('no_gcmc'):
            info("Skipping GCMC simulation")
            return
        elif not self.state['gcmc'] or 'gcmc' in self.options.args:
            # The dictionary is empty before any runs
            info("Starting gcmc step")
            jobids = self.run_fastmc()
            sys_argv_strip('gcmc')
            self.dump_state()

        for tp_point, jobid in jobids.items():
            if jobid is True:
                self.state['gcmc'][tp_point] = (FINISHED, False)
            elif jobid is False:
                self.state['gcmc'][tp_point] = (SKIPPED, False)
            else:
                info("FastMC job in queue. Jobid: %s" % jobid)
                self.state['gcmc'][tp_point] = (RUNNING, jobid)
                postrun_ids.append(jobid)
                # unfinished GCMCs
                end_after = True
        else:
            # when the loop completes write out the state
            self.dump_state()

        for tp_point in self.state['gcmc']:
            tp_state = self.state['gcmc'][tp_point]
            if tp_state[0] == RUNNING:
                new_state = self.job_handler.jobcheck(tp_state[1])
                if not new_state:
                    info("Queue reports GCMC %s finished" % (tp_point,))
                    # need to know we have finished to update below
                    tp_state = (FINISHED, False)
                    self.state['gcmc'][tp_point] = tp_state
                    self.dump_state()
                else:
                    info("GCMC %s still running" % (tp_point,))
                    # unfinished GCMC so exit later
                    end_after = True

            # any states that need to be updated should have been done by now
            if tp_state[0] == FINISHED:
                startdir = os.getcwd()
                # wooki seems slow to copy output files back
                # so we give them a few chances to appear
                max_attempts = 6
                for attempt_count in range(max_attempts):
                    time.sleep(attempt_count)
                    try:
                        self.structure.update_gcmc(tp_point, self.options)
                        self.state['gcmc'][tp_point] = (UPDATED, False)
                        self.dump_state()
                        break
                    except IOError:
                        os.chdir(startdir)
                else:
                    error('OUTPUT file never appeared')

        if postrun_ids:
            self.postrun(postrun_ids)

        if end_after:
            info("GCMC run has not finished completely")
            terminate(0)

    def step_properties(self):
        """Run the properties calculations if required."""

        if self.state['properties'][0] not in [UPDATED, SKIPPED]:
            if self.options.getbool('no_properties'):
                info("Skipping all properties calculations")
                self.state['properties'] = (SKIPPED, False)

        if self.state['properties'][0] == NOT_RUN or 'properties' in self.options.args:
            self.calculate_properties()
            self.state['properties'] = (UPDATED, False)
            self.dump_state()

    def import_old(self):
        """Try and import any data from previous stopped simulation."""
        job_name = self.options.get('job_name')
        job_dir = self.options.get('job_dir')
        try:
            self.structure.from_file(
                job_name,
                self.options.get('initial_structure_format'),
                self.options)
            self.state['init'] = (UPDATED, False)
        except IOError:
            info("No initial structure found to import")
        try:
            self.structure.update_pos(self.options.get('ff_opt_code'))
            self.state['ff_opt'] = (UPDATED, False)
        except IOError:
            info("No force field optimised structure found to import")
        try:
            self.structure.update_pos(self.options.get('dft_code'))
            self.state['dft'] = (UPDATED, False)
        except IOError:
            info("No optimized structure found to import")
        try:
            self.structure.update_charges(self.options.get('charge_method'),
                                          self.options)
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
                self.structure.update_gcmc(tp_point, self.options)
                self.state['gcmc'][tp_point] = (UPDATED, False)
            except (IOError, OSError):
                info("GCMC point %s not found" % str(tp_point))

        # Reset directory at end
        os.chdir(job_dir)

    def postrun(self, jobid):
        """Determine if we need the job handler to post submit itself."""
        # update the job tracker
        if jobid is not False and jobid is not True:
            if self.options.getbool('run_all'):
                debug('Submitting postrun script')
                os.chdir(self.options.get('job_dir'))
                self.job_handler.postrun(jobid)
                return True
            else:
                debug('Postrun script not submitted')
                return False
        else:
            return False


    def run_ff_opt(self):
        """Select correct method for running the dft/optim."""
        ff_opt_code = self.options.get('ff_opt_code')
        info("Running a %s calculation" % ff_opt_code)
        if ff_opt_code == 'gulp':
            jobid = self.run_optimise_gulp()
        else:
            critical("Unknown force field method: %s" % ff_opt_code)
            terminate(91)

        if jobid is True:
            # job run and finished
            self.state['ff_opt'] = (FINISHED, False)
        else:
            info("Running %s job in queue. Jobid: %s" % (ff_opt_code, jobid))
            self.state['ff_opt'] = (RUNNING, jobid)

        return jobid


    def run_dft(self):
        """Select correct method for running the dft/optim."""
        dft_code = self.options.get('dft_code')
        info("Running a %s calculation" % dft_code)
        if dft_code == 'vasp':
            jobid = self.run_vasp()
        elif dft_code == 'siesta':
            jobid = self.run_siesta()
        else:
            critical("Unknown dft method: %s" % dft_code)
            terminate(92)

        # update the job tracker
        #if jobid is False:
        #    self.state['dft'] = (NOT_SUBMITTED, False)
            # submission skipped
        if jobid is True:
            # job run and finished
            self.state['dft'] = (FINISHED, False)
        else:
            info("Running %s job in queue. Jobid: %s" % (dft_code, jobid))
            self.state['dft'] = (RUNNING, jobid)

        return jobid


    def run_charges(self):
        """Select correct charge processing methods."""
        chg_method = self.options.get('charge_method')
        info("Calculating charges with %s" % chg_method)
        if chg_method == 'repeat':
            # Get ESP
            self.esp_to_cube()
            jobid = self.run_repeat()
        elif chg_method == 'gulp':
            jobid = self.run_qeq_gulp()
        elif chg_method == 'egulp':
            jobid = self.run_qeq_egulp()
        else:
            critical("Unknown charge calculation method: %s" % chg_method)
            terminate(93)

        # update the job tracker
        if jobid is True:
            # job run and finished
            self.state['charges'] = (FINISHED, False)
        else:
            info("Running %s job in queue. Jobid: %s" % (chg_method, jobid))
            self.state['charges'] = (RUNNING, jobid)

        return jobid

    ## Methods for specific codes start here

    def run_optimise_gulp(self):
        """Run GULP to do a UFF optimisation."""
        job_name = self.options.get('job_name')
        optim_code = 'gulp'
        terse = self.options.getbool('gulp_terse')
        # put an opt in path to distinguish from the charge calculation
        optim_dir = path.join(self.options.get('job_dir'),
                                   'faps_%s_%s_opt' % (job_name, optim_code))
        mkdirs(optim_dir)
        os.chdir(optim_dir)
        debug("Running in %s" % optim_dir)

        filetemp = open('%s.gin' % job_name, 'w')
        filetemp.writelines(self.structure.to_gulp(optimise=True, terse=terse))
        filetemp.close()

        if 'GULP_LIB' not in os.environ:
            warning("gulp library directory not set; optimisation might fail")

        if self.options.getbool('no_submit'):
            info("GULP input files generated; skipping job submission")
            jobid = False
        else:
            jobid = self.job_handler.submit(optim_code, self.options,
                                            input_file='%s.gin' % job_name)

        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))
        return jobid


    def run_vasp(self):
        """Make inputs and run vasp job."""
        job_name = self.options.get('job_name')
        nproc = self.options.getint('vasp_ncpu')
        # Keep things tidy in a subdirectory
        dft_code = self.options.get('dft_code')
        vasp_dir = path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, dft_code))
        mkdirs(vasp_dir)
        os.chdir(vasp_dir)
        debug("Running in %s" % vasp_dir)
        info("Running on %i nodes" % nproc)

        filetemp = open("POSCAR", "w")
        filetemp.writelines(self.structure.to_vasp(self.options))
        filetemp.close()

        esp_grid = self.esp_grid

        filetemp = open("INCAR", "w")
        if self.esp_reduced:
            # Let VASP do the grid if we don't need to
            filetemp.writelines(mk_incar(self.options, esp_grid=esp_grid))
        else:
            filetemp.writelines(mk_incar(self.options))
        filetemp.close()

        filetemp = open("KPOINTS", "w")
        filetemp.writelines(mk_kpoints(self.options.gettuple('kpoints', int)))
        filetemp.close()

        potcar_types = unique(self.structure.types)
        filetemp = open("POTCAR", "w")
        potcar_dir = self.options.get('potcar_dir')
        for at_type in potcar_types:
            # Try and get the preferred POTCARS
            debug("Using %s pseudopotential for %s" %
                 (VASP_PSEUDO_PREF.get(at_type, at_type), at_type))
            potcar_src = path.join(potcar_dir,
                                      VASP_PSEUDO_PREF.get(at_type, at_type),
                                      "POTCAR")
            shutil.copyfileobj(open(potcar_src), filetemp)
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("Vasp input files generated; skipping job submission")
            # act as if job completed
            jobid = False
        else:
            self.job_handler.env(dft_code, options=self.options)
            jobid = self.job_handler.submit(dft_code, self.options)

        # Tidy up at the end and pass on job id
        os.chdir(self.options.get('job_dir'))
        return jobid

    def run_siesta(self):
        """Make siesta input and run job."""
        job_name = self.options.get('job_name')
        nproc = self.options.getint('siesta_ncpu')
        # Keep things tidy in a subdirectory
        dft_code = self.options.get('dft_code')
        siesta_dir = path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, dft_code))
        mkdirs(siesta_dir)
        os.chdir(siesta_dir)
        debug("Running in %s" % siesta_dir)
        info("Running on %i nodes" % nproc)

        filetemp = open('%s.fdf' % job_name, 'w')
        filetemp.writelines(self.structure.to_siesta(self.options))
        filetemp.close()

        psf_types = unique(self.structure.types)
        psf_dir = self.options.get('psf_dir')
        for at_type in psf_types:
            psf_atm = '%s.psf' % at_type
            psf_src = path.join(psf_dir, psf_atm)
            psf_dest = path.join(siesta_dir, psf_atm)
            try:
                if not path.exists(psf_atm):
                    os.symlink(psf_src, psf_dest)
            # symlinks not available pre 3.2 on windows
            except AttributeError:
                shutil.copy(psf_src, siesta_dir)
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("Siesta input files generated; skipping job submission")
            jobid = False
        else:
            # sharcnet does weird things for siesta
            self.job_handler.env(dft_code, options=self.options)
            jobid = self.job_handler.submit(dft_code, self.options,
                                            input_file='%s.fdf' % job_name)

        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))
        return jobid

    def run_qeq_gulp(self, fitting=False):
        """Run GULP to calculate charge equilibration charges."""
        job_name = self.options.get('job_name')
        qeq_code = 'gulp'
        if fitting:
            qeq_dir = path.join(self.options.get('job_dir'),
                                   'faps_%s_%s_fit' % (job_name, qeq_code))
        else:
            qeq_dir = path.join(self.options.get('job_dir'),
                                   'faps_%s_%s' % (job_name, qeq_code))
        mkdirs(qeq_dir)
        os.chdir(qeq_dir)
        debug("Running in %s" % qeq_dir)

        qeq_dict = parse_qeq_params(self.options.gettuple('qeq_parameters'))

        filetemp = open('%s.gin' % job_name, 'w')
        filetemp.writelines(self.structure.to_gulp(qeq_fit=fitting, qeq_dict=qeq_dict))
        filetemp.close()

        if self.options.getbool('no_submit'):
            info("GULP input files generated; skipping job submission")
            jobid = False
        elif fitting:
            jobid = self.job_handler.submit(qeq_code, self.options,
                                            input_file='%s.gin' % job_name)
            info("Running GULP fitting job in queue. Jobid: %s" % jobid)
            self.state['qeq_fit'] = (RUNNING, jobid)
        else:
            jobid = self.job_handler.submit(qeq_code, self.options,
                                            input_file='%s.gin' % job_name)

        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))
        return jobid

    def run_qeq_egulp(self):
        """Run EGULP to calculate charge equilibration charges."""
        job_name = self.options.get('job_name')
        qeq_code = 'egulp'
        qeq_dir = path.join(self.options.get('job_dir'),
                            'faps_%s_%s' % (job_name, qeq_code))
        mkdirs(qeq_dir)
        os.chdir(qeq_dir)
        debug("Running in %s" % qeq_dir)

        filetemp = open('%s.geo' % job_name, 'w')
        filetemp.writelines(self.structure.to_egulp())
        filetemp.close()

        # EGULP defaults to GULP parameters if not specified
        # TODO(tdaff): remove old terminaology
        try:
            egulp_parameters = self.options.gettuple('egulp_parameters')
            warning("egulp_parameters is deprecated use qeq_parameters instead")
        except AttributeError:
            egulp_parameters = self.options.gettuple('qeq_parameters')

        if egulp_parameters:
            info("Custom EGULP parameters selected")
            filetemp = open('%s.param' % job_name, 'w')
            filetemp.writelines(mk_egulp_params(egulp_parameters))
            filetemp.close()
            # job handler needs to know about both these files
            egulp_args = ['%s.geo' % job_name, '%s.param' % job_name]
        else:
            # Only geometry file is needed for this run
            egulp_args = ['%s.geo' % job_name]

        if self.options.getbool('no_submit'):
            info("EGULP input files generated; skipping job submission")
            jobid = False
        else:
            jobid = self.job_handler.submit(qeq_code, self.options,
                                            input_args=egulp_args)

        # Tidy up at the end
        os.chdir(self.options.get('job_dir'))
        return jobid

    def esp_to_cube(self):
        """Make the cube for repeat input."""
        job_name = self.options.get('job_name')
        esp_src = self.options.get('esp_src')
        repeat_dir = path.join(self.options.get('job_dir'),
                                  'faps_%s_repeat' % job_name)
        mkdirs(repeat_dir)
        src_dir = path.join(self.options.get('job_dir'),
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
                error("No cube files found; check vasp_to_cube output")
            # Move it to the repeat directory and give a proper name
            move_and_overwrite(cube_file,
                               path.join(repeat_dir, job_name + '.cube'))
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
        repeat_dir = path.join(self.options.get('job_dir'),
                                  'faps_%s_%s' % (job_name, charge_code))
        mkdirs(repeat_dir)
        os.chdir(repeat_dir)

        if self.options.getbool('symmetry'):
            mk_repeat(cube_name=job_name + '.cube', symmetry=True)
            mk_connectivity_ff(self.structure.symmetry_tree)
        else:
            mk_repeat(cube_name=job_name + '.cube', symmetry=False)
        if self.options.getbool('no_submit'):
            info("REPEAT input files generated; skipping job submission")
            jobid = False
        else:
            jobid = self.job_handler.submit(charge_code, self.options)

        os.chdir(self.options.get('job_dir'))
        return jobid

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

        gcmc_dir = path.join(self.options.get('job_dir'),
                                'faps_%s_%s' % (job_name, mc_code))
        mkdirs(gcmc_dir)
        os.chdir(gcmc_dir)

        config, field = self.structure.to_fastmc(self.options)

        filetemp = open("CONFIG", "w")
        filetemp.writelines(config)
        filetemp.close()

        filetemp = open("FIELD", "w")
        filetemp.writelines(field)
        filetemp.close()

        temps = self.options.gettuple('mc_temperature', float)
        presses = self.options.gettuple('mc_pressure', float)
        indivs = self.options.gettuple('mc_state_points', float)
        jobids = {}
        for tp_point in state_points(temps, presses, indivs, len(guests)):
            temp = tp_point[0]
            press = tp_point[1]
            info("Running GCMC: T=%.1f " % temp +
                 " ".join(["P=%.2f" % x for x in press]))
            tp_path = ('T%s' % temp +
                       ''.join(['P%.2f' % x for x in press]))
            mkdirs(tp_path)
            os.chdir(tp_path)
            try_symlink(path.join('..', 'CONFIG'),'CONFIG')
            try_symlink(path.join('..', 'FIELD'), 'FIELD')
            filetemp = open("CONTROL", "w")
            filetemp.writelines(mk_gcmc_control(temp, press, self.options,
                                                guests, self.structure.gcmc_supercell))
            filetemp.close()

            if self.options.getbool('no_submit'):
                info("FastMC input files generated; "
                     "skipping job submission")
                jobids[(temp, press)] = False
            else:
                jobid = self.job_handler.submit(mc_code, self.options)
                jobids[(temp, press)] = jobid
            os.chdir('..')

        os.chdir(self.options.get('job_dir'))
        return jobids

    def calculate_properties(self):
        """Calculate general structural properties."""

        job_name = self.options.get('job_name')
        job_dir = self.options.get('job_dir')
        props_dir = path.join(job_dir, 'faps_%s_properties' % job_name)
        mkdirs(props_dir)
        os.chdir(props_dir)

        self.structure.gen_neighbour_list()

        ##
        # Surface area calculations
        ##
        surf_probes = self.options.gettuple('surface_area_probe', dtype=float)
        for probe in surf_probes:
            if self.structure.surface_area(probe) is None:
                self.structure.surface_area(probe, value=self.calc_surface_area(probe))

        # Neighbour list makes .niss too big; remove them
        for atom in self.structure.atoms:
            atom.neighbours = None
            del atom.neighbours

        # Zeoplusplus gives fast access to many properties
        if self.options.getbool('zeo++'):
            try:
                self.calculate_zeo_properties()
            except (OSError, IOError):
                error("Error running zeo++; skipping")

        os.chdir(job_dir)


    def calculate_zeo_properties(self):
        """Run the zeo++ and update properties with no error trapping."""

        job_name = self.options.get('job_name')
        zeofiles = self.structure.to_zeoplusplus()

        filetemp = open("%s.cssr" % job_name, 'w')
        filetemp.writelines(zeofiles[0])
        filetemp.close()

        filetemp = open("%s.rad" % job_name, 'w')
        filetemp.writelines(zeofiles[1])
        filetemp.close()

        filetemp = open("%s.mass" % job_name, 'w')
        filetemp.writelines(zeofiles[2])
        filetemp.close()

        probes = set([1.0]) # Always have a helium probe
        for guest in self.structure.guests:
            if hasattr(guest, 'probe_radius'):
                probes.add(guest.probe_radius)

        zeo_exe = shlex.split(self.options.get('zeo++_exe'))
        zeo_exe += ['-mass', '%s.mass' % job_name, '-r', '%s.rad' % job_name]
        cssr_file = ['%s.cssr' % job_name]

        # incuded sphere, free sphere, included sphere along free path
        zeo_command = zeo_exe + ['-res'] + cssr_file
        info("Running zeo++ pore diameters")
        debug("Running command: '" + " ".join(zeo_command) + "'")
        zeo_process = subprocess.Popen(zeo_command, stdout=subprocess.PIPE)
        zeo_process.wait()

        res_file = open('%s.res' % job_name).read().split()
        self.structure.pore_diameter = tuple(float(x) for x in res_file[1:])

        atom_samples = '%i' % 2000
        volume_samples = '%i' % (20*self.structure.cell.volume)

        for probe in probes:
            zeo_command = zeo_exe + [
                '-chan', '%f' % probe,
                '-sa', '%f' % probe, '%f' % probe, atom_samples,
                '-vol', '%f' % probe, '%f' % probe, volume_samples] + cssr_file

            debug("Running command: '" + " ".join(zeo_command) + "'")
            zeo_process = subprocess.Popen(zeo_command, stdout=subprocess.PIPE)
            zeo_process.wait()
            # channel dimensionality
            channels = [int(x) for x in open('%s.chan' % job_name).read().split()[5:]]
            self.structure.sub_property('dimensionality', probe, channels)
            # surface area
            for line in open('%s.sa' % job_name):
                if 'A^2' in line:
                    self.structure.sub_property('zeo_surface_area', probe,
                                                value=float(line.split()[-1]))
            # accessible volume
            for line in open('%s.vol' % job_name):
                if 'A^3' in line:
                    self.structure.sub_property('void_volume', probe,
                                                value=float(line.split()[-1]))


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
        self._esp_reduced = False
        if memory_guess > vmem:
            warning("ESP at this resolution might need up to %.1f GB of "
                    "memory but calculation will only request %.1f" %
                    (memory_guess, vmem))
            resolution = resolution/pow(vmem/memory_guess, 1.0/3)
            esp_grid = tuple([int(4*np.ceil(x/(4*resolution)))
                              for x in self.structure.cell.params[:3]])
            warning("Reduced grid to %.2f A resolution to fit" % resolution)
            self._esp_reduced = True

        self._esp_grid = esp_grid
        return esp_grid

    @property
    def esp_reduced(self):
        """Has the esp been reduced to fit the memory requirements?"""
        if not hasattr(self, '_esp_reduced'):
            # generate the esp and check memory requirements
            self.esp_grid
        return self._esp_reduced

    def calc_surface_area(self, rprobe=0.0):
        """Accessible surface area by uniform or Monte Carlo sampling."""
        self.structure.gen_neighbour_list()
        xyz = []
        resolution = self.options.getfloat('surface_area_resolution')
        uniform = self.options.getbool('surface_area_uniform_sample')
        info("Calculating surface area: %.3f probe, %s points, %.3f res" %
             (rprobe, ("random","uniform")[uniform], resolution))
        total_area = 0.0
        negative_charge = 0.0
        positive_charge = 0.0
        cell = self.structure.cell.cell
        inv_cell = np.linalg.inv(cell.T)
        # Pre-calculate and organise the in-cell atoms
        atoms = [(atom.ipos(cell, inv_cell).tolist(),
                  atom.ifpos(inv_cell),
                  atom.vdw_radius+rprobe,
                  atom.neighbours,
                  atom) for atom in self.structure.atoms]

        # sigma is the vdw_radius plus distance to center of probe, which
        # gives accessible surface area;
        all_samples = []
        for a1_pos, a1_fpos, a1_sigma, neighbours, atom in atoms:
            surface_area = 4*pi*(a1_sigma**2)
            nsamples = int(surface_area/resolution)
            if not nsamples in all_samples:
                debug("Atom type with %i samples" % nsamples)
                all_samples.append(nsamples)
            ncount = 0

            if uniform:
                # uniform spiral sample of surface
                z_vals = np.linspace(1, -1, nsamples, endpoint=True)
                r_vals = sqrt(1-z_vals**2)
                t_vals = np.linspace(0, pi*(3-(5**0.5))*nsamples,
                                     nsamples, endpoint=False)
                points = array([r_vals*cos(t_vals),
                                r_vals*sin(t_vals),
                                z_vals]).transpose()*a1_sigma + a1_pos
            else:
                # random MC sampling
                phi = 2*np.random.random(nsamples)*pi
                costheta = np.random.random(nsamples)*2 - 1
                theta = arccos(costheta)
                points = array([sin(theta)*cos(phi),
                                sin(theta)*sin(phi),
                                cos(theta)]).transpose()*a1_sigma + a1_pos

            # All points are brought into the cell
            points = [dot(inv_cell, x) for x in points]
            fpoints = np.mod(points, 1.0)
            points = [list(dot(x, cell)) for x in fpoints]

            for point, fpoint in zip(points, fpoints):
            # Check for overlap
                for a2_dist, a2_idx in neighbours:
                    a2_pos = atoms[a2_idx][0]
                    a2_fpos = atoms[a2_idx][1]
                    a2_sigma = atoms[a2_idx][2]
                    if a2_dist > a1_sigma + a2_sigma:
                        # No more atoms within the radius, point valid
                        ncount += 1
                        xyz.append((atom.type, point, atom.charge))
                        if atom.charge < 0:
                            negative_charge += atom.charge
                        else:
                            positive_charge += atom.charge
                        break
                    elif vecdist3(point, a2_pos) < a2_sigma:
                        # Point collision
                        break
                    elif min_dist(point, fpoint, a2_pos, a2_fpos, cell) < a2_sigma:
                        # Periodic collision
                        break
                else:
                    # Loop over all atoms finished; point valid
                    ncount += 1
                    xyz.append((atom.type, point, atom.charge))
                    if atom.charge < 0:
                        negative_charge += atom.charge
                    else:
                        positive_charge += atom.charge

            # Fraction of the accessible surface area for sphere to real area
            total_area += (surface_area*ncount)/nsamples
        if self.options.getbool('surface_area_save'):
            job_name = self.options.get('job_name')
            xyz_out = open('%s-surf-%.2f.xyz' % (job_name, rprobe), 'w')
            xyz_out.write('%i\nResolution: %f Area: %f\n' %
                          (len(xyz), resolution, total_area))
            for ppt in xyz:
                xyz_out.write(('%-6s' % ppt[0]) +
                              ('%10.6f %10.6f %10.6f' % tuple(ppt[1])) +
                              ('%10.6f\n' % ppt[2]))
        info("Total positive charge: %f" % positive_charge)
        info("Total negative charge: %f" % negative_charge)
        return total_area

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
        self.space_group = None

    def from_file(self, basename, filetype, defaults):
        """Select the correct file parser."""
        if filetype in ['sql', 'sqlite']:
            # Makeshift selection method to select a mof
            # jobname is dbname.type.structure_id
            from backend.sql import AlchemyBackend
            # [db_name, sym or free, identity]
            db_params = self.name.split('.')
            reader = AlchemyBackend(db_params[0])
            cif_string = reader.start_cif(db_params[1], int(db_params[2]))
            self.from_cif(string=cif_string)
        elif filetype.lower() in ['pdb']:
            self.from_pdb(basename + '.' + filetype)
        elif filetype.lower() in ['pqr']:
            # Look for a pqr or just a pdb wih charges
            listdir = os.listdir('.')
            if (basename + '.pqr') in listdir:
                self.from_pdb(basename + '.pqr', charges=True)
            else:
                self.from_pdb(basename + '.pdb', charges=True)
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
        opt_path = path.join('faps_%s_%s' % (self.name, opt_code))
        info("Updating positions from %s" % opt_code)
        if opt_code == 'vasp':
            self.from_vasp(path.join(opt_path, 'CONTCAR'), update=True)
        elif opt_code == 'siesta':
            self.from_siesta(path.join(opt_path, '%s.STRUCT_OUT' % self.name))
        elif opt_code == 'gulp':
            opt_path = "%s_opt" % opt_path
            self.from_gulp_output(path.join(opt_path, '%s.grs' % self.name))
        else:
            error("Unknown positions to import %s" % opt_code)

    def update_charges(self, charge_method, options=None):
        """Select the method for updating charges."""
        charge_path = path.join('faps_%s_%s' % (self.name, charge_method))
        if charge_method == 'repeat':
            info("Updating charges from repeat")
            self.charges_from_repeat(
                path.join(charge_path, 'faps-%s.out' % self.name),
                options.getbool('symmetry'))
            # Cleanup of REPEAT files
            unneeded_files = options.gettuple('repeat_delete_files')
            remove_files(unneeded_files, charge_path)
            keep_files = options.gettuple('repeat_compress_files')
            compress_files(keep_files, charge_path)
        elif charge_method == 'gulp':
            info("Updating charges from GULP QEq")
            self.charges_from_gulp(
                path.join(charge_path, 'faps-%s.out' % self.name))
        elif charge_method == 'egulp':
            info("Updating charges from EGULP QEq")
            self.charges_from_egulp(path.join(charge_path, 'charges.dat'))
        else:
            error("Unknown charge method to import %s" % charge_method)

    def update_gcmc(self, tp_point, options):
        """Select the source for GCMC results and import."""
        gcmc_code = options.get('mc_code')
        gcmc_path = path.join('faps_%s_%s' % (self.name, gcmc_code))
        # Runs in subdirectories
        tp_path = path.join(gcmc_path, 'T%s' % tp_point[0] +
                               ''.join(['P%.2f' % x for x in tp_point[1]]))
        if gcmc_code == 'fastmc':
            info("Importing results from FastGCMC")
            self.fastmc_postproc(tp_path, tp_point, options)
        else:
            error("Unknown gcmc method to import %s" % gcmc_code)

    def from_pdb(self, filename, charges=False):
        """Read an initial structure from a pdb file."""
        info("Reading positions from pdb file: %s" % filename)
        filetemp = open(filename)
        pdb_file = filetemp.readlines()
        filetemp.close()
        # Build a local list before setting attribute
        newatoms = []
        for line in pdb_file:
            lline = line.lower()
            if lline.startswith('cryst1'):
                self.cell.from_pdb(line)
                self.space_group = line[55:56]
            elif lline.startswith('atom') or lline.startswith('hetatm'):
                newatom = Atom()
                newatom.from_pdb(line, charges=charges)
                newatoms.append(newatom)

        self.atoms = newatoms
        self.order_by_types()

    def from_cif(self, filename=None, string=None):
        """Genereate structure from a .cif file."""
        if filename is not None:
            info("Reading positions from cif file: %s" % filename)
            filetemp = open(filename)
            cif_file = filetemp.readlines()
            filetemp.close()
        elif string is not None:
            info("Positions from cif string")
            cif_file = string.splitlines()
        else:
            error("No source for cif file")
        cif_file = strip_blanks(cif_file)
        params = [None, None, None, None, None, None]
        atoms = []
        cif_bonds = {}
        symmetry = []
        loops = []
        idx = 0
        while idx < len(cif_file):
            # line text needs to be manageable; cif guidelines can be
            # permissive
            # Can no longer just check for _underscores in lines
            # as UFF types can have them and mess up parsing
            line = cif_file[idx].lower().strip()
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
            elif '_symmetry_space_group_name_h-m' in line:
                self.space_group = line.split()[1]
            elif 'loop_' in line:
                # loops for _atom_site, _symmetry and _geom
                heads = []
                body = []
                while line.startswith('_') or 'loop_' in line:
                    # must keep the loop_ line as this can still contain headers
                    heads.extend(line.split())
                    idx += 1
                    # don't lower these to keep atomic symbols
                    line = cif_file[idx].strip()
                while idx < len(cif_file) and not line.startswith('_') and not 'loop_' in line:
                    # shlex keeps 'quoted items' as one
                    # Some cifs seem to have primed atom symbols
                    # posix=False should help
                    # using .shlex instead of .split works with '#' comments too
                    split_line = shlex.shlex(line, posix=False)
                    split_line.whitespace_split = True
                    split_line = list(split_line)
                    body.extend([x.strip("'").strip('"') for x in split_line])
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
            if '_ccdc_geom_bond_type' in heads:
                while body:
                    bond_dict = dict(zip(heads, body))
                    # bond is sorted so there are no duplicates
                    # and tuple so it can be hashed
                    bond = (bond_dict['_geom_bond_atom_site_label_1'],
                            bond_dict['_geom_bond_atom_site_label_2'])
                    bond = tuple(sorted(bond))
                    # bond distance and type defualts to None if not specified
                    distance = bond_dict.get('_geom_bond_distance')
                    if distance is not None:
                        distance = float(distance)
                    bond_type = bond_dict.get('_ccdc_geom_bond_type')
                    cif_bonds[bond] = (distance, bond_type)
                    body = body[len(heads):]

        if not symmetry:
            debug('No symmetry found; assuming identity only')
            symmetry = [Symmetry('x,y,z')]

        newatoms = []
        for site_idx, atom in enumerate(atoms):
            for sym_op in symmetry:
                newatom = Atom(parent=self)
                newatom.from_cif(atom, self.cell.cell, sym_op, site_idx)
                newatoms.append(newatom)

        self.atoms = newatoms

        if len(symmetry) > 1:
            # can skip if just identity operation as it's slow for big systems
            # Some of pete's symmetrised mofs need a higher tolerence
            duplicate_tolerance = 0.2  # Angstroms
            self.remove_duplicates(duplicate_tolerance)
        self.order_by_types()

        bonds = {}
        # TODO(tdaff): this works for the one tested MOF; 0.1 was not enough
        # only check for bonds that are too long, not too short.
        bond_tolerence = 0.25
        # Assign bonds by index
        for bond, bond_data in cif_bonds.items():
            for first_index, first_atom in enumerate(self.atoms):
                if first_atom.site == bond[0]:
                    for second_index, second_atom in enumerate(self.atoms):
                        if second_atom is first_atom:
                            continue
                        if second_atom.site == bond[1]:
                            # TODO(tdaff): symmetry implementation for cif bonding
                            distance = min_distance(first_atom, second_atom)
                            bond_dist = bond_data[0]
                            if bond_dist is None:
                                bond_dist = first_atom.covalent_radius + second_atom.covalent_radius
                            if distance < (bond_dist + bond_tolerence):
                                # use the sorted index as bonds between the
                                # same type are doubly specified
                                bond_id = tuple(sorted((first_index, second_index)))
                                bonds[bond_id] = CCDC_BOND_ORDERS[bond_data[1]]
                                if first_atom.is_metal or second_atom.is_metal:
                                    first_atom.is_fixed = True
                                    second_atom.is_fixed = True

        self.bonds = bonds
        self.symmetry = symmetry

    def from_vasp(self, filename='CONTCAR', update=False):
        """Read a structure from a vasp [POS,CONT]CAR file."""
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

    def from_gulp_output(self, filename):
        """Update the structure from the gulp optimisation output."""
        info("Updating positions from file: %s" % filename)
        grs_out = open(filename)
        for line in grs_out:
            if line.strip() == 'cell':
                params = tuple(float(x) for x in grs_out.next().split()[:6])
                self.cell.params = params
            elif line.strip() == 'fractional':
                cell = self.cell.cell
                for atom, atom_line in zip(self.atoms, grs_out):
                    atom.pos = dot([gfloat(x) for x in atom_line.split()[2:5]], cell)

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

    def charges_from_repeat(self, filename, symmetry=False):
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
        if symmetry:
            tree = self.symmetry_tree
            if len(charges) != len(tree):
                error("Incorrect number of charge sets; check REPEAT output")
                terminate(97)
            for symm, charge in zip(sorted(tree.items()), charges):
                for at_idx in symm[1]:
                    self.atoms[at_idx].charge = charge[2]
        else:
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
        start_line = None
        failures = 0
        for line_num, line in enumerate(gout):
            if '  Final charges from QEq :' in line:
                # Take the last set of charges
                start_line = line_num + 7
            elif 'Failed to converge' in line:
                # This will occur twice in the output for complete failure
                failures += 1
        if start_line is None:
            error("Final charges not found in gulp output")
            terminate(184)
        elif failures > 1:
            warn("Gulp charges may not be converged")
        for atom, chg_line in zip(self.atoms, gout[start_line:]):
            atom.charge = float(chg_line.split()[2])

    def charges_from_egulp(self, filename):
        """Parse QEq charges from EGULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        # charges.dat file only has list of charges in it
        for atom, chg_line in zip(self.atoms, gout):
            atom.charge = float(chg_line.split()[2])
            if atom.charge != atom.charge:
                # These are 'nan'; should not continue or fastmc will mess up
                error("Egulp charges did not converge, check structure")
                terminate(107)
            elif abs(atom.charge) == INFINITY:
                error("Egulp gave infinite charges, check structure")
                terminate(108)
            elif abs(atom.charge) > 10:
                warn("Very high charge from egulp: %s %f"
                     (atom.site, atom.charge))


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

        if "In" in u_types:
            # semicore states need to be accounted for
            # TODO(tdaff): selective polarisation on the basis functions
            fdf.extend(["\n%block PAO.Basis\n",
                        "In   3\n",
                        "  n=5  0  2\n",
                        "    0.0 0.0\n",
                        "  n=5  1  2  P\n",
                        "    0.0 0.0\n",
                        "  n=4  2  2\n",
                        "    0.0 0.0\n",
                        "%endblock PAO.Basis\n"])

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

    def to_gulp(self, qeq_fit=False, optimise=False, terse=False, qeq_dict={}):
        """Return a GULP file to use for the QEq charges."""
        if qeq_fit:
            from elements import QEQ_PARAMS
            keywords = "fitting bulk_noopt qeq\n"
        elif optimise:
            return self.to_gulp_optimise(terse=terse)
        else:
            keywords = "single conp qeq\n"

        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# Initial file written by faps\n"
            "name %s\n" % self.name,
            "vectors\n"] + self.cell.to_vector_strings() + [
            "cartesian\n"]
        if qeq_fit:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                                 "%14.7f\n" % atom.charge])
            gin_file.append("\n# \n# Fitting variables:\n# \nobservables\n")
            for idx, atom in enumerate(self.atoms):
                gin_file.append("monopoleq\n%5s%11.6f\n" % (idx+1, atom.charge))
            gin_file.append("end\nqelectronegativity\n")
            for at_type in unique(self.types):
                gin_file.extend(
                    ["%-5s" % at_type,
                     "%9.5f %9.5f %9.5f 1 1 1\n" % QEQ_PARAMS[at_type]])
        else:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f\n" % tuple(atom.pos)])
            if qeq_dict:
                unique_types = unique(self.types)
                gin_file.append('\nqelectronegativity\n')
                for atom_type, params in qeq_dict.items():
                    if atom_type in unique_types:
                        gin_file.append('%-4s %f %f\n' % (atom_type, params[0], params[1]))

        gin_file.append("\ndump every %s.grs\nprint 1\n" % self.name)
        return gin_file


    def to_gulp_optimise(self, terse=False):
        """Return a GULP file to optimise with UFF."""
        # all atoms of the same forcefield type must have the same label
        # gulp fails if atoms of the same type have different labels!
        # this only crops up as an error in the 4.0 versions
        # 'decimal_only' stops fractions that cannot be parsed when reading
        # updated positions
        if terse:
            # Don't output the bonds
            keywords = "opti noautobond decimal_only\n"
        else:
            keywords = "opti noautobond bond decimal_only\n"
        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# UFF optimisation by gulp\n"
            "name %s\n" % self.name,
            # using an rfo minimiser with a preconditioned hessian from the
            # BFGS minimiser seems to be the most efficient
            "switch rfo gnorm 0.1\n",
            "vectors\n"] + self.cell.to_vector_strings() + [
            " 1 1 1\n 1 1 1\n 1 1 1\n",  # constant pressure relaxation
            "cartesian\n"]
        all_ff_types = {}
        #TODO(tdaff): warn for no ff types
        for at_idx, atom in enumerate(self.atoms):
            ff_type = atom.uff_type
            if not ff_type in all_ff_types:
                all_ff_types[ff_type] = atom.site
            if atom.is_fixed:
                fixed_flags = "0 0 0"
            else:
                fixed_flags = "1 1 1"
            gin_file.extend(["%-5s core " % all_ff_types[ff_type],
                             "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                             "%f " % atom.charge,
                             "%f " % 1.0,  # occupancy
                             fixed_flags, "\n"])

        #identify all the individual uff species for the library
        gin_file.append("\nspecies\n")
        for ff_type, species in all_ff_types.items():
            gin_file.append("%-6s core %-6s\n" % (species, ff_type))

        gin_file.append("\n")
        for bond in sorted(self.bonds):
            bond_type = GULP_BOND_ORDERS[self.bonds[bond]]
            gin_file.append("connect %6i %6i %s\n" % (bond[0] + 1, bond[1] + 1, bond_type))

        gin_file.append("\nlibrary uff\n")

        # Restart file is for final structure
        gin_file.append("\ndump every %s.grs\n" % self.name)
        # optimization movie useful for debugging mostly
        if terse:
            # These tell gulp to be quiet, but we also stop the massive arc
            # file being generated
            gin_file.append("\nterse inout structure\n"
                            "terse inout potentials\n"
                            "terse inout derivatives\n"
                            "#output movie arc %s\n" % self.name)
        else:
            gin_file.append("\noutput movie arc %s\n" % self.name)

        return gin_file


    def to_egulp(self):
        """Generate input files for Eugene's QEq code."""
        # bind cell locally for speed and convenience
        cell = self.cell.cell
        inv_cell = self.cell.inverse
        # format exactly as Eugene original script generates it
        geometry_file = ['%s\n' % self.name]
        geometry_file.extend(self.cell.to_vector_strings(fmt=' %15.12f'))
        geometry_file.append('%i\n' % self.natoms)
        geometry_file.extend([
            ('%6d ' % atom.atomic_number) +
            ('%12.7f %12.7f  %12.7f\n' % tuple(atom.ipos(cell, inv_cell)))
            for atom in self.atoms
        ])
        return geometry_file


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

    def to_gromacs(self):
        """Procedure:
        generate .gro atom positions
        generate .top lj topology
        generate .itp uff topology
        gererate .mdp parameters
        grompp -f mdp_file -c gro_file -p top_file -o output.tpr
        mdrun -s moffive -o traj-c moffive -nt 1
        """

        ##
        # .gro file
        ##
        # TITLE linw
        # number of atoms
        gro_file = ["%s\n" % self.name, " %i\n" % self.natoms]
        # TODO(tdaff): check residue in toplogy?
        residue = (1, self.name[:4].upper())
        for idx, atom in enumerate(self.atoms):
            # In the specification the atom position is given as
            # %8.3f in nm, but we probably would like more accuracy
            # but make sure there are 5 figures to left of decimal point
            pos_nm = tuple([x/10.0 for x in atom.pos])
            gro_file.append(("%5i%5s" % residue) +
                            ("%5s%5i" % (atom.uff_type, idx + 1)) +
                            ("%15.10f%15.10f%15.10f\n" % pos_nm))

        # cell is also in nm
        # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        box = self.cell.cell/10.0
        # gromacs needs these condtiions. Might need to rotate cell sometimes?
        if not (box[0][1] == box[0][2] == box[1][2] == 0):
            error("Gromacs can't handle this cell orientation")
        gro_file.append("%f %f %f %f %f %f %f %f %f\n" % (
            box[0][0], box[1][1], box[2][2],
            box[0][1], box[0][2],
            box[1][0], box[1][2],
            box[2][0], box[2][1]))

        ##
        # .top file
        ##
        comb_rule = 1 # 1 = LORENTZ_BERTHELOT; 2 = GEOMETRIC
        top_file = [
            "; Topology file for %s\n" % self.name,
            "[ defaults ]\n",
            "; nbfunc      comb-rule      gen-pairs       fudgeLJ    fudgeQQ\n",
            "1             %i              yes             1.0        1.0\n\n" % (comb_rule)]

        # define the forcefield for UFF
        unique_types = set()
        top_file.append("[ atomtypes ]\n")
        top_file.append("; name1 name2   mass     charge  ptype   sigma   epsilon\n")
        for atom in self.atoms:
            if atom.uff_type not in unique_types:
                uff_type = atom.uff_type
                # sigma = x1 * 2^(-1/6)
                # CONVERT TO nm !!
                sigma = 0.1*UFF_FULL[uff_type][2]*(2**(-1.0/6.0))
                # epsilon = D1 in kcal
                epsilon = UFF_FULL[uff_type][3] * KCAL_TO_KJ
                top_file.append("%-6s   %-6s   %9.4f   %9.4f   A %12.7f %12.7f\n" % (uff_type, uff_type, atom.mass, 0.0, sigma, epsilon))
                unique_types.add(uff_type)

        top_file.extend(["\n#include <%s.itp>\n\n" % self.name,
                         "[ system ]\n",
                         "UFF optimisation of %s\n\n" % self.name,
                         "[ molecules ]\n",
                         "%s  1\n" % residue[1]])

        ##
        # .itp file
        ##
        # exclude 3 neighbours
        itp_file = ["[ moleculetype ]\n",
                    "; molname nrexcl\n"
                    "%s 3\n" % residue[1],
                    "[ atoms ]\n",
                    "; nr type  resnr    residue    atom     cgnr    charge       mass \n"]

        # atoms
        for idx, atom in enumerate(self.atoms):
            uff_type = atom.uff_type
            # charge group is different for each atom as gromacs has max of 32 in a group
            itp_file.append("%-6i  %-6s  %i  %-6s  %-6s   %d  %9.4f  %9.4f\n" % (idx+1, uff_type, residue[0], residue[1], uff_type, idx+1, atom.charge, atom.mass))

        # bonds
        itp_file.append(" [ bonds ]\n")
        itp_file.append("; ai aj funct b0 kb\n")

        unique_bonds = {}
        bonding_table = dict([(x, {}) for x in range(self.natoms)])
        #FIXME(tdaff) bonds will have length in v2
        for bond, bondorder in self.bonds.items():
            idx_a = bond[0]
            idx_b = bond[1]
            atom_a = self.atoms[idx_a]
            atom_b = self.atoms[idx_b]
            uff_a = atom_a.uff_type
            uff_b = atom_b.uff_type
            typed_bond = tuple(sorted([uff_a, uff_b]) + [bondorder])

            bonding_table[idx_a][idx_b] = typed_bond
            bonding_table[idx_b][idx_a] = typed_bond

            # have we already calculated the parameters
            if not typed_bond in unique_bonds:
                ri = UFF_FULL[uff_a][0]
                rj = UFF_FULL[uff_b][0]
                chiI = UFF_FULL[atom_a.uff_type][8]
                chiJ = UFF_FULL[atom_b.uff_type][8]
                rbo = -0.1332*(ri+rj)*log(bondorder)

                ren = ri*rj*(((sqrt(chiI) - sqrt(chiJ))**2)) / (chiI*ri + chiJ*rj)
                r0 = 0.1*(ri + rj + rbo - ren)

                # force constant
                # parameters Z1
                kb = (KCAL_TO_KJ * 664.12 * UFF_FULL[uff_a][5] * UFF_FULL[uff_b][5])/(r0**3)

                unique_bonds[typed_bond] = (r0, kb)

            params = unique_bonds[typed_bond]
            bond_func = 1  # gromacs harmonic
            # add 1 to bond as 1 indexed
            itp_file.append('%-5i %-5i %1i %11.4f %11.4f ; %-5s %-5s %.2f\n' % (bond[0]+1, bond[1]+1, bond_func, params[0], params[1], typed_bond[0], typed_bond[1], bondorder))

        # angles
        itp_file.append("\n [ angles ]\n")

        for idx_a in sorted(bonding_table):
            bonded_atoms = bonding_table[idx_a]
            for l_idx in bonded_atoms:
                for r_idx in bonded_atoms:
                    if l_idx == r_idx:
                        continue
                    # FIXME(tdaff) special cases for 5 or more coord
                    central_atom = self.atoms[idx_a]
                    l_atom = self.atoms[l_idx]
                    r_atom = self.atoms[r_idx]
                    coordination = len(bonded_atoms)
                    # FIXME(tdaff) over/under coordinated atoms?
                    theta0 = UFF_FULL[central_atom.uff_type][1]
                    cosT0 = cos(theta0*DEG2RAD)
                    sinT0 = sin(theta0*DEG2RAD)
                    c2 = 1.0 / (4.0 * sinT0 * sinT0)
                    c1 = -4.0 * c2 * cosT0
                    c0 = c2*(2.0*cosT0*cosT0 + 1.0)
                    zi = UFF_FULL[l_atom.uff_type][5]
                    zk = UFF_FULL[r_atom.uff_type][5]
                    bond_ab = tuple(sorted((l_idx, idx_a)))
                    bond_ab = tuple(sorted([l_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_ab]])
                    bond_bc = tuple(sorted((r_idx, idx_a)))
                    bond_bc = tuple(sorted([r_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_bc]])
                    rab = unique_bonds[bond_ab][0]
                    rbc = unique_bonds[bond_bc][0]
                    rac = sqrt(rab*rab + rbc*rbc - 2.0 * rab*rbc*cosT0)

                    ka = (644.12 * KCAL_TO_KJ) * (zi * zk / (rac**5.0))
                    ka *= (3.0*rab*rbc*(1.0 - cosT0*cosT0) - rac*rac*cosT0)

                    thetamin = (pi-arccos(c1/(4.0*c2))/DEG2RAD)
                    kappa = ka * (16.0*c2*c2 - c1*c1) / (4.0* c2)
                    itp_file.append("%-5i %-5i %-5i %3i %f %f ; %-5s %-5s %-5s" % (l_idx + 1, idx_a + 1, r_idx + 1, 1, thetamin, kappa, l_atom.uff_type, central_atom.uff_type, r_atom.uff_type))

        with open("moffive.top", 'w') as tempfile:
            tempfile.writelines(top_file)
        with open("moffive.gro", 'w') as tempfile:
            tempfile.writelines(gro_file)
        with open("moffive.itp", 'w') as tempfile:
            tempfile.writelines(itp_file)


    def to_cssr(self, cartesian=False, no_atom_id=False):
        """
        Return a Cerius2 cssr file with coordinates and cell as a list of
        strings. Set no_atom_id to produce labels to work with Zeo++.

        """

        space_group = (1, "P 1")
        opt = 1
        cell = self.cell
        # spacers between all format strings ensures that there will always
        # be whitespace between components.
        cssr = [" "*38, "%8.3f %7.3f %7.3f\n" % cell.params[:3],
                " "*21, "%8.3f %7.3f %7.3f" % cell.params[3:], " "*4,
                "SPGR = %3i %-11s" % space_group, "OPT = %i\n" % opt,
                "%4i %3i %s\n" % (self.natoms, cartesian, self.name),
                "Structure file generated by faps\n"]

        for at_idx, atom in enumerate(self.atoms):
            # The default of TypeID does not work for Zeo++ as it
            # treats each label as a different type
            if no_atom_id:
                atom_name = "%s" % (atom.type)
            else:
                atom_name = "%s%i" % (atom.type, at_idx+1)
            # Zeo++ needs fractional coordinates
            if cartesian:
                string_pos = "%9.5f %9.5f %9.5f " % tuple(atom.pos)
            else:
                string_pos = "%9.5f %9.5f %9.5f " % tuple(atom.ifpos(cell.inverse))
            cssr.append("%4i %-4s  %s" % (at_idx+1, atom_name, string_pos) +
                        "   0"*8 + " %7.3f\n" % atom.charge)

        return cssr

    def to_zeoplusplus(self):
        """
        Return a tuple containing a cssr file, radii file and mass file.

        """

        radii = ["%-7s %-f\n" % (atom, UFF[atom][0]/2.0)
                 for atom in unique(self.types)]
        masses = ["%-7s %-f\n" % (atom, WEIGHT[atom])
                  for atom in unique(self.types)]

        return (self.to_cssr(no_atom_id=True), radii, masses)

    def fastmc_postproc(self, filepath, tp_point, options):
        """Update structure properties from gcmc OUTPUT."""
        startdir = os.getcwd()
        os.chdir(filepath)

        # Since Pete changed output strip indentation and blank lines
        filetemp = open('OUTPUT')
        output = strip_blanks(filetemp.readlines())
        filetemp.close()

        # Keep track of supercell so we can get unit cell values
        supercell_mult = prod(self.gcmc_supercell)
        # Still positional as we need multiple values simultaneously
        # and very old versions cahnged wording of heat of adsorption
        # and enthalpy of guest
        for idx, line in enumerate(output):
            # Assume that block will always start like this
            if 'final stats' in line:
                guest_id = int(line.split()[4]) - 1
                self.guests[guest_id].uptake[tp_point] = (
                    float(output[idx + 1].split()[-1]),
                    float(output[idx + 2].split()[-1]),
                    supercell_mult)
                # This will sometimes be NaN
                self.guests[guest_id].hoa[tp_point] = (
                    float(output[idx + 3].split()[-1]),
                    float(output[idx + 4].split()[-1]))
            elif 'total accepted steps' in line:
                counted_steps = int(line.split()[-1])
                if counted_steps < 10000:
                    warning("Number of accepted GCMC steps is very low; "
                            "only %i counted!" % counted_steps)

        fold = options.getbool('fold')
        find_maxima = options.getbool('find_maxima')
        prob_plot = options.getbool('mc_probability_plot')
        folded = False
        if prob_plot and (fold or find_maxima):
            folded = self.fold_and_maxima(fold, find_maxima, tp_point)

        if folded and not options.getbool('fastmc_keep_unfolded_cubes'):
            debug("Removing unfolded cube files")
            cubes = glob.glob("prob_guest??_prob_??.cube")
            remove_files(cubes)

        #TODO(tdaff): absl

        unneeded_files = options.gettuple('fastmc_delete_files')
        remove_files(unneeded_files)
        keep_files = options.gettuple('fastmc_compress_files')
        compress_files(keep_files)

        os.chdir(startdir)


    def fold_and_maxima(self, fold=True, find_maxima=True, tp_point=None):
        """Determine the positions of maxima and produce an xyz xyz file."""
        from cube import Cube
        folded = False
        if fold:
            fold = self.gcmc_supercell
        else:
            fold = None
        for guest_idx, guest in enumerate(self.guests):
            guest_locations = {}
            for site_idx, sites in enumerate(guest.probability):
                guest_cube = Cube("prob_guest%02i_prob_%02i.cube" %
                                  (guest_idx+1, site_idx+1), fold=fold)
                if fold is not None:
                    debug("Folded cube file: %s" % guest_cube.folded_name)
                    guest_cube.write_cube()
                    folded = True
                if find_maxima:
                    guest_locations[sites] = guest_cube.maxima()
            if guest_locations:
                if tp_point:
                    # We can keep them for later too, must create dict
                    # since it might not exist in old calculations
                    if not hasattr(guest, 'guest_locations'):
                        guest.guest_locations = {tp_point: guest_locations}
                    else:
                        guest.guest_locations[tp_point] = guest_locations
                maxima = []
                for sites in guest_locations:
                    atom_name = name_from_types(sites, guest)
                    for atom, magnitude in guest_locations[sites]:
                        maxima.append((magnitude, "%-6s" % atom_name +
                                      "%10.6f %10.6f %10.6f " % tuple(atom) +
                                      "%10.6f\n" % magnitude))
                locations = open("%s-%s.xyz" % (self.name, guest.ident), 'w')
                locations.write(" %i\nBinding sites at %r\n" %
                                (len(maxima), tp_point))
                locations.writelines([x[1] for x in sorted(maxima, reverse=True)])
                locations.close()

        return folded


    def remove_duplicates(self, tolerance=0.02):
        """Find overlapping atoms and remove them."""
        uniq_atoms = []
        found_atoms = []
        for atom in self.atoms:
            for uniq_atom in uniq_atoms:
                if atom.type != uniq_atom.type:
                    continue
                elif min_distance(atom, uniq_atom) < tolerance:
                    break
            # else excutes when not found here
            else:
                uniq_atoms.append(atom)
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
        """
        Iterate over all the atoms of supercell where scale is an integer
        to scale uniformly or triplet with scale factors for each direction.

        """
        # Beware supercells larger than 2147483647 are not supported in
        # python 2
        if isinstance(scale, int):
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

    def gen_neighbour_list(self, force=False):
        """All atom pair distances."""
        # This can be expensive so skip if already calcualted
        if not force:
            for atom in self.atoms:
                if not hasattr(atom, 'neighbours'):
                    break
            else:
                # finished loop over all atoms
                debug("Neighbour list already calculated")
                return

        debug("Calculating neighbour list.")

        cell = self.cell.cell
        inv_cell = self.cell.inverse
        cpositions = [atom.ipos(cell, inv_cell) for atom in self.atoms]
        fpositions = [atom.ifpos(inv_cell) for atom in self.atoms]
        cell = cell.tolist()

        # loop over all pairs to find minimum periodic distances
        for atom, a_cpos, a_fpos in zip(self.atoms, cpositions, fpositions):
            neighbours = []
            for o_idx, o_cpos, o_fpos in zip(count(), cpositions, fpositions):
                sep = min_dist(a_cpos, a_fpos, o_cpos, o_fpos, cell)
                neighbours.append((sep, o_idx))
            # First one is self == 0
            # save in incresaing distance order
            atom.neighbours = sorted(neighbours)[1:]

        ## Neighbourlist printed in VASP style
        #for idx, atom in enumerate(self.atoms):
        #    print("%4i" % (idx+1) +
        #          "%7.3f%7.3f%7.3f" % tuple(atom.ifpos(inv_cell)) +
        #          "-" +
        #          "".join("%4i%5.2f" % (y+1, x) for x, y in atom.neighbours if x<2.5))

    def surface_area(self, probe=None, value=None, delete=False):
        """
        Helper:
          Return all {probe:area} if no arguments given
          Return the area or None for a given probe
          Set area if value given
          Delete value if delete is True
        Areas in A^2
        """
        surface_areas = self.properties.get('surface_area', {})
        if value is not None:
            surface_areas[probe] = value
            self.properties['surface_area'] = surface_areas
        elif delete:
            # Set it to None to avoid KeyErrors
            surface_areas[probe] = None
            del surface_areas[probe]
            self.properties['surface_area'] = surface_areas
        elif probe is not None:
            return surface_areas.get(probe, None)
        else:
            return surface_areas

    def sub_property(self, name, probe=None, value=None, delete=False):
        """
        Helper:
          Return all {probe:value} if no arguments given
          Return the value or None for a given probe
          Set area if value given
          Delete value if delete is True
        Units are based on Angstrom
        """
        property_data = self.properties.get(name, {})
        if value is not None:
            property_data[probe] = value
            self.properties[name] = property_data
        elif delete:
            # Set it to None to avoid KeyErrors
            property_data[probe] = None
            del property_data[probe]
            self.properties[name] = property_data
        elif probe is not None:
            return property_data.get(probe, None)
        else:
            return property_data

    def void_volume(self):
        """Estimate the void volume based on VdW radii."""
        initial_resolution = 0.2
        params = self.cell.params
        cell = self.cell.cell
        inv_cell = np.linalg.inv(cell.T)

        grid_size = [int(ceil(x/initial_resolution)) for x in params[:3]]
        print(grid_size)
        grid_resolution = [params[0]/grid_size[0],
                           params[1]/grid_size[1],
                           params[2]/grid_size[2]]
        print(grid_resolution)
        grid = np.zeros(grid_size, dtype=bool)
        atoms = [(atom.ipos(cell), inv_cell.tolist(),
                  atom.ifpos(inv_cell),
                  atom.vdw_radius) for atom in self.atoms]

        for x_idx in range(grid_size[0]):
            print(x_idx)
            for y_idx in range(grid_size[1]):
                print(y_idx)
                print(grid.sum())
                for z_idx in range(grid_size[2]):
                    grid_pos = [x_idx*grid_resolution[0],
                                y_idx*grid_resolution[1],
                                z_idx*grid_resolution[2]]
                    grid_fpos = [float(x_idx)/grid_size[0],
                                 float(y_idx)/grid_size[1],
                                 float(z_idx)/grid_size[2]]
                    for pos, fpos, radius in atoms:
                        dist = min_dist(pos, fpos, grid_pos, grid_fpos, cell)
                        if dist < radius:
                            grid[x_idx, y_idx, z_idx] = 1
                            break

        print(grid.sum())

    @property
    def types(self):
        """Ordered list of atom types."""
        return [atom.type for atom in self.atoms]

    @property
    def atomic_numbers(self):
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

    @property
    def symmetry_tree(self):
        """Tree of atoms that are symmetrically equivalent."""
        tree = {}
        for atom_id, atom in enumerate(self.atoms):
            if atom.site in tree:
                tree[atom.site].append(atom_id)
            else:
                tree[atom.site] = [atom_id]
        if len(tree) == 1 and None in tree:
            return dict((i, [i]) for i in range(self.natoms))
        else:
            return tree

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
        self._inverse = None

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
        # If the supercell is more than 2147483647 in any direction this
        # will fail in python 2, but 'long' removed for py3k forward
        # compatibility
        if isinstance(scale, int):
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
    def crystal_system(self):
        """Return the IUCr designation for the crystal system."""
        #FIXME(tdaff): must be aligned with x to work
        if self.alpha == self.beta == self.gamma == 90:
            if self.a == self.b == self.c:
                return 'cubic'
            elif self.a == self.b or self.a == self.c or self.b == self.c:
                return 'tetragonal'
            else:
                return 'orthorhombic'
        elif self.alpha == self.beta == 90:
            if self.a == self.b and self.gamma == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.alpha == self.gamma == 90:
            if self.a == self.c and self.beta == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.beta == self.gamma == 90:
            if self.b == self.c and self.alpha == 120:
                return 'hexagonal'
            else:
                return 'monoclinic'
        elif self.a == self.b == self.c and self.alpha == self.beta == self.gamma:
            return 'trigonal'
        else:
            return 'triclinic'

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
        # Class internally expects an array
        self._cell = array(value).reshape((3,3))
        self.__mkparam()
        self._inverse = np.linalg.inv(self.cell.T)

    # Property so that params are updated when cell is set
    cell = property(get_cell, set_cell)

    def get_params(self):
        """Get the six parameter cell representation as a tuple."""
        return tuple(self._params)

    def set_params(self, value):
        """Set cell and params from the cell parameters."""
        self._params = value
        self.__mkcell()
        self._inverse = np.linalg.inv(self.cell.T)

    # Property so that cell is updated when params are set
    params = property(get_params, set_params)

    @property
    def inverse(self):
        """Inverted cell matrix for converting to fractional coordinates."""
        try:
            if self._inverse is None:
                self._inverse = np.linalg.inv(self.cell.T)
        except AttributeError:
            self._inverse = np.linalg.inv(self.cell.T)
        return self._inverse

    @property
    def a(self):
        """Magnitude of cell a vector."""
        return self.params[0]

    @property
    def b(self):
        """Magnitude of cell b vector."""
        return self.params[1]

    @property
    def c(self):
        """Magnitude of cell c vector."""
        return self.params[2]

    @property
    def alpha(self):
        """Cell angle alpha."""
        return self.params[3]

    @property
    def beta(self):
        """Cell angle beta."""
        return self.params[4]

    @property
    def gamma(self):
        """Cell angle gamma."""
        return self.params[5]

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

    def __init__(self, at_type=False, pos=False, parent=None, **kwargs):
        """Accept arbritary kwargs as attributes."""
        if parent is not None:
            self._parent = parent
        self.type = at_type
        self.pos = pos
        self.charge = 0.0
        self.idx = False
        self.site = None
        self.mass = 0.0
        self.molecule = None
        self.uff_type = None
        self.is_fixed = False
        # Sets anything else specified as an attribute
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __str__(self):
        return "%s %f %f %f" % tuple([self.type] + list(self.pos))

    def __repr__(self):
        return "Atom(%r,%r)" % (self.type, self.pos)

    def fpos(self, inv_cell):
        """Fractional position within a given cell."""
        return dot(inv_cell, self.pos).tolist()

    def ifpos(self, inv_cell):
        """In cell fractional position."""
        return [i % 1 for i in self.fpos(inv_cell)]

    def ipos(self, cell, inv_cell):
        """In cell cartesian position."""
        return dot(self.ifpos(inv_cell), cell)

    def from_cif(self, at_dict, cell, symmetry=None, idx=None):
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
        if re.match('[0-9]', self.site) and idx is not None:
            debug("Site label may not be unique; appending index")
            self.site = "%s%i" % (self.site, idx)
        if '_atom_type_description' in at_dict:
            self.uff_type = at_dict['_atom_type_description']
        if '_atom_type_parital_charge' in at_dict:
            self.charge = float(at_dict['_atom_type_parital_charge'])

    def from_pdb(self, line, charges=False):
        """
        Parse the ATOM line from a pdb file.
        Occupancy field may be used to specify the charge as in a '.pqr' file.
        """
        # pdb is defined with fixed width fields rather than splitting
        self.idx = try_int(line[6:11])
        self.site = line[12:16].strip()
        self.molecule = try_int(line[22:26])
        at_pos = float(line[30:38]), float(line[38:46]), float(line[47:54])
        self.pos = at_pos
        self.type = line[76:78].strip()
        self.mass = WEIGHT[self.type]
        if charges:
            self.charge = float(line[54:60])

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
        split_line = line.split()
        self.pos = [float(x) for x in split_line[1:4]]
        if len(split_line) > 4:
            try:
                self.charge = float(split_line[4])
            except ValueError:
                # Assume a comment so skip it
                pass
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
                name = name[:-1]
        return atomic_number

    def set_atomic_number(self, value):
        """Set the atom type based on the atomic number."""
        self.type = ATOMIC_NUMBER[value]

    atomic_number = property(get_atomic_number, set_atomic_number)

    def get_fractional_coordinate(self):
        """Retrieve the fractional coordinates or calculate from the parent."""
        try:
            # Hopefully the attribute is just set
            return self._fractional
        except AttributeError:
            # Not set yet
            try:
                self._fractional = self.fpos(self._parent.cell.inverse)
                return self._fractional
            except AttributeError:
                return None

    def set_fractional_coordinate(self, value):
        """Set the position using the fractional coordinates."""
        fractional = [x % 1.0 for x in value]
        self._fractional = fractional
        self.pos = dot(fractional, self._parent.cell.cell)

    def del_fractional_coordinate(self):
        """Remove the fractional coordinate; run after updating cell."""
        del self._fractional

    fractional = property(get_fractional_coordinate, set_fractional_coordinate, del_fractional_coordinate)

    @property
    def element(self):
        """Guess the element from the type, fall back to type."""
        name = self.type
        while name:
            if name in ATOMIC_NUMBER:
                return name
            else:
                name = name[:-1]
        return self.type

    @property
    def vdw_radius(self):
        """Get the vdw radius from the UFF parameters."""
        return UFF[self.type][0]/2.0

    @property
    def covalent_radius(self):
        """Get the covalent radius from the library parameters."""
        return COVALENT_RADII[self.type]

    @property
    def is_metal(self):
        """Return True if element is in a predetermined set of metals."""
        return self.atomic_number in METALS


class Guest(object):
    """Guest molecule and properties."""
    def __init__(self, ident=None, guest_path=None):
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
            self.load_guest(ident, guest_path=None)

    def load_guest(self, ident, guest_path=None):
        """Look in guests.lib in submit directory and default."""
        # Ident set here to keep consistent
        self.ident = ident
        # Need the different directories
        if guest_path is None:
            guest_path = os.getcwd()
        if __name__ != '__main__':
            script_dir = path.dirname(__file__)
        else:
            script_dir = path.abspath(sys.path[0])
        dot_faps_dir = path.join(path.expanduser('~'), '.faps')
        # A parser for each location
        job_guests = configparser.SafeConfigParser()
        dot_faps_guests = configparser.SafeConfigParser()
        lib_guests = configparser.SafeConfigParser()
        # Try and find guest in guests.lib
        job_guests.read(path.join(guest_path, 'guests.lib'))
        dot_faps_guests.read(path.join(dot_faps_dir, 'guests.lib'))
        lib_guests.read(path.join(script_dir, 'guests.lib'))
        # Job dir and user defined have priority
        if job_guests.has_section(ident):
            debug("%s found in job dir" % ident)
            self._parse_guest(job_guests.items(ident))
        elif dot_faps_guests.has_section(ident):
            debug("%s found in dot_faps" % ident)
            self._parse_guest(dot_faps_guests.items(ident))
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
                self._molar_volume = float(val)
            elif key == 'probe radius':
                self.probe_radius = float(val)
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
        # TODO(tdaff): should we translate into cell?
        new_pos = [x%1.0 for x in new_pos]
        return new_pos


def mk_repeat(cube_name='REPEAT_ESP.cube', symmetry=False):
    """Standard REPEAT input file."""
    if symmetry:
        symmetry_flag = 1
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
        "Max distance cut-off scaling factor (1000.0 default)\n",
        "2.0\n"]

    filetemp = open('REPEAT_param.inp', 'w')
    filetemp.writelines(repeat_input)
    filetemp.close()


def mk_connectivity_ff(sym_tree):
    """Write connectivity.ff for input tree."""
    out_tree = ["%i\n" % len(sym_tree)]
    # Always sort to get consistent order
    for idx in sorted(sym_tree):
        out_tree.extend([
            "# %s tree\n" % str(idx),
            "%i\n" % len(sym_tree[idx]),
            " ".join(["%i" % (i+1) for i in sym_tree[idx]]),
            "\n"
        ])

    filetemp = open('connectivity.ff', 'w')
    filetemp.writelines(out_tree)
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
        incar.extend(["ENCUT   = 520\n",
                      "IBRION  = 2\n",
                      "NSW     = 800\n",
                      "ISIF    = 3\n"])
    elif optim_all or optim_h:
        # Just move positions
        incar.extend(["IBRION  = 2\n",
                      "NSW     = 600\n",
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

    # VASP recommends, for best performance:
    # NPAR = 4 - approx SQRT( number of cores)
    vasp_ncpu = options.getint('vasp_ncpu')
    npar = vasp_ncpu  # recommended up to 8 cpus
    if vasp_ncpu > 8:
        npar = 4*max(int((vasp_ncpu**0.5)/4.0), 1)
    incar.append("NPAR    = %i\n" % npar)

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


def mk_gcmc_control(temperature, pressures, options, guests, supercell=None):
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

    if options.getbool('fold') and supercell is not None:
        control.append('\ngrid factors %i %i %i\n' % supercell)
    elif supercell is not None:
        control.append('\n# grid factors %i %i %i\n' % supercell)

    resolution = options.getfloat('mc_probability_plot_spacing')
    if resolution != FASTMC_DEFAULT_GRID_SPACING:
        control.append('\ngrid spacing %f\n' % resolution)
    else:
        control.append('\n#grid spacing %f\n' % FASTMC_DEFAULT_GRID_SPACING)


    control.append("\nfinish\n")
    return control


def parse_qeq_params(param_tuple):
    """Convert an options tuple to a dict of values."""
    # group up into ((atom, electronegativity, 0.5*hardness), ... )
    param_tuple = subgroup(param_tuple, 3)
    param_dict = {}
    for param_set in param_tuple:
        try:
            atom_type = int(param_set[0])
            atom_type = ATOMIC_NUMBER[atom_type]
        except ValueError:
            # assume it is an element symbol instead
            atom_type = param_set[0]
        try:
            param_dict[atom_type] = (float(param_set[1]), float(param_set[2]))
        except IndexError:
            warn("Cannot read parameters for %s" % atom_type)

    return param_dict


def mk_egulp_params(param_tuple):
    """Convert an options tuple to an EGULP parameters file filling."""
    # group up into ((atom, electronegativity, 0.5*hardness), ... )
    param_dict = parse_qeq_params(param_tuple)
    # first line is the number of parametr sets
    #TODO(tdaff): try adding some error checking
    egulp_file = ["%s\n" % len(param_dict)]
    for atom_type, params in param_dict.items():
        egulp_file.append("%-4i %f %f\n" % (ATOMIC_NUMBER.index(atom_type),
                                            params[0], params[1]))

    return egulp_file



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
    if not path.exists(directory):
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
    if path.exists(dest):
        if path.isdir(dest):
            dest_full = path.join(dest, path.basename(src))
            if path.exists(dest_full):
                if path.isfile(dest_full):
                    os.remove(dest_full)
                    shutil.move(src, dest)
                else:
                    raise OSError("Directory %s already exists" % dest_full)
            else:
                shutil.move(src, dest)
        elif path.isfile(dest):
            os.remove(dest)
            shutil.move(src, dest)
        else:
            raise OSError("%s is not a folder or file" % dest)
    else:
        shutil.move(src, dest)


def try_symlink(src, dest):
    """Delete an existing dest file, symlink a new one if possible."""
    if path.lexists(dest):
        os.remove(dest)
    try:
        os.symlink(src, dest)
    except AttributeError:
        shutil.copy(src, dest)


def sys_argv_strip(argument):
    """Remove an argument from the sys.argv if it is there."""
    while argument in sys.argv:
        sys.argv.remove(argument)
    return 0


def same_guests(base, other):
    """Test if the guests are the same index and order."""
    return [item.ident for item in base] == [item.ident for item in other]


def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub('\(.*\)', '', text))


def gfloat(text):
    """Parse a gulp output float, where it could also be a rational fraction."""
    # TODO(tdaff): use Fraction in 2.7+?
    if "/" in text:
        text = text.split('/')
        return float(text[0])/float(text[1])
    else:
        return float(text)


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


def dot3(vec1, vec2):
    """Calculate dot product for two 3d vectors."""
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]


def min_distance(first_atom, second_atom, cell=None):
    """Helper to find mimimum image criterion distance."""
    if cell is None:
        cell = first_atom._parent.cell.cell
    return min_dist(first_atom.pos,
                    first_atom.fractional,
                    second_atom.pos,
                    second_atom.fractional,
                    cell)


def min_dist(c_coa, f_coa, c_cob, f_cob_in, box):
    """Calculate the closest distance assuming fractional, in-cell coords."""
    f_cob = f_cob_in[:]
    fdx = f_coa[0] - f_cob[0]
    if fdx < -0.5:
        f_cob[0] -= 1
    elif fdx > 0.5:
        f_cob[0] += 1
    fdy = f_coa[1] - f_cob[1]
    if fdy < -0.5:
        f_cob[1] -= 1
    elif fdy > 0.5:
        f_cob[1] += 1
    fdz = f_coa[2] - f_cob[2]
    if fdz < -0.5:
        f_cob[2] -= 1
    elif fdz > 0.5:
        f_cob[2] += 1
    if f_cob == f_cob_in:
        # if nothing has changed, use initial values
        return vecdist3(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return vecdist3(c_coa, new_b)


def vecdist3(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    vec = [coord2[0] - coord1[0],
           coord2[1] - coord1[1],
           coord2[2] - coord1[2]]

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5


def remove_files(files, directory='.'):
    """
    Delete any of the files if they exist using standard globbing,
    remove empty directories, or ignore silently if not found.

    """
    del_list = []
    for file_name in files:
        del_list.extend(glob.glob(path.join(directory, file_name)))
    for del_name in del_list:
        debug("deleting %s" % del_name)
        try:
            os.remove(del_name)
        except OSError:
            try:
                os.rmdir(del_name)
            except OSError:
                debug("Directory not empty: %s?" % del_name)


def compress_files(files, directory='.'):
    """Gzip any big files to keep."""
    zip_list = []
    for file_name in files:
        zip_list.extend(glob.glob(path.join(directory, file_name)))
    for zip_name in zip_list:
        debug("compressing %s" % zip_name)
        gzip_command = ['gzip', '-f', zip_name]
        subprocess.call(gzip_command)
        #TODO(tdaff): internal gzip


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


def name_from_types(sites, guest):
    """Generate a string that gives the atoms represented in sites."""
    stypes = []
    for site_ident in sites:
        if site_ident is 0:
            stypes.append('COM')
        else:
            stypes.append(guest.atoms[site_ident-1].element)
    stypes = unique(stypes)
    if len(stypes) > 1:
        site_name = "-".join(stypes)
    else:
        site_name = stypes[0]
    return site_name


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
    job_name = main_options.get('job_name')
    niss_name = "%s.niss" % job_name
    tar_name = "%s.tar" % job_name
    if path.exists(niss_name):
        info("Existing simulation found: %s; loading..." % niss_name)
        load_niss = open(niss_name, 'rb')
        my_simulation = pickle.load(load_niss)
        load_niss.close()
        my_simulation.re_init(main_options)
    else:
        if not main_options.getbool('quiet'):
            welcome()
        info("Starting a new simulation...")
        my_simulation = PyNiss(main_options)

    # Should we extract a previous simulation?
    if main_options.getbool('tar_extract_before') and path.exists(tar_name):
        info("Extracting files from %s" % tar_name)
        faps_tar = tarfile.open(tar_name, 'r')
        faps_tar.extractall()
        faps_tar.close()
        debug("Deleting extracted file %s" % tar_name)
        os.remove(tar_name)

    # run requested jobs
    my_simulation.job_dispatcher()

    # Should we bundle the files after use?
    if main_options.getbool('tar_after'):
        info("Bundling calculation into %s" % tar_name)
        # We don't overwrite an existing tar
        faps_tar = tarfile.open(tar_name, 'a')
        # Similar file names with _underscores_ might match each other so
        # cannot just glob for them
        for suffix in FOLDER_SUFFIXES:
            directory = "faps_%s_%s/" % (job_name, suffix)
            if path.isdir(directory):
               debug("Adding directory %s" % directory)
               faps_tar.add(directory)
               shutil.rmtree(directory)
        faps_tar.close()

    my_simulation.dump_state()
    terminate(0)


if __name__ == '__main__':
    main()
