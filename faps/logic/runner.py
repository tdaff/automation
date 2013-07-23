"""
Runner provides the code that runs the traditional faps simulation.

Processes that should take place include:
  * Initialise the settings and user input
  * Derive the dependency graph for the simulation
  * Decide what state the simulation is in
  * Initialise the required calculators
  * Write out the state

"""

from faps.logic.resolver import resolve
from faps.settings import config, path
from faps.managers import JobManager
from faps.io import StructureReader
from faps.backend.file import FileBackend

class Runner(object):
    """
    A complete faps simulation manager. Pulls in all modules and is
    sufficient to do any simulation. The faps script just needs to instance
    this and execute the run method.

    """
    def __init__(self):
        """
        (Re-)initialize the simulation and draw the dependency graph.

        """
        self.simulation = FileBackend.load(config.job_name)
        self.depends = resolve(self.simulation)

    def run(self):
        pass


class PyNiss(object):
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
            if hasattr(guest, 'c_v') and guest.c_v:
                #TODO(tdaff): Make standard in 2.0
                # makes sure that c_v is there and not empty
                cv_header = "C_v,stdev,"
            else:
                cv_header = ""
            # Generate headers separately
            csv = ["#T/K,p/bar,molc/uc,mmol/g,stdev,",
                   "v/v,stdev,wt%,stdev,hoa/kcal/mol,stdev,",
                   guest_excess, he_excess, cv_header,
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
                if cv_header:
                    csv.append("%f,%f," % (guest.c_v[tp_point]))
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
        typed_atoms = self.options.getbool('egulp_typed_atoms')

        mkdirs(qeq_dir)
        os.chdir(qeq_dir)
        debug("Running in %s" % qeq_dir)

        filetemp = open('%s.geo' % job_name, 'w')
        filetemp.writelines(self.structure.to_egulp(typed_atoms))
        filetemp.close()

        # EGULP defaults to GULP parameters if not specified
        # TODO(tdaff): remove old terminaology
        try:
            egulp_parameters = self.options.gettuple('egulp_parameters')
            warning("egulp_parameters is deprecated use qeq_parameters instead")
        except AttributeError:
            egulp_parameters = self.options.gettuple('qeq_parameters')

        if not egulp_parameters:
            # parameters are mandatory in new egulp
            egulp_parameters = ('H', QEQ_PARAMS['H'][0], QEQ_PARAMS['H'][1])
        else:
            info("Custom EGULP parameters selected")

        filetemp = open('%s.param' % job_name, 'w')
        filetemp.writelines(mk_egulp_params(egulp_parameters))
        filetemp.close()

        filetemp = open('%s.ini' % job_name, 'w')
        filetemp.writelines(mk_egulp_ini(self.options))
        filetemp.close()

        egulp_args = ['%s.geo' % job_name,
                      '%s.param' % job_name,
                      '%s.ini' % job_name]

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
        zeo_process = subprocess.Popen(zeo_command, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        zeo_process.wait()

        zeo_stderr = " ".join(x.strip() for x in zeo_process.stderr.readlines())
        print zeo_stderr
        if "Voronoi volume check failed" in zeo_stderr:
            warning("Structure is likely bad; zeo++ is unable to complete")
            warning(zeo_stderr)
            self.structure.bad_structure = True

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
        hydrophilic_area = 0.0
        # gromacs default of 0.2 seems very constrained
        hydrophilic_threshold = 0.3
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

            # Fraction of the accessible surface area for sphere to real area
            if abs(atom.charge) > hydrophilic_threshold:
                hydrophilic_area += (surface_area*ncount)/nsamples
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
        info("Hydrophilic area (A^2) and fraction (probe: %f): %f, %f" %
             (rprobe, hydrophilic_area, hydrophilic_area/total_area))
        return total_area



def welcome():
    """Print any important messages."""
    print(LOGO)
    print(("faps %s" % __version__).rjust(79))


def main():
    """Do a standalone calculation when run as a script."""
    main_options = Options()
    info("Starting faps version %s" % __version__)
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
