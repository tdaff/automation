=======
Options
=======

.. program:: faps


Sensible defaults are provided by faps for most options, however it also aims to
offer many options for customising calculations. Defaults for all options can be
found in the :file:`defaults.ini`, which is located in the source code
directory; site-wide options such as the executable locations, directories and
queueing system and directories should be set in the :ref:`site.ini <site-ini>`.
When using the same options for many jobs the custom job types should be set up
in your :file:`~/faps/$JOB_TYPE.fap` and referenced with the
:option:`-j` commandline option. Settings for an individual job can be
set in the :file:`$JOBNAME.fap` file or on the command line with :option:`-o`.

.. _config-files:

-------------------------
Configuration file format
-------------------------

Config files use a simple ini format of ``option = value``, or ``option:
value``. Values will requre a spcific type:

``str``
  Any string value, such as text or a path.

``bool``
  A case-insensitive true/false value; ``1``, ``yes``, ``true``, and ``on``
  will evaluate to ``True``, and ``0``, ``no``, ``false``, and ``off``
  will evaluate to ``False``.

``int``
  An integer value, such as ``1234``, ``0``, ``-12``.

``float``
  A decimal value, such as ``12.5``, ``7``, ``-33.98``.

Many options will take a ``list`` of several values. Items can be separated by
commas and whitespace, and brackets are ignored. In some cases the values
also require a specifc type, this will be noted. The following are all valid
examples, ``float``: ``0.1 0.2 0.3``; ``str``: ``(CO2, CH4)``; ``int``:
``[1,2,1]``.

Lines starting with ``#`` or ``;`` are comments and will be ignored.

Options should be one per line. If you wish to span multiple lines, subsequent
ones should be indented. Options are evaluated hierarchically and will silently
fall back to defaults in the order (use the :option:`-v` option to show which
options are used):

* Command line arguments or option from :option:`-o`
* Job specific option from :file:`$JOBNAME.fap`
* Job type :file:`~/.faps/JOB_TYPE.fap` specified with :option:`-j`
* Site specific from :file:`site.ini`
* Default from :file:`defaults.ini`



The following is an automatically generated list of most options. The default
value for each option is given here. For the most up-to-date list, see the
:file:`defaults.ini` file.


.. envvar:: charge_method

  Default: repeat

  yes, on, 1 and False, no, off, 0
  Method for calculating charge. [str] {repeat, gulp, egulp}

.. envvar:: dedicated_queue

  Default:

  Use a a specific queue when submitting jobs. Leave blank for none [str]
  e.g. NRAP12345 on sharcnet...

.. envvar:: default_cell

  Default: (10.0, 10.0, 10.0, 90, 90, 90)

  Cell parameters in A to use only when they are not specified in the
  input structure (e.g. for .xyz). Use either (a, b, c, alpha, beta, gamma) or
  ax, ay, az,
  bx, by, bz,
  cx, cy, cz
  notation. [float, list]

.. envvar:: dft_code

  Default: vasp

  Method to use for dft and/or optimization step. [str] {vasp, siesta}

.. envvar:: dispersion

  Default: True

  Turn on empirical dispersion corrections in dft codes that
  support it. [bool]

.. envvar:: egulp_exe

  Default: egulp

  Location of Eugene's QEq code. [str]

.. envvar:: egulp_typed_atoms

  Default: False

  Identify Sulphonic acid and nitro oxygen as separate types in egulp. [bool]

.. envvar:: esp_resolution

  Default: 0.1

  Target resolution of the esp grid (A). [float]

.. envvar:: esp_src

  Default: vasp

  Source for the ESP. This will usually be the same as the dft_code,
  but not always. [str] {vasp, siesta}

.. envvar:: fastmc_compress_files

  Default: \*.cube

  files to keep and compress after a successful fastmc job [str, list]

.. envvar:: fastmc_delete_files

  Default:

  files to delete after a successful fastmc job [str, list]

.. envvar:: fastmc_exe

  Default: fastmc

  Location of fastmc executable, must be the full path or be in the
  user's $PATH. [str]

.. envvar:: fastmc_keep_unfolded_cubes

  Default: False

  Should the original cube files be kept after folding (or deleted)? [bool]

.. envvar:: fastmc_ncpu

  Default: 1

  Number of CPUs to run fastmc on. Make sure that you use the
  correct fastmc_exe for parallel runs. [int]

.. envvar:: ff_opt_code

  Default: gulp

  Method to use to do force field optimisations. [str]

.. envvar:: find_maxima

  Default: True

  Calculate the location of the guests from the probability cube [bool]

.. envvar:: fold

  Default: True

  Fold probability cube into the unit cell [bool]

.. envvar:: guests

  Default: CO2

  Guest(s) to use in GCMC. [str, list] {see guests.lib}

.. envvar:: gulp_terse

  Default: False

  Reduce the standard output of GULP to a minimum with no movie files [bool]

.. envvar:: gulp_exe

  Default: gulp

  Location of GULP exe. [str]

.. envvar:: import

  Default: False

  Try to read in data from a previous calculation. [bool]

.. envvar:: initial_structure_format

  Default: cif

  Filetype for input structure file. [str] {pdb, cif, vasp, xyz}

.. envvar:: interactive

  Default: False

  Enable interactive interface. [bool]

.. envvar:: kpoints

  Default: (1, 1, 1)

  Kpoint grid size for dft calculations. Ensure that gamma-point only
  exe is not used for >1 kpoint. [(int, int, int)]

.. envvar:: mc_code

  Default: fastmc

  Method to use for Monte Carlo calculations. [str] {fastmc}

.. envvar:: mc_cutoff

  Default: 12.5

  Potential cutoff in A to use in GCMC. This will also be used to
  determine the minimum supercell size. [float]

.. envvar:: mc_eq_steps

  Default: 1000000

  GCMC equilibration steps. Use negative values for cycles. [int]

.. envvar:: mc_history_freq

  Default: 0

  How often to write the fastmc history. [int]

.. envvar:: mc_jobcontrol

  Default: False

  Add the 'jobcontrol' directive with fastmc so that GCMC must be
  stopped manually. [bool]

.. envvar:: mc_numguests_freq

  Default: 1000

  How often to write the fastmc numguests.out file. [int]

.. envvar:: mc_pressure

  Default: 0.15

  GCMC pressure(s) (bar). For multiple pressure points and guests use
  nested lists ((g1p1, g2p1, ...), (g1p2, g2p2, ...), ...), these are
  all run at every temperature to generate isotherms. [float, list]

.. envvar:: mc_probability_plot

  Default: True

  Turn on probability plots in GCMC. [bool]

.. envvar:: mc_probability_plot_spacing

  Default: 0.1

  Distance between grid points (resolution) in A for GCMC
  probability plot. [float]

.. envvar:: mc_prod_steps

  Default: 10000000

  GCMC production steps. Use negative values for cycles. [int]

.. envvar:: mc_state_points

  Default:

  Individual state points to run gcmc simulations; not combined with
  temperature/pressure isotherms. Specify points (bar/Kelvin) as:
  (T1, (g1p1, g2p1, ...)), (T2, (g1p2, g2p2, ...), ... [float, list]

.. envvar:: mc_supercell

  Default: (1, 1, 1)

  Minimum supercell to use for GCMC. These values will only be used
  if the individual dimenstions are larger than the supercell calculated
  from the cutoff. [(int, int, int)]

.. envvar:: mc_temperature

  Default: 298

  Temperature(s) to use in GCMC (Kelvin) combined with pressures to
  collect isotherms. [float, list]

.. envvar:: mc_zero_charges

  Default: False

  Do not use charges in the GCMC even if they have been calculated. [bool]

.. envvar:: no_charges

  Default: False

  Skip the charge calculation step; If charges are not set in the
  input file or by other means they will all be zero. [bool]

.. envvar:: no_dft

  Default: False

  Skip the dft/optimization step; structure is not dft optimized and charge
  calculation may fail if it depends on the esp from this step. [bool]

.. envvar:: no_force_field_opt

  Default: True

  Do not pre-optimise with a force field. Ususally requires topology
  defined in the input file. [bool]

.. envvar:: no_gcmc

  Default: False

  Skip the gcmc step. [bool]

.. envvar:: no_properties

  Default: False

  Skip the property calculations. [bool]

.. envvar:: no_submit

  Default: False

  Do not submit jobs; just create input files. [bool]

.. envvar:: optim_h

  Default: True

  Optimize positions of hydrogens in dft/optimization step. [bool]

.. envvar:: optim_all

  Default: False

  Optimize all atom positions in dft/optimization step. [bool]

.. envvar:: optim_cell

  Default: False

  Optimize cell vectors in dft/optimization step. [bool]

.. envvar:: plain

  Default: False

  Do not colourise output. Ignored here; use commandline. [bool]

.. envvar:: potcar_dir

  Default: vasp_pseudopotentials/

  Location of VASP POTCARs; each element in a folder. [str]

.. envvar:: psf_dir

  Default: siesta_pseuodpotentials/

  Location of siesta psf pseudopotentials. [str]

.. envvar:: qeq_fit

  Default: False

  Fit charge equilibration parameters to calculated charges. [bool]

.. envvar:: qeq_parameters

  Default:

  custom parameter sets for QEq as (atom, electronegativity, 0.5\*hardness)
  [(int/str, float, float), list]

.. envvar:: queue

  Default: wooki

  Queuing system to use. [str] {wooki, sharcnet}

.. envvar:: quiet

  Default: False

  Only output errors. This will be ignored here; set on commandline. [bool]

.. envvar:: repeat_compress_files

  Default: \*.cube

  files to keep and compress after a successful REPEAT job. [str, list]

.. envvar:: repeat_delete_files

  Default: ESP_real_coul.dat fort.30 fort.40 REPEAT_param.inp

  files to delete after a successful REPEAT job. [str, list]

.. envvar:: repeat_exe

  Default: repeat.x

  Location of REPEAT executable. [str]

.. envvar:: repeat_ncpu

  Default: 1

  Cpus to use for REPEAT calculation. Ensure that repeat_exe points to a
  parallel version if using more than one CPU. [int]

.. envvar:: run_all

  Default: True

  Run all the steps without stopping. [bool]

.. envvar:: serial_memory

  Default: 2.5

  Maximum memory that can be used for serial calculations (GB).
  If using serial repeat this can restrict the resolution of the ESP. [float]

.. envvar:: siesta_accuracy

  Default: med

  General acucracy setting for siesta calcualtions. [str] {low, med, high}

.. envvar:: siesta_compress_files

  Default:

  Files to keep and compress after a successful SIESTA job. [str, list]

.. envvar:: siesta_delete_files

  Default: \*.ion \*.xml INPUT_TMP\* \*.DM

  Files to delete after a successful SIESTA job. [str, list]

.. envvar:: siesta_exe

  Default: siesta

  Location of siesta executable. [str]

.. envvar:: siesta_ncpu

  Default: 1

  Number of CPUs to use for siesta. [str]

.. envvar:: siesta_to_cube

  Default: siesta2repeat

  Command to convert siesta ESP to .cube file. [str]

.. envvar:: silent

  Default: False

  Only emit critical messages. Ignored here; use commandline. [bool]

.. envvar:: spin

  Default: False

  Turn on spin polarization in dft. [bool]

.. envvar:: surface_area_probe

  Default:

  Radius of probe for calculating surface areas. A probe of radius 0.0 will
  generate the VdW surface typical values for probe molecules are 1.42 (H2),
  1.72 (CO2) or 1.82 (N2) (A). [float, list]

.. envvar:: surface_area_resolution

  Default: 0.03

  Approximate area per point when subdividing accessible surface
  areas (A^2). [float]

.. envvar:: surface_area_save

  Default: False

  Save the valid points on the surface to a file. [bool]

.. envvar:: surface_area_uniform_sample

  Default: False

  Use points with a uniform spacing? (or do Monte Carlo sampling) [bool]

.. envvar:: symmetry

  Default: False

  Treat symmetrical atoms as equivalent for charges. [bool]

.. envvar:: tar_after

  Default: False

  Bundle all the output into an archive when the job is finished. [bool]

.. envvar:: tar_extract_before

  Default: False

  Extract all the tarred files before starting (useful for import). [bool]

.. envvar:: threaded_codes

  Default: repeat

  Codes that run with openmp threads, not mpi. [str, list]

.. envvar:: threaded_memory

  Default: 12

  Maximum memory to use for threaded calculations (GB). [float]

.. envvar:: update_opts

  Default: True

  Re-read options on restart. [bool]

.. envvar:: vasp_compress_files

  Default: LOCPOT CHGCAR vasprun.xml

  files to keep and compress after a successful VASP job. [str, list]

.. envvar:: vasp_delete_files

  Default: WAVECAR CHG DOSCAR EIGENVAL POTCAR PCDAT IBZKPT XDATCAR KPOINTS

  files to delete after a successful VASP job. [str, list]

.. envvar:: vasp_exe

  Default: vasp

  Name (location) of vasp executable. [str]

.. envvar:: vasp_ncpu

  Default: 8

  Number of cpus to run vasp on. [int]

.. envvar:: vasp_to_cube

  Default: vasp_to_cube

  Command to convert LOCPOT to .cube for REPEAT. [str]

.. envvar:: verbose

  Default: False

  Print debugging information. This will be ignored here; set on commandline.

.. envvar:: zeo++

  Default: True

  Run zeo++ on the structure? [bool]

.. envvar:: zeo++_exe

  Default: network

  Command to run for zeo++. Faps will generate the command lines. [str]

.. envvar:: daemon

  Default: False


  Function switch specific options

  Run fapswitch as a service and wait for line-by-line input.
  See also commandline options. [bool]

.. envvar:: fapswitch_backends

  Default: file

  Backends to store the output structures. [str, list] {file, sqlite}

.. envvar:: fapswitch_connectivity

  Default: openbabel

  Where to get the connectivity information from [str] {openbabel, file}

.. envvar:: fapswitch_custom_strings

  Default:

  Make functionalisations with the set of {.freeform.srings.} and
  [symm@try.strings]. [str, list]

.. envvar:: fapswitch_full_random_count

  Default: 0

  Number of completely randomised structures to make. [int]

.. envvar:: fapswitch_max_different

  Default: 0

  Maximum number of groups that will be used simultaneously. [int]

.. envvar:: fapswitch_port

  Default: 0

  Socket port to run the server mode on; leave as zero to pick random
  available port as two instances cannot share a port. [int]

.. envvar:: fapswitch_replace_all_sites

  Default: False

  Should fapswitch produce all group@site combinations? [bool]

.. envvar:: fapswitch_replace_groups

  Default:

  Only use the specified groups in systematic functionalisations. [list]

.. envvar:: fapswitch_replace_only

  Default:

  Only replace the listed sites in systematic functionalisations. [list]

.. envvar:: fapswitch_site_random_count

  Default: 0

  Number of symmetry based randomised structures to make. [int]

.. envvar:: fapswitch_unfunctionalised_probability

  Default: 0.5

  Probability that a site will have no functionalisation in random switching
  scheme. [float]


.. _commandline-options:

-------------------
Commandline options
-------------------

A list of commandline options may be obtained by running :command:`faps -h` at
any time. Most options will be set in the :ref:`config files <config-files>` but
all options can be set at runtime with a commandline switch. The most useful
flags are described here.


.. option:: -v, --verbose

  Print additional debugging information to the terminal and the
  :file:`$JOBNAME.flog` file.

.. option:: -q, --quiet

  Only output errors and warnings to the terminal. All normal output is still
  logged to the :file:`$JOBNAME.flog` file.

.. option:: -s, --silent

  Do not produce any terminal output except for critical errors. All normal
  output is still logged to the :file:`$JOBNAME.flog` file.

.. option:: -p, --plain

  Do not colorise or wrap the terminal output. Default is to colour the
  information and wrap the text at 80 characters. File output is always plain.

.. option:: -j <$JOB_TYPE>, --job-type=<$JOB_TYPE>

  Use the :file:`~/.faps/$JOB_TYPE.fap` file to set options for the current job.
  This will override defaults but options will still be overridden by
  :file:`$JOBNAME.fap` and options set on the commandline.

.. option:: -m, --import

  Faps will try to import data from an old or broken simulation and continue
  from there.

.. option:: -n, --no-submit

  Faps will create input files but not submit any jobs. As steps may depend on
  each other, calculations may need to run to continue the simulations.

.. option:: -o <CMDOPTS>, --option=<CMDOPTS>

  Allows any option from the :ref:`config file <config-files>` to be specified
  for a single job or step. These will override all other config files. Options
  should be specified as ``key=value`` pairs with no spaces or boolean values
  are set to true when they appear on the commandline. For example
  :command:`faps -o vasp_ncpu=24 -o spin -o optim_h=false $JOBNAME`, will
  override the number of vasp CPUs, turn on spin and turn off hydrogen
  optimisation.

.. option:: -i, --interactive

  After loading any previous simulation, faps will immediately enter the
  *expert only* interactive mode. This is probably only for debugging and
  fixing calculations. No support for this.

.. option:: -d, --daemon

  Fapswitch only: Start the program as a service; see the fapswitch
  documentation for how to use this.
