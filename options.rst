================
Customising faps
================

Faps aims to provide sensible defaults for most options, however it also aims
to offer many options for customising calculations. Defaults for all options
can be found in the ``defaults.ini``; site-wide options such as the exe names
and directories should be set in the :ref:`site.ini <site-ini>`.

.. _config-files:

------------
Config files
------------

Config files use a simple ini format of ``option = value``, or ``option:
value``. Values will requre a spcific type:

 * ``str``
      Any string value, such as a path or text.

 * ``bool``
      A case-insensitive true/false value; ``1``, ``yes``, ``true``, and ``on``
      will evaluate to ``True``, and ``0``, ``no``, ``false``, and ``off``
      will evaluate to ``False``.

 * ``int``
      An integer value, such as ``1234``, ``0``, ``-12``.

 * ``float``
      A decimal value, such as ``12.5``, ``7``, ``-33.98``.

Many options will take a ``list`` of several values. Items can be separated by
commas and whitespace, and brackets are ignored. In some cases the values
also require a specifc type, this will be noted. The following are all valid
examples, ``float``: ``0.1 0.2 0.3``; ``str``: ``(CO2, CH4)``; ``int``:
``[1,2,1]``.

Lines starting with ``#`` or ``;`` are comments and will be ignored.

The following is a list of most options; for the most up-to-date list, see the
``defaults.ini`` file.

.. program:: faps


.. envvar:: charge_method = repeat

   Method for calculating charge. [str] {repeat}

.. envvar:: dft_code = vasp

   Method to use for dft and/or optimization step. [str] {vasp, siesta}

.. envvar:: dispersion = True

   Turn on empirical dispersion corrections in dft codes that
   support it. [bool]

.. envvar:: esp_resolution = 0.1

   Resolution, in Angstrom, of the esp grid. [float]

.. envvar:: esp_src = vasp

   Source for the ESP. This will usually be the same as the dft_code,
   but not always. [str] {vasp, siesta}

.. envvar:: fastmc_exe = fastmc

   Location of fastmc executable, must be the full path or be in the
   user's $PATH. [str]

.. envvar:: fastmc_ncpu = 1

   Number of CPUs to run fastmc on. Make sure that you use the
   correct fastmc_exe for parallel runs. [int]

.. envvar:: guests = CO2

   Guest(s) to use in GCMC. [str, list] {see guests.lib}

.. envvar:: import = False

   Try to read in data from a previous calculation. [bool]

.. envvar:: initial_structure_format = pdb

   Filetype for input structure file. [str] {pdb, cif}

.. envvar:: interactive = False

   Enable interactive interface. [bool]

.. envvar:: kpoints = (1, 1, 1)

   Kpoint grid size for dft calculations. Ensure that gamma-point only
   exe is not used for >1 kpoint. [(int, int, int)]

.. envvar:: mc_code = fastmc

   Method to use for Monte Carlo calculations. [str] {fastmc}

.. envvar:: mc_cutoff = 12.5

   Potential cutoff to use in GCMC. This will also be used to determine
   the minimum supercell size. [float]

.. envvar:: mc_eq_steps = 1000000

   GCMC equilibration steps. [int]

.. envvar:: mc_history_freq = 1000

   How often to write the fastmc history. [int]

.. envvar:: mc_jobcontrol = False

   Add the 'jobcontrol' directive with fastmc so that GCMC must be
   stopped manually. [bool]

.. envvar:: mc_numguests_freq = 1000

   How often to write the fastmc numguests. [int]

.. envvar:: mc_pressure = 1.0

   GCMC pressure(s) (bar). For multiple pressure points and guests use
   nested lists ((g1p1, g2p1, ...), (g1p2, g2p2, ...), ...) [float, list]

.. envvar:: mc_probability_plot = True

   Turn on probability plots in GCMC. [bool]

.. envvar:: mc_prod_steps = 10000000

   GCMC production steps. [int]

.. envvar:: mc_supercell = (1, 1, 1)

   Supercell to use for GCMC. These values will only be used if the
   individual dimenstions are larger than the supercell calculated from
   the cutoff. [(int, int, int)]

.. envvar:: mc_temperature = 273

   Temperature(s) to use in GCMC (Kelvin). [float, list]

.. envvar:: no_charges = False

   Skip the charge calculation step; Charges will all be zero. [bool]

.. envvar:: no_dft = False

   Skip the dft/optimization step; structure is not optimized and charge
   calculation may fail if it depends on this step. [bool]

.. envvar:: no_gcmc = False

   Skip the gcmc step. [bool]

.. envvar:: no_submit = False

   Do not submit jobs; just create input files. [bool]

.. envvar:: optim_h = True

   Optimize positions of hydrogens in dft/optimization step. [bool]

.. envvar:: optim_all = False

   Optimize all atom positions in dft/optimization step. [bool]

.. envvar:: optim_cell = False

   Optimize cell vectors in dft/optimization step. [bool]

.. envvar:: potcar_dir = vasp_pseudopotentials/

   Location of VASP POTCARs; each element in a folder. [str]

.. envvar:: psf_dir = siesta_pseuodpotentials/

   Location of siesta psf pseudopotentials. [str]

.. envvar:: queue = wooki

   Queuing system to use. [str] {wooki, sharcnet}

.. envvar:: quiet = False

   Silence stdout. This will be ignored here; set on commandline. [bool]

.. envvar:: repeat_exe = repeat.x

   Location of REPEAT executable. [str]

.. envvar:: repeat_ncpu = 1

   Cpus to use for REPEAT calculation. Ensure that repeat_exe points to a
   parallel version if using more than one CPU. [int]

.. envvar:: run_all = True

   Run all the steps without stopping. [bool]

.. envvar:: serial_memory = 2.5

   Maximum memory that can be used for serial calculations (GB). [float]

.. envvar:: siesta_accuracy = med

   General acucracy setting for siesta calcualtions. [str] {low, med, high}

.. envvar:: siesta_exe = siesta

   Location of siesta executable. [str]

.. envvar:: siesta_ncpu = 1

   Number of CPUs to use for siesta. [str]

.. envvar:: siesta_to_cube = siesta2repeat

   Command to convert siesta ESP to .cube file. [str]

.. envvar:: spin = False

   Turn on spin polarization in dft. [bool]

.. envvar:: threaded_codes = repeat

   Codes that run with openmp threads, not mpi. [str, list]

.. envvar:: threaded_memory = 12

   Maximum memory to use for threaded calculations (GB). [float]

.. envvar:: update_opts = True

   Re-read options on restart. [bool]

.. envvar:: vasp_exe = vasp

   Name (location) of vasp executable. [str]

.. envvar:: vasp_ncpu = 8

   Number of cpus to run vasp on. [int]

.. envvar:: vasp_to_cube = vasp2cube 2

   Command to convert LOCPOT to .cube for REPEAT [str]

.. envvar:: verbose = False

   Print debugging information. This will be ignored here; set on commandline.



.. _commandline-options:

-------------------
Commandline options
-------------------

A list of commandline options may be obtained by running ``faps -h`` at any
time. Most options will be set in the :ref:`config files <config-files>` but
all options can be set at runtime with a commandline switch. The most useful
flags are described here.


.. option:: -v, --verbose

   Print additional debugging information to the terminal and the
   ``$JOBNAME.flog`` file.

.. option:: -q, --quiet

   Do not produce any terminal output. All normal output is still logged to the
   ``$JOBNAME.flog`` file.

.. option:: -i, --interactive

   After loading any previous simulation, faps will immediately enter the
   *expert only* interactive mode. This is probably only for debugging and
   fixing calculations. No support for this.

.. option:: -m, --import

   Faps will try to import data from an old or broken simulation and continue
   from there.

.. option:: -n, --no-submit

   Faps will create input files but not submit any jobs.

.. option:: -o, --option

   Allows any option from the :ref:`config file <config-files>` to be specified
   for a single job or step. These will override all other config files.
   Options should be specified as ``key=value`` pairs with no spaces or boolean
   values are set to true when they appear on the commandline. For example
   ``faps -o vasp_ncpu=24 -o spin -o optim_h=false $JOBNAME``, will override
   the number of vasp CPUs, turn on spin and turn off hydrogen optimisation.
