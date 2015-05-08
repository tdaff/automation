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


-----------------------------
List of all available options
-----------------------------


The following is an automatically generated list of most options. The default
value for each option is given here. For the most up-to-date list, see the
:file:`defaults.ini` file.


.. include:: options.inc


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
