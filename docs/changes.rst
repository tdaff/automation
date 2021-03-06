=========
Changelog
=========

-----------
Development
-----------

----------
1.5 series
----------

  * **New** Faps web interface! The fapweb.py starts a webserver that can be
    used to construct input files and submit jobs.
  * ABSL binding sites are optimied with CG instead of zero K MD.
  * Many ABSL tweaks and fixes.
  * Option to email results at the end of a job.
  * Properties run before the GCMC step.
  * Console output changed slightly to be less intense, but shows more
    information on options in ``-v`` verbose mode.
  * Some fixes for the gromacs UFF implementation.

----------
1.4 series
----------

  * **breaking change:** Atoms in the structure are no longer re-ordered when
    reading in the input files so that atom numbering is the same for input and
    output files. This will break ordering when using ``--import`` with previous
    runs. Since VASP doesn't label the POSCAR, this can lead to long lines in
    files that are fixed transparently, but may be missed by a user.
  * **breaking change:** VASP uses DFT-D3 by default instead of DFT-D2.
  * **breaking change:** new directory naming scheme for pressures below 0.1
    bar, allows for higher resolution at low pressure.
  * GROMACS integrated as a force-field optimisation code (GULP is still
    the default)
  * ``.cif`` files with correct bonding informnation (for example, from
    Materials Studio can be automatically typed for force field optimisation
    (OpenBabel permitting).
  * Option to use the Peng-Robinson ``equation_of_state`` to calculate the
    fugacities of gasses on the fly.
  * Probability plots are turned off by default, must be explicitly turned on
    for processing binding sites.
  * ``Structure`` classes have a ``to_cif()`` method.
  * ``.cif`` files are generated for the final structure in the properties
    directory.
  * Use PLATON to generate a powder XRD pattern.
  * High energy and ``nan`` binding sites are pruned from the output.
  * ABSL runs DL_POLY jobs in serial through a single script.
  * MEPO_QEq parameters now accessible by name in ``qeq_parameters``.
  * ``egulp_parameters`` option removed, use ``qeq_parameters`` instead.
  * Add any custom options to a VASP calculation with ``vasp_custom_incar``.


-----------------
previous versions
-----------------

  * See commits on bitbucket: https://bitbucket.org/tdaff/automation/commits/all
