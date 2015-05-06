============
Running faps
============

Each job is uniquely identified by the ``$JOBNAME``. To run faps you will
usually just issue the command ``faps $JOBNAME``. A structure file is required,
and a configuration file is optional. By default, all calculation step will run
in sequence and ``$JOBNAME-$GUEST.csv`` files with the adsoprtion data are
produced. A directory ``faps_$JOBNAME_properties`` will contain output from
extra structure related analysis codes.

.. _structure-files:

---------------
Structure files
---------------

Structures can be read from ``CIF``, ``PDB``, ``PQR``, ``VASP``, ``CSSR``
``SQL`` or ``xyz`` files. One file from this list is **required**. By default
the code will look for ``$JOBNAME.pdb``, other file formats may be set in the
config files or on the commandline with the ``initial_structure_format`` option.

The charges can be specified with the structure for the ``CIF``, ``PQR``,
``CSSR`` and ``xyz`` types. Details are given below.

.. object:: $JOBNAME.cif

  The Crystallogrphic Information File (CIF) is the preferred initial strcuture
  format and the default in faps. CIF is a very complex self describing format
  but the parser has been thoroughly tested and is capable of extracting atom
  positions, applying symmetry operations with a check made to remove duplicated
  atoms. The cell must be given and atom positions specified in fractional
  cooridinates.

  The cif parser implements some non-standard fields;
  `_atom_type_partial_charge` will be read in as the atomic charge.
  `_atom_site_description` is used to determine the force-field atom type.

  Bonding is also read in if specified and expanded according to the site
  labels. These are required for forec field optimisations.

  .. warning::

    Be careful with complex structures, symmetry rules and the duplicate atom
    check might not be perfect.


.. object:: $JOBNAME.pdb

   The structure from any standard Protein Data Bank (PDB) structure file,
   ``$JOBNAME.pdb``, can be read by faps. These are fixed column formats and
   faps may not work if your pdb file is broken. The cell is taken from the
   ``cryst1`` line. The position, in columns 31 to 54, and the atom type, at
   column 77, are required for ``ATOM`` or ``HETATM`` records.

   If ``pqr`` input is requested then the extension ``.pqr`` may be used and
   initial charges are also read in from the occupancy column in position 54.
   These will be replaced if a charge calulation is carried out.

   .. warning::

      A .pdb file generated from a cif in Materials Studio can cause problems!
      In the lattice settings ensure that the cristal is re-orineted with `'A
      along X and B in XY'`, and in the peridicity make the crtystal `P1`
      symmetry.


.. object:: {$JOBNAME.{poscar|contcar}|{POSCAR|CONTCAR}}

   Structures that are compatible with VASP 5 can be read by faps. The code
   will read a standard POSCAR type file provided the line with atom types is
   given (i.e. files written by VASP 4 with only the atom counts will not
   work).


.. object:: $JOBNAME.xyz

   Initial coordinates may be given to faps in plain xyz format. The file
   should be structured as one line containing only the number of atoms, a
   comment line (ignored), followed by ``atom_type x.x y.y z.z`` cartesian
   coordinates for each atom. As faps requires periodicity the structure is
   placed in the box defined by the ``default_cell`` option.

   To specify charges in the ``.xyz`` put them in a fifth field on the atom
   line, ``atom_type x.x y.y z.z cha.rge``. Values that cannot be parsed as
   a number will fail silently and have an initial charge of 0.0.

------------
Config files
------------

Configuration files are optional. Faps will try to use sensible defaults for
everything, most installations will require customisation of the ``site.ini``.

.. object:: $JOBNAME.fap

   The per-job settings. This file should be placed in the running directory
   with the same basename as the structure file. This is not a conventional
   *'input file'*, as it is only required if non-defalut options are needed.
   Options are set in standard ``.ini`` format, as described in :ref:`config
   files <config-files>`.


.. object:: ~/.faps/$JOB_TYPE.fap

   Each user may have a directory with standard job types that are called with
   the --job-type commandline option. As with the ``.fap`` file, options are
   set in standard ``.ini`` format, as described in :ref:`config files
   <config-files>`. These options will override the defaults but be overridden
   by a per-job ``.fap`` file.


.. object:: site.ini

   .. _site-ini:

   The ``site.ini`` is located in the code directory and can be used to
   override any of the options set in the ``default.ini`` that is found in the
   same directory. Usually this file will be used to set all configuration for
   a particular machine (e.g. binary or pseudopotential locations), or a set of
   calaultions (desired state points for all high throughput structures).

.. _library-files:

-------------
Library files
-------------

.. object:: guests.lib

   Predefined guests are stored here. The library file is in standard ``.ini``
   format. Faps will search for guests in the working directory first, then in
   the ``~/.faps/`` directory and finally the standard guests library
   distributed with the code. The ``guests.lib`` provided with the code can be
   used as a template for new guests in a :ref:`custom guests.lib
   <custom-guests>` but do not modify this file directly as it will be
   overwritten on updates.
