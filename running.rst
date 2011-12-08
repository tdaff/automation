============
Running faps
============

Each job is uniquely identified by the ``$JOBNAME``. To run faps you will
usually just issue the command ``faps $JOBNAME``. A structure file is required,
and a configuration file is optional. By default, all calculation step will run
in sequence and ``$JOBNAME-$GUEST.csv`` files with the adsoprtion data are
produced.

.. _structure-files:

---------------
Structure files
---------------

Structures can be read from ``CIF``, ``PDB``, ``VASP`` or ``xyz`` files. One
file from this list is **required**. By default the code will look for
``$JOBNAME.pdb``, other file formats may be set in the config files or on the
commandline.

.. object:: $JOBNAME.pdb

   The structure from any standard Protein Data Bank (PDB) structure file,
   ``$JOBNAME.pdb``, can be read by faps. These are fixed column formats and
   faps may not work if your pdb file is broken. The cell is taken from the
   ``cryst1`` line. The position, in colums 31 to 54, and the atom type, at
   column 77, are required for ``ATOM`` or ``HETATM`` records.

   .. warning::

      A .pdb file generated from a cif in Materials Studio can cause problems!
      In the lattice settings ensure that the cristal is re-orineted with `'A
      along X and B in XY'`, and in the peridicity make the crtystal `P1`
      symmetry.


.. object:: $JOBNAME.cif

   Rudimentary Crystallogrphic Information File (CIF) file processing is
   available in faps; atom positions are extracted, symmetry operations are
   applied and a check is made to remove duplicated atoms. The cell must be
   given and atom positions specified in fractional cooridinates.

   .. warning::

      Be careful with complex structures, symmetry rules and the duplicate atom
      check might not be perfect.


.. object:: {$JOBNAME.{poscar|contcar}|{POSCAR|CONTCAR}}

   Structures that are compatible with VASP 5 can be read by faps. The code
   will read a standard POSCAR type file provided the line with atom types is
   given (i.e. files written by VASP 4 will not work)


.. object:: $JOBNAME.xyz

   Initial coordinates may be given to faps in plain xyz format. The file
   should be structured as one line containing only the number of atoms, a
   comment line (ignored), followed by ``atom_type x.x y.y z.z`` cartesian
   coordinates for each atom. As faps requires periodicity the structure is
   placed in the box defined by the ``default_cell`` option.


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
   format and can be used as a template for new guests in a :ref:`custom
   guests.lib <custom-guests>` in the working directory. Do not modify this
   file directly.
