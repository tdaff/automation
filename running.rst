============
Running faps
============

Each job is uniquely identified by the ``$JOBNAME``. To run faps you will
usually just issue the command ``faps $JOBNAME``. A structure file is required,
and a configuration file is optional.

On ``wooki`` the job will stop after each step and can be continued by
re-issuing the command. Jobs on ``orca`` will run to completion by default.

.. _structure-files:

---------------
Structure files
---------------

One file from this list is **required**. By default the code will use a
``$JOBNAME.pdb``, other file formats may be set in the config files or on the
commandline.

.. object:: $JOBNAME.pdb

   Faps will take the structure from any standard Protein Data Bank (PDB)
   structure file, ``$JOBNAME.pdb``.

   .. warning::

      A .pdb file generated from a cif in Materials Studio can cause problems!
      In the lattice settings ensure that the cristal is re-orineted with `'A
      along X and B in XY'`, and in the peridicity make the crtystal `P1`
      symmetry.


.. object:: $JOBNAME.cif

   Faps has rudimentary Crystallogrphic Information File (CIF) file processing
   capability; atom positions are extracted and symmetry operations are
   applied.

   .. warning::

      Faps will use symmetry rules defined in a cif file to generate a crystal,
      however no checking is done for duplicated atoms. Be careful with complex
      structures. These will be fine if only the identity operation is used.


.. _config-files:

------------
Config files
------------

These files are optional. Faps will try to use sensible defaults for
everything, most installations will require customisation of the ``site.ini``.

.. object:: $JOBNAME.fap

   The per-job settings. This file should be placed in the running directory
   with the same basename as the structure file. This is not a conventional
   *'input file'*, as it is only required if non-defalut options are needed.
   Options are set as


.. object:: site.ini

   .. _site-ini:

   The ``site.ini`` is located in the code directory and can be used to
   override any of the options set in ``default.ini`` that is found in the same
   directory. Usually this file will be used to set all configuration for a
   particular machine (e.g. binary or pseudopotential locations), or a set of
   calaultions (desired state points for all high throughput structures).

.. _library-files:

-------------
Library files
-------------

.. object:: guests.lib

   Predefined guests
