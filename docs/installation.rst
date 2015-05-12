.. _installation:

============
Installation
============

There are a number of mehthods to install a copy of faps on your system
which depend on your use-case. For the large majority of users the system
installed faps should be sufficient. Only in the case that you need to modify
the source or you need a different version you can obtain a personal copy.

------------
Requirements
------------

Faps requires Python 2.7 or 3.3+. Some older versions were compatible
with Python 2.4, but support for older versions has now been dropped.
If you need a version compatible with older Pythons, contact a developer.

Some features require third party modules to function correctly. Most
optional:

- ``numpy`` : required for most numerical operations
- ``sqlalchemy`` : (optional) reading and writing to the fapswitch database
  format
- ``scipy`` : (optional) required for binding site location
- ``flask`` : (optional) only needed for the web interface, with either one
  of ``gevent`` or ``tornado``.


---------------------
System installed faps
---------------------

Provided that faps is in your ``$PATH`` you can run the system installed
version with the command :command:`faps`. If it is not in your path, you can
either :command:`/use/path/to/faps.py` or make a symbolic link somewhere in your
path :command:`ln -s /path/to/faps.py ~/bin/faps`.


---------------
Mercurial clone
---------------

The faps code has has been developed using the `mercurial
<http://mercurial.selenic.com/>`_ source control manager. If you have
read access to the `bitbucket repository
<https://bitbucket.org/tdaff/automation>`_ the most current version can
be checked out with the :command:`hg clone` command or downloaded as a
:file:`.tar.gz`. You may also clone an existing repository that you have read
access to.

To keep your local copy up-to-date issue an :command:`hg pull` followed by an
:command:`hg up`. Please see the mercurial manual for more details.

.. warning::
  Updating your copy will potentially overwrite modified files. Use the
  :file:`site.ini` file for any customisations as this wil not be
  overwritten.


========================
Setting up web interface
========================

The script :command:`web/fapweb.py` will run the web interface. To ensure that
the jsmol rendering works you will need to add jsmol files to the static files
of the web interface. The following commands unzip put the jsmol files into
a directory ``web/static/jsmol/``.

.. code-block:: console

   $ cd web/static
   $ curl -O -J -L http://sourceforge.net/projects/jmol/files/latest/download\?source\=files
   $ unzip Jmol-*.zip '*/jsmol.zip
   $ unzip Jmol*/jsmol.zip


===========================================
Setting up executables and additional files
===========================================

-----------
Executables
-----------

In order to run calculations, faps needs access to the required
computational chemistry executables. On wooki, the set of
``${PROG}submit-faps`` scripts will take care of deciding on the best
executables. For other systems, the configuration flies must specify the
full path or an executable that is in the user's ``$PATH``. These must
corespond to the mode you will be running in; e.g. parallel compiles
must be specified if the calculation will request more than one
processor and threaded codes identified by the ``threaded`` option.
Gamma-point only codes should not be used with k-point calculations.

----------------
Pseudopotentials
----------------

Pseudopotentials for ``vasp`` must follow the layout in which they are
provided in the vasp download, i.e. ``${ELEMENT}${_SEMICORE}/POTCAR``.
The vasp guide gives recommendations for the best pseudopotential and is
used to automatically select the ones used in a calculation. The
``POTCAR`` will normally be deleted after a successful run to save
space, add it to the ``compress_files`` list of you wish to keep it.

Pseudopotentials for ``siesta`` must all be in a single directory and
named ``${ELEMENT}.psf``. If the OS allows it, the files will be
symlinked otherwise a copy is made. They are not deleted by default.

---------------------------
Custom job types and guests
---------------------------

Each user may have a dot_faps directory, ``${HOME}/.faps/``, to store
``.fap`` files for common custom job types that are accessed with the
``--job-type`` commandline option, and for a user maintained
``guests.lib`` with a library of common guests.
