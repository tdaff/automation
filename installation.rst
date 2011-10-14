============
Getting faps
============

There are a number of mehthods to install a copy of faps on your system
which depend on your use-case. For a single simulation or a set of
simulations with different parameters it is easiest to use the system
installed version of the code with an individual input file for each.
For a large number of identical jobs it is best to obtain a local copy
of the code and modify the settings for all jobs in that copy.

---------------------
System installed faps
---------------------

Provided that faps is in your ``$PATH`` you can run the system installed
version with the command ``faps``.


-------------
Personal copy
-------------

Simply copy the entire faps directy into your own directory and edit the
``site.ini`` as appropriate.


---------------
Mercurial clone
---------------

If you have read access to the `bitbucket repository
<https://bitbucket.org/tdaff/automation>`_ a current version can be checked
out with the ``hg clone`` command.

To keep your local copy up-to-date issue an ``hg pull`` followed by an
``hg up``.

.. warning::
   Updating your copy will potentially overwrite modified files. Use the
   `site.ini` file for any customisations as this wil not be
   overwritten.


================================
Executables and additional files
================================

-----------
Executables
-----------

On wooki, the set of ``${PROG}submit-faps`` scripts will take care of
deciding on the best executables. For other systems, the configuration
flies must specify the full path or an executable that is in the user's
``$PATH``. These must corespond to the mode you will be running in; e.g.
parallel compiles must be specified if the calculation will request more
than one processor and gamma-point only codes should not be used with
k-point calculations.

----------------
Pseudopotentials
----------------

``vasp`` pseudopotentials must follow the layout in which they are
provided in the vasp download, i.e. ``${ELEMENT}${_SEMICORE}/POTCAR``.
The guidelines in the vasp guide are used to select the ones used in a
calculation. The ``POTCAR`` will normally be deleted after a successful
run.

``siesta`` pseudopotenails must all be in the same folder and named
``${ELEMENT}.psf``.
