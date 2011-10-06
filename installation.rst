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
