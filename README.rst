fa\ :sup:`3`\ ps
----------------

fa\ :sup:`3`\ ps - Fully automated adsorption analysis in porous solids.

Automated high throughput structure adsorption analysis and library of
python code for further customising simulations.

This code was written by Tom Daff as a member of the Woo Lab
(http://titan.chem.uottawa.ca) between 2010 and 2015. Since then
the development has continued in the Lab, so please contact someone
there for more information as this version will lack any updates since
then.

Although the key Monte Carlo code must also be obtained from the Woo Lab,
this code contains interfaces to several codes (VASP, SIESTA, GULP,
GROMACS, CIF files, ...), integration with batch queue systems, analysis
methods for porous solids that can easily be used on their own or
transplanted into someone else's code.

The code is unfortunately not a complete showcase of best development
practices (i.e. not ``pip install`` compliant, thoroughly tested,
or fully modularised, ...), but there is a ``faps2.0`` branch that shows the
what the future could have been if time was available in unlimited quantities!

Documentation
=============

Documentation is available online at: https://faps.readthedocs.io

Complete documentation is included in the ``docs`` directory and must be
built with ``sphinx`` before use. The documentation requires the
``sphinxcontrib.bibtex`` and ``numpydoc`` extensions to compile successfully
(``pip install sphinxcontrib-bibtex --user`` and ``pip install numpydoc
--user`` will install these for the current user). Any missing faps
dependencies may also result in missing parts of the API documentation.

Running the command ``make html`` will generate a ``_build`` directory that
contains a file ``index.html`` which can be viewed in any web browser.

If this has been installed on your machine, the documentation may be
available with the commands ``faps-doc`` or ``faps-man``.

Installation
============

Faps requires Python 2.7 or 3.3+. Some older versions were compatible
with Python 2.4, but support for ancient Pythons has now been dropped.

Obtain a copy of the package, put it somewhere and make a link to the
``faps.py`` file in your ``$PATH``. Edit 'site.ini' in the installation
directory to customize all simulations for the machine you will be
running on (an example file is included).

For the web interface, you can also create a link to the file
``web/fapweb.py``.

Some features require third party modules to function correctly. Most
optional:

- ``numpy`` : required for most numerical operations
- ``sqlalchemy`` : (optional) reading and writing to the fapswitch database
  format
- ``scipy`` : (optional) required for binding site location
- ``flask`` : (optional) only needed for the web interface, with either one
  of ``gevent`` or ``tornado``.

Running
=======

Simulations require a structure. Executing the script in the working
directory of the structure will try to automatically run all steps for a
complete analysis.

Disclaimer
==========

This software was designed for our own research purposes and we make no claims
about the suitability for any other purpose.

See the license (BSD 3-clause) for further information.
