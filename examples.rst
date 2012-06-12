====================
Example calculations
====================

For all the exmaple calculations it is assumed that :ref:`site.ini <site-ini>`
has been correcty set up for the default paths (executables; pseudopotentials)
but all other options are defaults. These customised options can be combined in
the ``structure.fap`` file to tune the calculation.


-----------------------------
Single point |CO2| adsorption
-----------------------------

Assuming that the structure file ``structure.pdb`` running the command::

   faps structure

Will generate, in sequence, a DFT hydrogen-optimized structure, REPEAT charges,
and the GCMC uptake of |CO2| at 1 bar pressure in a ``structure-CO2.csv`` file.


--------------------
Isotherm calculation
--------------------

Specify a number of pressures and temperatures and the uptake will be
calcaulted for every combination.

.. code-block:: ini

   # structure.fap
   mc_pressure = 0.01 0.1 0.2 0.4 0.8 1.2
   mc_temperature = 263 273


---------------
Uncharged guest
---------------

The single site methane model does not have any charged sites, so the charge
calculation can be skipped altogether (any dft optimisation is also skipped).
Charges are automatically initialsed to 0.


.. code-block:: ini

   # structure.fap
   no_dft = true
   no_charges = true
   guests = CH4-TraPPE


--------------
Multiple guest
--------------

Mixtures can be run by specifying multiple guests. This calculation will run
three simulations:

=========== ===========
p(|CO2|)    p(|CH4|)
=========== ===========
0.6         0.2
0.5         0.3
0.4         0.4
=========== ===========

.. code-block:: ini

   # structure.fap
   guests = CO2 CH4-TraPPE
   mc_pressure = (0.6, 0.2), (0.5, 0.3), (0.4, 0.4)


------------------
Siesta calculation
------------------

The default dft package in faps is VASP. Siesta can be used to perform the DFT
geometry optimization and to generate the ESP. Dispersion corrections have not
been inplemented for Siesta.

.. code-block:: ini

   # structure.fap
   dft_code = siesta
   esp_src = siesta
   siesta_accuracy = high
   optim_all = True
   optim_cell = True


--------------------
Charge equilibration
--------------------

For fast charge derivation faps can use the charge equilibration method in
EGULP, which requires no dft and completes within minutes even for 1000+ atom
structures compared to hours or days of CPU time for DFT charges. Charges are
likely to be less accurate and the structure cannot be optimised. If parameters
is blank then the defaults are used.

.. code-block:: ini

   # structure.fap
   no_dft = True
   charge_method = egulp
   egulp_parameters =
       C   5.87730000   5.23176667
       8   9.61510000   7.08292000
      Zn   4.59540000   3.85650000


-------------------------
GULP Charge equilibration
-------------------------

Fast charge equilibration in faps was originally implemented with GULP. This can
still be used, but EGULP is preferred and allows better manipulation of the
parameters. The qeq_fit option can be used to generate a file that will use gulp
to fit the parameters.

.. code-block:: ini

   # structure.fap
   no_dft = True
   charge_method = gulp


-----------------------
Accessible surface maps
-----------------------

By default faps will not calculate the structure properties, such as the
surface area. To skip straight to the surface area step set the following
options:

.. code-block:: ini

   # structure.fap
   # Skip all the calculations
   no_dft = True
   no_charges = True
   no_gcmc = True
   # Parameters for surface calcaultions
   # probes for VdW surface, H2, CO2, and N2
   surface_area_probe = 0.0, 1.42, 1.72, 1.82
   # approximate area per point on the surface
   surface_area_resolution = 0.03
   # write out all the points on the surface (off by default)
   surface_area_save = True
   # Use a spiral point generation algorithm rather than random points (MC)
   surface_area_uniform_sample = True



.. |H2O| replace:: H\ :sub:`2`\ O

.. |CO2| replace:: CO\ :sub:`2`

.. |CH4| replace:: CH\ :sub:`4`
