====================
Example calculations
====================

For all the exmaple calculations it is assumed than :ref:`site.ini <site-ini>`
has been correcty set up for the default paths but all other options are
defaults.


---------------------------
Single point CO2 adsorption
---------------------------

Assuming that the structure file ``structure.pdb`` running the command::

   faps structure

Will generate, in sequence, a DFT hydrogen-optimized structure, REPEAT charges,
and the GCMC uptake of |CO2| at 1 bar pressure in a ``structure-CO2.csv`` file.


--------------------
Isotherm calculation
--------------------



.. code-block:: ini

   # structure.fap
   mc_pressure = 0.01 0.1 0.2 0.4 0.8 1.2
   mc_temperature = 263 273






.. |H2O| replace:: H\ :sub:`2`\ O

.. |CO2| replace:: CO\ :sub:`2`
