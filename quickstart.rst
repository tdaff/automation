===================
Quickstart examples
===================

Several example calculations are provided to familiarise the user with the most
frequently used features. These assume that you are working on a system where
faps has been installed and set up. If this is not the case, please consult the
:ref:`installation guide <installation>`.

----------------------
Your first calculation
----------------------

Download :download:`this example structure file <_static/MIL-47.cif>` here and
save it where you will be running faps as :file:`MIL-47.cif`.

Faps is designed to have sensible defaults for most options, so it is possible
to run a simulation without any further input. The usual way to run faps is
using a ``JOBNAME``, in this case this is the name of the cif without ``.cif``.
All the output files and interim calcualtions will use this jobname as the
basename. You can fap the structure by running:

  :command:`faps MIL-47`

.. code-block:: text

  >> Starting faps version 0.999-r367
          .;StX;8@%%8%...
       .XX;::88SXXXX@.8t8.              S@@
      .%t888X8@888@88SXXt@;.            S .X.:
     ;%8@@XXSS%8888@888@@St..             %8XX8S
    .;XXSXX@@X%@8888XS@@88@:8    %X@.%8 :@X:.  X:                    S8:@S8
   .8S8XX88S%SS@8:XS;XXSS88X8. :@ %8: 8t:X@8 . S           .8@8;   S8: t: :.
   @S88@X88@@%X@ .X888888t%:.   :.       :@:@@         .S8S:S8;.8 S.t%
  .%;8SXX@888XSX%8S XX  8:S:;            8;t.         SX S.S  :: :% S.S88S.
  .:@X%@8888888t88t 8@8tX@ @             %;  ::@:     :  S   t  8 .; X8@; ::S
   ;;X%%8@8:X. 88XX@888@8@8%       @.X%8X   .8. tt     @ @ 88;@          ;8 t
   .8Xt8 .  S8.888X%%%S8;%;SS   SXS:;8 @:   .  % 8@    .:.8.: .        t;8:%
    :;Xt%.8;:;X8X8@XSX8%  t@    t 8:        @ 8 ..;%    t 8;       X;;8.@8S
   ..:@XX8 : tS8888888S    :    X ;  8t 8  . @@;:8X.@;  8 S        :S8;8
      :S:tt8: tX888888tX% %:     % :X.XX@ :8 :.;S %  ;Xt8:%
       ;@;%88% ; .:;;.t@8Xt      t  S:    @ @       XSX; .    88:
     t@.%@XX88@ 88@@XS    t      @        @:@               ;  .S
   ; @SS%88 :88t 8%    8;8.      % :%      .               88St       :.
  ;8@S%;;;;:;t88; @88X;           .t:                     8S  S ;@8S: .;
  8@Stt;::::::.tS@ S t                                    8: ;@:;X8Xt.
  8X%t;;:::...::;S8X                                       t;%.;tS
  XS%t;;:::.....:;t8                                          S% ;
                                                                  faps 0.999-r367
  >> Starting a new simulation...
  >> Reading in structure
  >> Reading positions from cif file: MIL-47.cif
  >> Writing state file, MIL-47.niss.
  >> Skipping force field optimisation
  >> Running a vasp calculation
  >> Running on 16 nodes
  >> Optimizing hydrogen positions
  >> Dispersion correction will be used
  >> Running vasp job in queue. Jobid: 36961790
  36961791.wooki.cluster
  >> Writing state file, MIL-47.niss.
  >> Faps terminated normally


If you are running on a cluster, you will see the Faps logo, some output
information and the program will end. The job should now be running through your
queueing system and subsequent sections will automatically run through the
queue. This will generate, in sequence, a DFT hydrogen-optimized structure,
REPEAT charges, and the GCMC uptake of |CO2| at 298 K and 0.15 bar pressure in a
:file:`MIL-47-CO2.csv` file. You can check the status of the job by running
:command:`faps MIL-47` again.

.. code-block:: text

  >> Starting faps version 0.999-r367
  >> Existing simulation found: MIL-47.niss; loading...
  >> Using new options.
  >> Previous system state reported from .niss file (running jobs may have
     already finished):
  >>  * State of esp: Not run
  >>  * State of charges: Not run
  >>  * State of dft: Running, jobid: 36961790
  >>  * State of GCMC: Not run
  >>  * State of init: Processed
  >>  * State of ff_opt: Skipped
  >>  * State of properties: Not run
  >> DFT still in progress
  >> Faps terminated normally

The DFT calculation might take a while so you can have a cup of tea and come
back later to check your results.

All faps output is recorded in the :file:`MIL-47.flog` file, which will also
show any warnings or erros that occure. The file :file:`MIL-47.niss` stores the
processed outputs of all the calculations and the current state of the system;
this file is required to continue a calculation so should be preserved, it can
also be copied to duplicate the structure in different simulation conditions.


---------------------------
Repeating calcualtion parts
---------------------------

Since the `.niss` file stores the state of the system calcaultions that have
already been done will not be repeated. If, however, you change some options you
may need to redo parts of your simulation. The easiest way to do this is to
specify the parts on the commandline before the jobname.

To repeat the gcmc section of the previous simulation you would type:

  :command:`faps gcmc MIL-47`

You can also specify multiple parts:

  :command:`faps dft charges gcmc MIL-47`


--------------------
Isotherm calculation
--------------------

Faps comes with sensible defaults for everything, but offers a lot of
customisability. One way to customise your calculations is with a configuration
file in the directory where your structure is stored, the `.fap` file. The faps
file follows standard ini format with ``option = choice`` syntax. We could
create a :file:`MIL-47.fap` file to calculate a complete isotherm. Simply
specify a number of pressures and temperatures and the uptake will be calculated
for every combination.

.. code-block:: ini

  # MIL-47.fap
  mc_pressure = 0.01 0.1 0.2 0.4 0.8 1.2
  mc_temperature = 263 273


---------------
Uncharged guest
---------------

Sometimes it is more convenient to have the same settings for several jobs, in
this case we can create a centralised jobfile. These need to be stored in your
home directory in a directory called `.faps`. This directory is also used for
other settings and should be used to store descriptive fap files for all your
job types. In this example, the single site methane model does not have any
charged sites, so the charge calculation can be skipped altogether (any dft
optimisation is also skipped). Charges are automatically initialised to 0.


.. code-block:: ini

   # ~/.faps/methane.fap
   # Skip the dft and charges as we don't need
   # them for methane.
   no_dft = true
   no_charges = true
   guests = CH4-TraPPE

To make use of our new jobfile, we run the command:

  :command:`faps -j methane MIL-47`


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
