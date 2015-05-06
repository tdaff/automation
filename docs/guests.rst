======
Guests
======

Faps is distributed with a number of standard guests available in the
:file:`guests.lib` file in the faps directory. It is also possible to add
custom guests for your simulation either in :file:`~/.faps/guests.lib`
that is queried for all simulations or a :file:`$JOBDIR/guests.lib` that
is queried for all jobs running in that directory.

----------------
Available guests
----------------

These guests are available as standard with the identifier with the
:envvar:`guests` option.

============ ===========================================================
Identifier   Guest
============ ===========================================================
CO2          Carbon Dioxide for Adsorption in Zeolites :cite:`Garcia-Sanchez2009`
CO2-TraPPE   TraPPE Carbon Dioxide :cite:`Potoff2001`
CO-MOF       Carbon Monoxide Optimised for CuBTC MIL-47 and IRMOF1 :cite:`Martin-Calvo2012`
CO-UFF       Carbon Monoxide with UFF Parameters and Mulliken Charges
CO-Graphite  Carbon Monoxide Adsorbed on Graphite
CO-Ligand    Carbon Monoxide in Myoglobin
TIP5P        TIP5P Water
TIP5P-Ew     TIP5P-Ew Water for Ewald sums
TIP4P-Ew     TIP4P-Ew Four-Site Water Model for Biomolecular Simulations
CH4-TraPPE   TraPPE Methane; Single Site
CH4-5s       Methane 5 Site
N2           Three Site Nitrogen Based on TraPPE
N2-2site     Two Site, Chargeless Nitrogen
N2-Etters    Nitrogen Based on Etters Potential
H2S          Hydrogen Sulphide for Vapor-Liquid Equilibria
N2O-PE       Nitrous Oxide for Phase Equilibria with Oxygen
N2O-Sol      Nitrous Oxide for Solubility in Fluorinated Liquids
NO2          Nitrogen Dioxide
N2O4         Dinitrogen Tetroxide
NO           Nitric Oxide for Diffusion in Water
O2-TraPPE    Oxygen for Liquid Vapour Phase Equilibria
O2-PE        Oxygen for Phase Equilibria with Nitrous Oxide
SO2          Sulphur Dioxide for Vapor-Liquid Equilibria
H2           Anisotropic Hydrogen Diatom for Condensed Phase
I2           Iodine Diatom with UFF parameters
============ ===========================================================


.. _custom-guests:

------------------
Customising guests
------------------

To use custom guests create a ``guests.lib`` in the job working
directory to add your own. The file follows standard ``.ini`` syntax.
The following example is explained below.


.. code-block:: ini

    [CO2-GS]
    name = Carbon Dioxide (Garcia Sanchez)
    source = http://dx.doi.org/10.1021/jp810871f
    atoms =
        Cx  12.0107  0.65120  0.0    0.0  0.0
        Ox  15.9994 -0.32560  1.149  0.0  0.0
        Ox  15.9994 -0.32560 -1.149  0.0  0.0
    potentials =
        Cx  2.745  0.05948
        Ox  3.017  0.17023
    probability =
        0
        2  3
    probe radius = 1.71


The name in square brackets ``[CO2-GS]`` is the identifier to use in the
``guests`` option in your ``.fap`` file. Try to be short but descriptive,
but this must be unique.

``name`` is the full name and should differentiate it from similar
guests.

Include a reference in ``source`` so that the results are always
traceable.

The ``atoms`` block has the molecule's structure. The block **must be
indented** so that all the atoms are read in. Each line must follow:
``atom_name atom_mass atom_charge atom_x atom_y atom_z``. If
``atom_name`` is the element then standard UFF parameters will be used,
use unique identifiers, like ``Cx``, with the ``potentials`` to use
custom potential parameters.

The ``potentials`` block specifies custom potentials to use when running
simulations with this molecule, in the form ``atom_name atom_sigma
atom_epsilon``. The vdW parameters are standard sigma and epsilon values
for Lorentz-Berthelot mixing. Any number of atom types can be specified
here; If you set parameters for standard element
names, like ``C``, they will be used instead of the standard UFF values.
Unspecified atoms will be set to 0.

Goups of atoms can be chosen for probability plotting in the
``probability`` block. Each line is the atom identifiers for a single
group (1 is the first atom), 0 is the centre of mass. In this case the
second plot is for the oxygen atoms.

.. bibliography:: library.bib
