==================
Customising guests
==================

Faps is distributed with a number of standard guests available, the
``guests.lib`` file in the faps directory has the standard guests, but
it is possible to add custom guests for your simulation. Create a
``guests.lib`` in the job working directory to add your own. The file
follows standard ``.ini`` syntax. The following example is explained
below.


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

The name in square brackets ``[CO2-GS]`` is the identifier to use in the
``guests`` option in the ``.fap`` file. Try to be short but descriptive,
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
