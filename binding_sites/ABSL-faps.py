"""

Binding site location algorithm
===============================

General procedure is to:

 * use the most central point in the molecule to first place the guest
 * align the guest with nearby probability peaks
 * report the outcome as site + orientation sets

"""

import operator
import pickle
import sys

from numpy import dot

sys.path.append('..')
import faps
from config import Options
from faps import PyNiss, Structure, Cell, Atom, Guest, Symmetry
from faps import vecdist3, min_dist


def mk_dl_poly_control(options):
    """CONTROL file for binding site energy calculation."""
    control = [
        "# minimisation\n",
        "zero\n",
        "steps 1000\n",
        "timestep 0.001 ps\n",
        "ensemble nvt hoover 0.1\n",
        "cutoff %f angstrom\n" % options.getfloat('mc_cutoff'),
        "delr 1.0 angstrom\n",
        "ewald precision 1d-6\n",
        "job time 199990 seconds\n",
        "close time 2000 seconds\n",
        "stats  200\n",
        "#traj 1,100,2\n"
        "finish\n"]

    return control


### TESTING

niss_file = open("BaTPPEt.niss", 'rb')
my_simulation = pickle.load(niss_file)
niss_file.close()

guest = my_simulation.structure.guests[0]

cell = my_simulation.structure.cell
guest_locations = guest.guest_locations[(298.0, (0.15, ))]
# Generate fractional coordinates for all guest sites.
for pdist in guest_locations:
    guest_locations[pdist] = [(loc[0], loc[1], dot(cell.inverse, loc[0]).tolist()) for loc in guest_locations[pdist]]


## guest geometry
probability = guest.probability

print(probability)

# Work out which atoms to use
# Only use atom positions for site location
guest_atom_distances = []

for idx, atom in enumerate(guest.atoms):
    distance = vecdist3(guest.com, atom.pos)
    for p_idx, probability in enumerate(guest.probability):
        # idx for atoms starts from 1; don't lose track!
        if (idx + 1) in probability:
            guest_atom_distances.append((distance, idx, p_idx))

# from COM outwards
# (distance, atom_idx, prob_idx)
guest_atom_distances.sort()
print(guest_atom_distances)

if len(guest_atom_distances) == 1:
    # There is only one site, use it
    print("One atom")
elif len(guest_atom_distances) == 2:
    # Two points to align to, use both
    print("Two atoms")
else:
    print("More than two atoms, take the closest three")

    #linear molecules

# Atom descriptions (probability index,

distance_0_1 = 1.149
distance_0_2 = 1.149
distance_1_2 = 2.298

# Eugene's tolerence was 0.2
overlap = 0.2

binding_sites = []

out_file = open("bs_loc.xyz".format(overlap), 'w')
#out_file = sys.stdout
for c_atom in sorted(guest_locations[(1, )], key=operator.itemgetter(1), reverse=True):
    attached = []
    for o_atom in sorted(guest_locations[(2, 3)], key=operator.itemgetter(1), reverse=True):
        # don't forget periodic boundaries!
        #separation_0_1 = vecdist3(c_atom[0], o_atom[0])
        separation_0_1 = min_dist(c_atom[0], c_atom[2], o_atom[0], o_atom[2], cell.cell)
        if abs(separation_0_1 - distance_0_1) < overlap:
            for other_o_atom in sorted(guest_locations[(2, 3)], key=operator.itemgetter(1), reverse=True):
                separation_0_2 = min_dist(c_atom[0], c_atom[2], other_o_atom[0], other_o_atom[2], cell.cell)
                separation_1_2 = min_dist(o_atom[0], o_atom[2], other_o_atom[0], other_o_atom[2], cell.cell)
                #separation_0_2 = vecdist3(c_atom[0], other_o_atom[0])
                #separation_1_2 = vecdist3(o_atom[0], other_o_atom[0])
                if abs(separation_0_2 - distance_0_2) < overlap and abs(separation_1_2 - distance_1_2) < 2*overlap:
                    attached.append(o_atom)
                    attached.append(other_o_atom)
#    print("{}".format(len(attached) + 1))
    out_file.write("9\n")
    out_file.write("CO2 guest\n")
    out_file.write("C {0[0]} {0[1]} {0[2]}\n".format(c_atom[0]))
    for idx in range(8):
        if idx < len(attached):
            out_file.write("O {0[0]} {0[1]} {0[2]} {1}\n".format(attached[idx][0],attached[idx][1]))
        else:
            out_file.write("O 0.0 0.0 0.0\n")

my_options = Options()

include_guests = {"CO2": [((0, 0, 0), (0, 1.149, 0), (0, -1.149, 0))]}

with open("CONFIG", "w") as config:
    with open("FIELD", "w") as field:
        dlp_files = my_simulation.structure.to_config_field(my_options, include_guests=include_guests, dummy=True)
        config.writelines(dlp_files[0])
        field.writelines(dlp_files[1])

with open("CONTROL", "w") as control:
    control.writelines(mk_dl_poly_control(my_options))

