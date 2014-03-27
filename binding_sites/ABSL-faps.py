"""

Binding site location algorithm
===============================

General procedure is to:

 * use the most central point in the molecule to first place the guest
 * align the guest with nearby probability peaks
 * report the outcome as site + orientation sets

"""

import operator
import os
import pickle
import sys

from numpy import dot

sys.path.append('..')
import faps
from config import Options
from faps import PyNiss, Structure, Cell, Atom, Guest, Symmetry
from faps import vecdist3, min_dist, mkdirs


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
    for probability in guest.probability:
        # idx for atoms starts from 1; don't lose track!
        if (idx + 1) in probability:
            guest_atom_distances.append((distance, idx, probability))

# from COM outwards
# (distance, atom_idx, prob_key)
guest_atom_distances.sort()
print(guest_atom_distances)

if len(guest_atom_distances) == 1:
    # There is only one site, use it
    print("One atom")
    distance_0_1 = None
    distance_0_2 = None
    distance_1_2 = None
elif len(guest_atom_distances) == 2:
    # Two points to align to, use both
    print("Two atoms")
    distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]],
                            guest.atoms[guest_atom_distances[1][1]])
    distance_0_2 = None
    distance_1_2 = None
else:
    print("More than two atoms, take the closest three")
    distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                            guest.atoms[guest_atom_distances[1][1]].pos)
    distance_0_2 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                            guest.atoms[guest_atom_distances[2][1]].pos)
    distance_1_2 = vecdist3(guest.atoms[guest_atom_distances[1][1]].pos,
                            guest.atoms[guest_atom_distances[2][1]].pos)

    #linear molecules

# Eugene's tolerence was 0.2
overlap_tol = 0.2

binding_sites = []

out_file = open("bs_loc.xyz", 'w')
#out_file = sys.stdout

# Atom closest to the COM
origin_key = guest_atom_distances[0][2]
for origin_atom in sorted(guest_locations[origin_key],
                          key=operator.itemgetter(1), reverse=True):
    if distance_0_1 is None:
        binding_sites.append(origin_atom)  # add the un-oriented guest
        continue

    # We have more than one atom to align
    attached = []
    align_key = guest_atom_distances[1][2]
    align_closest = (999.9, None)
    for align_atom in sorted(guest_locations[align_key],
                             key=operator.itemgetter(1), reverse=True):
        # don't forget periodic boundaries!
        separation_0_1 = min_dist(origin_atom[0], origin_atom[2],
                                  align_atom[0], align_atom[2], cell.cell)
        align_overlap = abs(separation_0_1 - distance_0_1)
        if align_overlap < align_closest[0]:
            align_closest = (align_overlap, align_atom)
        if align_overlap < overlap_tol:
            # We fit the alignment so it's good
            if distance_0_2 is None:
                binding_sites.append((origin_atom, align_atom))
                continue

            # we can try and fit all three
            orient_key = guest_atom_distances[1][2]
            orient_closest = (999.9, None)
            found_site = False
            for orient_atom in sorted(guest_locations[orient_key],
                                      key=operator.itemgetter(1), reverse=True):
                separation_0_2 = min_dist(origin_atom[0], origin_atom[2],
                                          orient_atom[0], orient_atom[2],
                                          cell.cell)
                separation_1_2 = min_dist(align_atom[0], align_atom[2],
                                          orient_atom[0], orient_atom[2],
                                          cell.cell)
                overlap_0_2 = abs(separation_0_2 - distance_0_2)
                overlap_1_2 = abs(separation_1_2 - distance_1_2)

                if overlap_0_2 + 0.5*overlap_1_2 < orient_closest[0]:
                    orient_closest = (overlap_0_2 + 0.5*overlap_1_2, orient_atom)
                if overlap_0_2 < overlap_tol and overlap_1_2 < 2*overlap_tol:
                    binding_sites.append((origin_atom, align_atom, orient_atom))
                    found_site = True

            if not found_site:
                #TODO(tdaff): nothing within overlap, use closest
                pass
    else:
        #TODO(tdaff): nothing within overlap, use closest
        pass

    out_file.write("9\n")
    out_file.write("CO2 guest\n")
    out_file.write("C {0[0]} {0[1]} {0[2]}\n".format(origin_atom[0]))
    for idx in range(8):
        if idx < len(attached):
            out_file.write("O {0[0]} {0[1]} {0[2]} {1}\n".format(attached[idx][0],attached[idx][1]))
        else:
            out_file.write("O 0.0 0.0 0.0\n")

my_options = Options()

for bs_idx, binding_site in enumerate(binding_sites):

    bs_directory = "%s_bs_%04d" % (guest.ident, bs_idx)

    startdir = os.getcwd()
    mkdirs(bs_directory)
    os.chdir(bs_directory)

    include_guests = {guest.ident: [(binding_site[0][0], binding_site[1][0], binding_site[2][0])]}

    with open("CONFIG", "w") as config:
        with open("FIELD", "w") as field:
            dlp_files = my_simulation.structure.to_config_field(my_options, include_guests=include_guests, dummy=True)
            config.writelines(dlp_files[0])
            field.writelines(dlp_files[1])

    with open("CONTROL", "w") as control:
        control.writelines(mk_dl_poly_control(my_options))

    os.chdir(startdir)
