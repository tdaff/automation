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
import sys

from numpy import dot
from numpy.linalg import norm

sys.path.append('..')
import faps
from config import Options
from faps import PyNiss, Structure, Cell, Atom, Guest, Symmetry
from faps import vecdist3, mkdirs


def mk_dl_poly_control(options, dummy=False):
    """CONTROL file for binding site energy calculation."""
    if dummy:
        stats = 1
    else:
        stats = 200
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
        "stats  %i\n" % stats,
        "#traj 1,100,2\n"
        "finish\n"]

    return control


def min_vect(c_coa, f_coa, c_cob, f_cob_in, box):
    """Calculate the closest distance assuming fractional, in-cell coords."""
    f_cob = f_cob_in[:]
    fdx = f_coa[0] - f_cob[0]
    if fdx < -0.5:
        f_cob[0] -= 1
    elif fdx > 0.5:
        f_cob[0] += 1
    fdy = f_coa[1] - f_cob[1]
    if fdy < -0.5:
        f_cob[1] -= 1
    elif fdy > 0.5:
        f_cob[1] += 1
    fdz = f_coa[2] - f_cob[2]
    if fdz < -0.5:
        f_cob[2] -= 1
    elif fdz > 0.5:
        f_cob[2] += 1
    if f_cob == f_cob_in:
        # if nothing has changed, use initial values
        return direction3d(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return direction3d(c_coa, new_b)


def direction3d(source, target):
    """Return the vector connecting two 3d points."""
    return [target[0] - source[0],
            target[1] - source[1],
            target[2] - source[2]]


def calculate_binding_sites(guest, tp_point, cell):
    ### TESTING

    #guest = structure.guests[0]

    #cell = structure.cell

    guest_locations = guest.guest_locations[tp_point]

    # Generate fractional coordinates for all guest sites.
    for pdist in guest_locations:
        guest_locations[pdist] = [
            (loc[0], loc[1], dot(cell.inverse, loc[0]).tolist())
            for loc in guest_locations[pdist]]

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

    if len(guest_atom_distances) == 1:
        # There is only one site, use it
        print("One atom")
        distance_0_1 = None
        distance_0_2 = None
        distance_1_2 = None
        linear_guest = True
    elif len(guest_atom_distances) == 2:
        # Two points to align to, use both
        print("Two atoms")
        distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]],
                                guest.atoms[guest_atom_distances[1][1]])
        distance_0_2 = None
        distance_1_2 = None
        linear_guest = True
    else:
        print("More than two atoms, take the closest three")
        distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                                guest.atoms[guest_atom_distances[1][1]].pos)
        distance_0_2 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                                guest.atoms[guest_atom_distances[2][1]].pos)
        distance_1_2 = vecdist3(guest.atoms[guest_atom_distances[1][1]].pos,
                                guest.atoms[guest_atom_distances[2][1]].pos)
        linear_guest = guest.is_linear(guest_atom_distances[0][1],
                                       guest_atom_distances[1][1],
                                       guest_atom_distances[2][1])

    # Eugene's tolerance was 0.2
    overlap_tol = 0.20

    binding_sites = []

    # Atom closest to the COM
    origin_key = guest_atom_distances[0][2]
    for origin_atom in sorted(guest_locations[origin_key],
                              key=operator.itemgetter(1), reverse=True):
        if distance_0_1 is None:
            # add the un-oriented guest
            binding_sites.append([(guest_atom_distances[0][1], origin_atom[0])])
            continue

        # We have more than one atom to align
        align_key = guest_atom_distances[1][2]
        align_closest = (999.9, None)
        align_found = False
        for align_atom in sorted(guest_locations[align_key],
                                 key=operator.itemgetter(1), reverse=True):
            # don't forget periodic boundaries!
            vector_0_1 = min_vect(origin_atom[0], origin_atom[2],
                                  align_atom[0], align_atom[2], cell.cell)
            separation_0_1 = norm(vector_0_1)
            align_overlap = abs(separation_0_1 - distance_0_1)
            if align_overlap < align_closest[0]:
                align_closest = (align_overlap, align_atom)
            if align_overlap < overlap_tol:
                # We fit the alignment so it's good
                align_found = True
                if linear_guest or distance_0_2 is None:
                    binding_sites.append([
                        (guest_atom_distances[0][1], origin_atom[0]),
                        (guest_atom_distances[1][1], vector_0_1)])
                    continue

                # we can try and fit all three
                orient_key = guest_atom_distances[1][2]
                orient_closest = (999.9, None)
                found_site = False
                for orient_atom in sorted(guest_locations[orient_key],
                                          key=operator.itemgetter(1),
                                          reverse=True):

                    # need all the separations
                    vector_0_2 = min_vect(origin_atom[0], origin_atom[2],
                                          orient_atom[0], orient_atom[2],
                                          cell.cell)
                    separation_0_2 = norm(vector_0_2)
                    vector_1_2 = min_vect(align_atom[0], align_atom[2],
                                          orient_atom[0], orient_atom[2],
                                          cell.cell)
                    separation_1_2 = norm(vector_1_2)

                    overlap_0_2 = abs(separation_0_2 - distance_0_2)
                    overlap_1_2 = abs(separation_1_2 - distance_1_2)

                    if overlap_0_2 + 0.5*overlap_1_2 < orient_closest[0]:
                        orient_closest = (overlap_0_2 + 0.5*overlap_1_2,
                                          orient_atom)
                    if overlap_0_2 < overlap_tol and overlap_1_2 < 2*overlap_tol:
                        binding_sites.append([
                            (guest_atom_distances[0][1], origin_atom[0]),
                            (guest_atom_distances[1][1], vector_0_1),
                            (guest_atom_distances[1][1], vector_0_2)])
                        found_site = True

                if not found_site:
                    #TODO(tdaff): nothing within overlap, use closest
                    pass
        else:
            #TODO(tdaff): nothing within overlap, use closest
            pass


for bs_idx, binding_site in enumerate(binding_sites):

    bs_directory = "%s_bs_%04d" % (guest.ident, bs_idx)

    startdir = os.getcwd()
    mkdirs(bs_directory)
    os.chdir(bs_directory)

    #raise SystemExit

    include_guests = {guest.ident: [guest.aligned_to(*binding_site)]}
    print include_guests

    with open("CONFIG", "w") as config:
        with open("FIELD", "w") as field:
            dlp_files = my_simulation.structure.to_config_field(my_options, include_guests=include_guests)
            config.writelines(dlp_files[0])
            field.writelines(dlp_files[1])

    with open("CONTROL", "w") as control:
        control.writelines(mk_dl_poly_control(my_options))

    os.chdir(startdir)
