"""

Binding site location algorithm
===============================

General procedure is to:

 * use the most central point in the molecule to first place the guest
 * align the guest with nearby probability peaks
 * report the outcome as site + orientation sets

"""

import operator
from logging import debug

from numpy import dot
from numpy.linalg import norm


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


def vecdist3(coord1, coord2):
    """Calculate vector between two 3d points."""
    vec = [coord2[0] - coord1[0],
           coord2[1] - coord1[1],
           coord2[2] - coord1[2]]

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5


def calculate_binding_sites(guest, tp_point, cell):
    """
    Return the sets of points in the maxima that identify binding sites.
    Each set minimally contains a tuple (index of first atom, position of
    second atom), plus tuples of the index and vectors to a second and third
    atom.

    Parameters
    ----------

    guest : Guest
        The guest instance containing the TABASCO maxima.
    tp_point : tuple
        The identifier for the state point to calculate the binding sites for.
    cell : Cell
        The unit cell of the system.

    Returns
    -------

    binding_sites : list
        A list of lists, each containing:
        (index_0, position_0, magnitude_0),
        (index_1, vector_0_1), # [optional]
        (index_2, vector_0_2)  # [optional]

    """

    guest_locations = guest.guest_locations[tp_point]

    # Generate fractional coordinates for all guest sites.
    for pdist in guest_locations:
        guest_locations[pdist] = [
            (loc[0], loc[1], dot(cell.inverse, loc[0]).tolist())
            for loc in guest_locations[pdist]]

    # Work out which atoms to use
    # Only use atom positions for site location
    guest_atom_distances = []

    # Some guests previously only specified the COM probability
    # This is not supported, but fix by saying that COM is the
    # only atom
    if len(guest.atoms) == 1 and guest.probability == [(0, )]:
        guest.probability = [(1, )]
        guest_locations[(1, )] = guest_locations[(0,)]

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
        debug("One atom")
        distance_0_1 = None
        distance_0_2 = None
        distance_1_2 = None
        linear_guest = True
        reversible_guest = True
    elif len(guest_atom_distances) == 2:
        # Two points to align to, use both
        debug("Two atoms")
        distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                                guest.atoms[guest_atom_distances[1][1]].pos)
        distance_0_2 = None
        distance_1_2 = None
        linear_guest = True
        reversible_guest = guest.is_reversible(guest_atom_distances[0][1],
                                               guest_atom_distances[1][1])
    else:
        debug("More than two atoms, take the closest three")
        distance_0_1 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                                guest.atoms[guest_atom_distances[1][1]].pos)
        distance_0_2 = vecdist3(guest.atoms[guest_atom_distances[0][1]].pos,
                                guest.atoms[guest_atom_distances[2][1]].pos)
        distance_1_2 = vecdist3(guest.atoms[guest_atom_distances[1][1]].pos,
                                guest.atoms[guest_atom_distances[2][1]].pos)
        linear_guest = guest.is_linear(guest_atom_distances[0][1],
                                       guest_atom_distances[1][1],
                                       guest_atom_distances[2][1])
        reversible_guest = guest.is_reversible(guest_atom_distances[0][1],
                                               guest_atom_distances[1][1],
                                               guest_atom_distances[2][1])

    # Eugene's tolerance was 0.2
    # Using 0.3 fits the width of the neighbourhood
    overlap_tol = 0.30

    binding_sites = []

    # Atom closest to the COM
    origin_key = guest_atom_distances[0][2]

    # keep a list of sets of atoms that are already included so we don't
    # use them twice
    reversible_sets = []

    for origin_atom in sorted(guest_locations[origin_key],
                              key=operator.itemgetter(1), reverse=True):
        if distance_0_1 is None:
            # add the un-oriented guest
            binding_sites.append([(guest_atom_distances[0][1], origin_atom[0],
                                   origin_atom[1])])
            continue

        # We have more than one atom to align
        align_key = guest_atom_distances[1][2]
        align_closest = (999.9, None)
        align_found = False
        for align_atom in sorted(guest_locations[align_key],
                                 key=operator.itemgetter(1), reverse=True):
            if align_atom == origin_atom:
                continue
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
                if distance_0_2 is None:
                    # found to points to align a two atom guest. Use them!
                    if reversible_guest:
                        this_set = sorted([origin_atom, align_atom])
                        if this_set in reversible_sets:
                            # already have this the other way round
                            # (from the lower occupancy pair)
                            continue
                        else:
                            # first occurrance of this set, use it
                            binding_sites.append([
                                (guest_atom_distances[0][1], origin_atom[0],
                                 origin_atom[1]),
                                (guest_atom_distances[1][1], vector_0_1)])
                            # make sure it's not used again
                            reversible_sets.append(this_set)
                            continue

                    else:
                        # add the guest in this orientation, it's unique
                        binding_sites.append([
                            (guest_atom_distances[0][1], origin_atom[0],
                             origin_atom[1]),
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
                        # Found three points to align to! Binding site here
                        found_site = True
                        # Check we are not duplicating guests that have symmetry
                        if reversible_guest:
                            this_set = sorted([origin_atom, align_atom,
                                               orient_atom])
                            if this_set in reversible_sets:
                                # We already have this one, don't add again
                                continue
                            elif linear_guest:
                                # just add the two sites
                                binding_sites.append([
                                    (guest_atom_distances[0][1], origin_atom[0],
                                     origin_atom[1]),
                                    (guest_atom_distances[1][1], vector_0_1)])
                                # these atom are now taken
                                reversible_sets.append(this_set)
                            else:
                                binding_sites.append([
                                    (guest_atom_distances[0][1], origin_atom[0],
                                     origin_atom[1]),
                                    (guest_atom_distances[1][1], vector_0_1),
                                    (guest_atom_distances[1][1], vector_0_2)])
                                found_site = True
                        elif linear_guest:
                            # just add the two sites
                            binding_sites.append([
                                (guest_atom_distances[0][1], origin_atom[0],
                                 origin_atom[1]),
                                (guest_atom_distances[1][1], vector_0_1)])
                        else:
                            # Add all three sites
                            binding_sites.append([
                                (guest_atom_distances[0][1], origin_atom[0],
                                 origin_atom[1]),
                                (guest_atom_distances[1][1], vector_0_1),
                                (guest_atom_distances[1][1], vector_0_2)])

                if not found_site:
                    if linear_guest:
                        # don't care about reversible guests since no third
                        # site is found, but we can still make the guest with
                        # two sites alone
                        # just add the two sites
                        binding_sites.append([
                            (guest_atom_distances[0][1], origin_atom[0],
                             origin_atom[1]),
                            (guest_atom_distances[1][1], vector_0_1)])
                    #TODO(tdaff): nothing within overlap, use closest (3 sites)
        else:
            if distance_0_2 is None and align_closest[0] > distance_0_1:
                # Very isolated atom, not within 2 distances of any others
                # treat as isolated point atom and still make a guest
                binding_sites.append([(guest_atom_distances[0][1],
                                       origin_atom[0], origin_atom[1])])
            #TODO(tdaff): nothing within overlap, use closest

    return binding_sites

