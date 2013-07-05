"""

Reader for VASP coordinate files (Crystallographic information files).

Parses symmetry and some custom properties.

"""
import shlex
from logging import debug, info, warning, error

from faps.core.structure import Structure
from faps.core.symmetry import Symmetry
from faps.core.cell import Cell
from faps.util.text import strip_blanks, ufloat
from . import StructureReader

class GulpWriter(StructureWriter):
    """
    CIF file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def to_gromacs(self):
        """Procedure:
        generate .gro atom positions
        generate .top lj topology
        generate .itp uff topology
        gererate .mdp parameters
        grompp -f mdp_file -c gro_file -p top_file -o output.tpr
        mdrun -s moffive -o traj-c moffive -nt 1
        """

        ##
        # .gro file
        ##
        # TITLE linw
        # number of atoms
        gro_file = ["%s\n" % self.name, " %i\n" % self.natoms]
        # TODO(tdaff): check residue in toplogy?
        residue = (1, self.name[:4].upper())
        for idx, atom in enumerate(self.atoms):
            # In the specification the atom position is given as
            # %8.3f in nm, but we probably would like more accuracy
            # but make sure there are 5 figures to left of decimal point
            pos_nm = tuple([x/10.0 for x in atom.pos])
            gro_file.append(("%5i%5s" % residue) +
                            ("%5s%5i" % (atom.uff_type, idx + 1)) +
                            ("%15.10f%15.10f%15.10f\n" % pos_nm))

        # cell is also in nm
        # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        box = self.cell.cell/10.0
        # gromacs needs these condtiions. Might need to rotate cell sometimes?
        if not (box[0][1] == box[0][2] == box[1][2] == 0):
            error("Gromacs can't handle this cell orientation")
        gro_file.append("%f %f %f %f %f %f %f %f %f\n" % (
            box[0][0], box[1][1], box[2][2],
            box[0][1], box[0][2],
            box[1][0], box[1][2],
            box[2][0], box[2][1]))

        ##
        # .top file
        ##
        comb_rule = 1 # 1 = LORENTZ_BERTHELOT; 2 = GEOMETRIC
        top_file = [
            "; Topology file for %s\n" % self.name,
            "[ defaults ]\n",
            "; nbfunc      comb-rule      gen-pairs       fudgeLJ    fudgeQQ\n",
            "1             %i              yes             1.0        1.0\n\n" % (comb_rule)]

        # define the forcefield for UFF
        unique_types = set()
        top_file.append("[ atomtypes ]\n")
        top_file.append("; name1 name2   mass     charge  ptype   sigma   epsilon\n")
        for atom in self.atoms:
            if atom.uff_type not in unique_types:
                uff_type = atom.uff_type
                # sigma = x1 * 2^(-1/6)
                # CONVERT TO nm !!
                sigma = 0.1*UFF_FULL[uff_type][2]*(2**(-1.0/6.0))
                # epsilon = D1 in kcal
                epsilon = UFF_FULL[uff_type][3] * KCAL_TO_KJ
                top_file.append("%-6s   %-6s   %9.4f   %9.4f   A %12.7f %12.7f\n" % (uff_type, uff_type, atom.mass, 0.0, sigma, epsilon))
                unique_types.add(uff_type)

        top_file.extend(["\n#include <%s.itp>\n\n" % self.name,
                         "[ system ]\n",
                         "UFF optimisation of %s\n\n" % self.name,
                         "[ molecules ]\n",
                         "%s  1\n" % residue[1]])

        ##
        # .itp file
        ##
        # exclude 3 neighbours
        itp_file = ["[ moleculetype ]\n",
                    "; molname nrexcl\n"
                    "%s 3\n" % residue[1],
                    "[ atoms ]\n",
                    "; nr type  resnr    residue    atom     cgnr    charge       mass \n"]

        # atoms
        for idx, atom in enumerate(self.atoms):
            uff_type = atom.uff_type
            # charge group is different for each atom as gromacs has max of 32 in a group
            itp_file.append("%-6i  %-6s  %i  %-6s  %-6s   %d  %9.4f  %9.4f\n" % (idx+1, uff_type, residue[0], residue[1], uff_type, idx+1, atom.charge, atom.mass))

        # bonds
        itp_file.append(" [ bonds ]\n")
        itp_file.append("; ai aj funct b0 kb\n")

        unique_bonds = {}
        bonding_table = dict([(x, {}) for x in range(self.natoms)])
        #FIXME(tdaff) bonds will have length in v2
        for bond, bondorder in self.bonds.items():
            idx_a = bond[0]
            idx_b = bond[1]
            atom_a = self.atoms[idx_a]
            atom_b = self.atoms[idx_b]
            uff_a = atom_a.uff_type
            uff_b = atom_b.uff_type
            typed_bond = tuple(sorted([uff_a, uff_b]) + [bondorder])

            bonding_table[idx_a][idx_b] = typed_bond
            bonding_table[idx_b][idx_a] = typed_bond

            # have we already calculated the parameters
            if not typed_bond in unique_bonds:
                ri = UFF_FULL[uff_a][0]
                rj = UFF_FULL[uff_b][0]
                chiI = UFF_FULL[atom_a.uff_type][8]
                chiJ = UFF_FULL[atom_b.uff_type][8]
                rbo = -0.1332*(ri+rj)*log(bondorder)

                ren = ri*rj*(((sqrt(chiI) - sqrt(chiJ))**2)) / (chiI*ri + chiJ*rj)
                r0 = 0.1*(ri + rj + rbo - ren)

                # force constant
                # parameters Z1
                kb = (KCAL_TO_KJ * 664.12 * UFF_FULL[uff_a][5] * UFF_FULL[uff_b][5])/(r0**3)

                unique_bonds[typed_bond] = (r0, kb)

            params = unique_bonds[typed_bond]
            bond_func = 1  # gromacs harmonic
            # add 1 to bond as 1 indexed
            itp_file.append('%-5i %-5i %1i %11.4f %11.4f ; %-5s %-5s %.2f\n' % (bond[0]+1, bond[1]+1, bond_func, params[0], params[1], typed_bond[0], typed_bond[1], bondorder))

        # angles
        itp_file.append("\n [ angles ]\n")

        for idx_a in sorted(bonding_table):
            bonded_atoms = bonding_table[idx_a]
            for l_idx in bonded_atoms:
                for r_idx in bonded_atoms:
                    if l_idx == r_idx:
                        continue
                    # FIXME(tdaff) special cases for 5 or more coord
                    central_atom = self.atoms[idx_a]
                    l_atom = self.atoms[l_idx]
                    r_atom = self.atoms[r_idx]
                    coordination = len(bonded_atoms)
                    # FIXME(tdaff) over/under coordinated atoms?
                    theta0 = UFF_FULL[central_atom.uff_type][1]
                    cosT0 = cos(theta0*DEG2RAD)
                    sinT0 = sin(theta0*DEG2RAD)
                    c2 = 1.0 / (4.0 * sinT0 * sinT0)
                    c1 = -4.0 * c2 * cosT0
                    c0 = c2*(2.0*cosT0*cosT0 + 1.0)
                    zi = UFF_FULL[l_atom.uff_type][5]
                    zk = UFF_FULL[r_atom.uff_type][5]
                    bond_ab = tuple(sorted((l_idx, idx_a)))
                    bond_ab = tuple(sorted([l_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_ab]])
                    bond_bc = tuple(sorted((r_idx, idx_a)))
                    bond_bc = tuple(sorted([r_atom.uff_type, central_atom.uff_type]) + [self.bonds[bond_bc]])
                    rab = unique_bonds[bond_ab][0]
                    rbc = unique_bonds[bond_bc][0]
                    rac = sqrt(rab*rab + rbc*rbc - 2.0 * rab*rbc*cosT0)

                    ka = (644.12 * KCAL_TO_KJ) * (zi * zk / (rac**5.0))
                    ka *= (3.0*rab*rbc*(1.0 - cosT0*cosT0) - rac*rac*cosT0)

                    thetamin = (pi-arccos(c1/(4.0*c2))/DEG2RAD)
                    kappa = ka * (16.0*c2*c2 - c1*c1) / (4.0* c2)
                    itp_file.append("%-5i %-5i %-5i %3i %f %f ; %-5s %-5s %-5s" % (l_idx + 1, idx_a + 1, r_idx + 1, 1, thetamin, kappa, l_atom.uff_type, central_atom.uff_type, r_atom.uff_type))

        #with open("moffive.top", 'w') as tempfile:
        #    tempfile.writelines(top_file)
        #with open("moffive.gro", 'w') as tempfile:
        #    tempfile.writelines(gro_file)
        #with open("moffive.itp", 'w') as tempfile:
        #    tempfile.writelines(itp_file)
