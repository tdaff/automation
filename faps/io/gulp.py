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
from . import StructureReader, StructureWriter


class GulpReader(StructureReader):
    def from_gulp_output(self, filename):
        """Update the structure from the gulp optimisation output."""
        info("Updating positions from file: %s" % filename)
        grs_out = open(filename)
        for line in grs_out:
            if line.strip() == 'cell':
                params = tuple(float(x) for x in grs_out.next().split()[:6])
                self.cell.params = params
            elif line.strip() == 'fractional':
                cell = self.cell.cell
                for atom, atom_line in zip(self.atoms, grs_out):
                    atom.pos = dot([gfloat(x) for x in atom_line.split()[2:5]], cell)
                    # FIXME(tdaff) fractionals need to change automatically
                    # when the position gets updated (or the cell!)
                    del atom.fractional

        # Make sure everything is good from here
        if self.check_close_contacts(covalent=1.0):
            warning("Structure might have atom overlap, check gulp output!")
            self.bad_structure = True

        if self.bond_length_check():
            warning("Structure might have strained bonds, check gulp output!")
            self.bad_structure = True



    def charges_from_gulp(self, filename):
        """Parse QEq charges from GULP output."""
        info("Getting charges from file: %s" % filename)
        filetemp = open(filename)
        gout = filetemp.readlines()
        filetemp.close()
        start_line = None
        failures = 0
        for line_num, line in enumerate(gout):
            if '  Final charges from QEq :' in line:
                # Take the last set of charges
                start_line = line_num + 7
            elif 'Failed to converge' in line:
                # This will occur twice in the output for complete failure
                failures += 1
        if start_line is None:
            error("Final charges not found in gulp output")
            terminate(184)
        elif failures > 1:
            warning("Gulp charges may not be converged")
        for atom, chg_line in zip(self.atoms, gout[start_line:]):
            atom.charge = float(chg_line.split()[2])



class GulpWriter(StructureWriter):
    """
    CIF file reader with some symmetry and custom properties.

    """

    filetypes = ['*.cif', 'cif']

    @staticmethod
    def to_gulp(self, qeq_fit=False, optimise=False, terse=False, qeq_dict={}):
        """Return a GULP file to use for the QEq charges."""
        if qeq_fit:
            keywords = "fitting bulk_noopt qeq\n"
        elif optimise:
            return self.to_gulp_optimise(terse=terse)
        else:
            keywords = "single conp qeq\n"

        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# Initial file written by faps\n"
            "name %s\n" % self.name,
            "vectors\n"] + self.cell.to_vector_strings() + [
            "cartesian\n"]
        if qeq_fit:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                                 "%14.7f\n" % atom.charge])
            gin_file.append("\n# \n# Fitting variables:\n# \nobservables\n")
            for idx, atom in enumerate(self.atoms):
                gin_file.append("monopoleq\n%5s%11.6f\n" % (idx+1, atom.charge))
            gin_file.append("end\nqelectronegativity\n")
            for at_type in unique(self.types):
                gin_file.extend(
                    ["%-5s" % at_type,
                     "%9.5f %9.5f %9.5f 1 1 0\n" % QEQ_PARAMS[at_type]])
        else:
            for atom in self.atoms:
                gin_file.extend(["%-5s core " % atom.type,
                                 "%14.7f %14.7f %14.7f\n" % tuple(atom.pos)])
            if qeq_dict:
                unique_types = unique(self.types)
                gin_file.append('\nqelectronegativity\n')
                for atom_type, params in qeq_dict.items():
                    if atom_type in unique_types:
                        gin_file.append('%-4s %f %f\n' % (atom_type, params[0], params[1]))

        gin_file.append("\ndump every %s.grs\nprint 1\n" % self.name)
        return gin_file



    def to_gulp_optimise(self, terse=False):
        """Return a GULP file to optimise with UFF."""
        # all atoms of the same forcefield type must have the same label
        # gulp fails if atoms of the same type have different labels!
        # this only crops up as an error in the 4.0 versions
        # 'decimal_only' stops fractions that cannot be parsed when reading
        # updated positions
        if terse:
            # Don't output the bonds
            keywords = "opti noautobond decimal_only conj\n"
        else:
            keywords = "opti noautobond bond decimal_only conj\n"
        gin_file = [
            "# \n# Keywords:\n# \n",
            keywords,
            "# \n# Options:\n# \n",
            "# UFF optimisation by gulp\n"
            "name %s\n" % self.name,
            # using an rfo minimiser with a preconditioned hessian from the
            # BFGS minimiser seems to be the most efficient
            "switch rfo gnorm 0.3\n",
            "vectors\n"] + self.cell.to_vector_strings() + [
            " 1 1 1\n 1 1 1\n 1 1 1\n",  # constant pressure relaxation
            "cartesian\n"]

        all_ff_types = {}

        for at_idx, atom in enumerate(self.atoms):
            ff_type = atom.uff_type
            if not ff_type in all_ff_types:
                all_ff_types[ff_type] = atom.site
                # Sanity check in case types are not given in the input
                if ff_type not in UFF_FULL:
                    error("Atom %s has unknown type %s" % (atom, ff_type))
            if atom.is_fixed:
                fixed_flags = "0 0 0"
            else:
                fixed_flags = "1 1 1"
            gin_file.extend(["%-5s core " % all_ff_types[ff_type],
                             "%14.7f %14.7f %14.7f " % tuple(atom.pos),
                             "%f " % atom.charge,
                             "%f " % 1.0,  # occupancy
                             fixed_flags, "\n"])


        #identify all the individual uff species for the library
        gin_file.append("\nspecies\n")
        for ff_type, species in all_ff_types.items():
            gin_file.append("%-6s core %-6s\n" % (species, ff_type))

        gin_file.append("\n")
        for bond in sorted(self.bonds):
            bond_type = GULP_BOND_ORDERS[self.bonds[bond]]
            gin_file.append("connect %6i %6i %s\n" % (bond[0] + 1, bond[1] + 1, bond_type))

        gin_file.append("\nlibrary uff\n")

        gin_file.append("\nstepmx opt 0.05\n")

        # Restart file is for final structure
        gin_file.append("\ndump every %s.grs\n" % self.name)
        # optimization movie useful for debugging mostly
        if terse:
            # These tell gulp to be quiet, but we also stop the massive arc
            # file being generated
            gin_file.append("\nterse inout structure\n"
                            "terse inout potentials\n"
                            "terse inout derivatives\n"
                            "#output movie arc %s\n" % self.name)
        else:
            gin_file.append("\noutput movie arc %s\n" % self.name)

        return gin_file


