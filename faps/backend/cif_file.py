"""
Cif file wriing backend for fapswitch.

"""

import hashlib


class CifFileBackend(object):
    """Abstraction for writing cif files in a pluggable manner."""

    def add_symmetry_structure(self, base_structure, functions, cif_file):
        """
        Write out the cif file with a name derived from the base structure
        and the functionalisations.

        """
        new_mof_name = ".".join(["@".join(x) for x in functions])
        cif_filename = '%s_func_%s.cif' % (base_structure, new_mof_name)
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)

    def add_freeform_structure(self, base_structure, functions, cif_file):
        """
        Write a cif file with an md5 fixed-length name based on the
        functionalisation.

        """
        unique_name = hashlib.md5(str(functions)).hexdigest()
        cif_filename = '%s_free_%s.cif' % (base_structure, unique_name)
        with open(cif_filename, 'w') as output_file:
            output_file.writelines(cif_file)
