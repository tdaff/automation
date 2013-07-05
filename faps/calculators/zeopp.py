    def to_zeoplusplus(self):
        """
        Return a tuple containing a cssr file, radii file and mass file.

        """

        radii = ["%-7s %-f\n" % (atom, UFF[atom][0]/2.0)
                 for atom in unique(self.types)]
        masses = ["%-7s %-f\n" % (atom, WEIGHT[atom])
                  for atom in unique(self.types)]

        return (self.to_cssr(no_atom_id=True), radii, masses)

