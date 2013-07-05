    def fastmc_postproc(self, filepath, tp_point, options):
        """Update structure properties from gcmc OUTPUT."""
        startdir = os.getcwd()
        os.chdir(filepath)

        # Since Pete changed output strip indentation and blank lines
        filetemp = open('OUTPUT')
        output = strip_blanks(filetemp.readlines())
        filetemp.close()

        # Keep track of supercell so we can get unit cell values
        supercell_mult = prod(self.gcmc_supercell)
        # Still positional as we need multiple values simultaneously
        # and very old versions changed wording of heat of adsorption
        # and enthalpy of guest
        # TODO(tdaff, r2.0): deprecate reading older fastmc files
        # and put Cv in the guest definition
        for line in output[::-1]:
            if "+/-" in line:
                # This is version 1.3 of fastmc
                debug("NEW OUTPUT")
                line_offset = 5
                read_cv = True
                # In future this should be assumed to exist
                for guest in self.guests:
                    if not hasattr(guest, 'c_v'):
                        guest.c_v = {}
                break
        else:
            # +/- not found, assume old style output
            debug("OLD OUTPUT")
            line_offset = 0
            read_cv = False
        for idx, line in enumerate(output):
            # Assume that block will always start like this
            if 'final stats' in line:
                idx += line_offset
                guest_id = int(line.split()[4]) - 1
                self.guests[guest_id].uptake[tp_point] = (
                    float(output[idx + 1].split()[-1]),
                    float(output[idx + 2].split()[-1]),
                    supercell_mult)
                # This will sometimes be NaN
                self.guests[guest_id].hoa[tp_point] = (
                    float(output[idx + 3].split()[-1]),
                    float(output[idx + 4].split()[-1]))
                if read_cv:
                    self.guests[guest_id].c_v[tp_point] = (
                    float(output[idx + 5].split()[-1]),
                    float(output[idx + 6].split()[-1]))
            elif 'total accepted steps' in line:
                counted_steps = int(line.split()[-1])
                if counted_steps < 10000:
                    warning("Number of accepted GCMC steps is very low; "
                            "only %i counted!" % counted_steps)

        fold = options.getbool('fold')
        find_maxima = options.getbool('find_maxima')
        prob_plot = options.getbool('mc_probability_plot')
        folded = False
        if prob_plot and (fold or find_maxima):
            folded = self.fold_and_maxima(fold, find_maxima, tp_point)

        if folded and not options.getbool('fastmc_keep_unfolded_cubes'):
            debug("Removing unfolded cube files")
            cubes = glob.glob("prob_guest??_prob_??.cube")
            remove_files(cubes)

        #TODO(tdaff): absl

        unneeded_files = options.gettuple('fastmc_delete_files')
        remove_files(unneeded_files)
        keep_files = options.gettuple('fastmc_compress_files')
        compress_files(keep_files)

        os.chdir(startdir)


    def fold_and_maxima(self, fold=True, find_maxima=True, tp_point=None):
        """Determine the positions of maxima and produce an xyz xyz file."""
        from cube import Cube
        folded = False
        if fold:
            fold = self.gcmc_supercell
        else:
            fold = None
        for guest_idx, guest in enumerate(self.guests):
            guest_locations = {}
            for site_idx, sites in enumerate(guest.probability):
                guest_cube = Cube("prob_guest%02i_prob_%02i.cube" %
                                  (guest_idx+1, site_idx+1), fold=fold)
                if fold is not None:
                    debug("Folded cube file: %s" % guest_cube.folded_name)
                    guest_cube.write_cube()
                    folded = True
                if find_maxima:
                    guest_locations[sites] = guest_cube.maxima()
            if guest_locations:
                if tp_point:
                    # We can keep them for later too, must create dict
                    # since it might not exist in old calculations
                    if not hasattr(guest, 'guest_locations'):
                        guest.guest_locations = {tp_point: guest_locations}
                    else:
                        guest.guest_locations[tp_point] = guest_locations
                maxima = []
                for sites in guest_locations:
                    atom_name = name_from_types(sites, guest)
                    for atom, magnitude in guest_locations[sites]:
                        maxima.append((magnitude, "%-6s" % atom_name +
                                      "%10.6f %10.6f %10.6f " % tuple(atom) +
                                      "%10.6f\n" % magnitude))
                locations = open("%s-%s.xyz" % (self.name, guest.ident), 'w')
                locations.write(" %i\nBinding sites at %r\n" %
                                (len(maxima), tp_point))
                locations.writelines([x[1] for x in sorted(maxima, reverse=True)])
                locations.close()

        return folded

