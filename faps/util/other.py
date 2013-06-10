#!/usr/bin/env python

"""
Small utility functions for general use.

"""

from logging import warning, debug, error, info, critical


def mk_repeat(cube_name='REPEAT_ESP.cube', symmetry=False):
    """Standard REPEAT input file."""
    if symmetry:
        symmetry_flag = 1
    else:
        symmetry_flag = 0
    repeat_input = [
        "Input ESP file name in cube format\n",
        "%s\n" % cube_name,
        "Fit molecular(0) or periodic(1:default) system?\n",
        "1\n",
        "van der Waals scaling factor (default = 1.0)\n",
        "1.00000\n",
        "Apply RESP penalties?, no(0:default), yes(1)\n",
        "0\n",
        "Read cutoff radius? no(0), yes(1:default)\n",
        "1\n",
        "If flag above=1 provide R_cutoff next (in Bohrs)\n",
        "20.00000\n",
        "Apply symmetry restrain? no(0:default), yes(1)\n",
        "%i\n" % symmetry_flag,
        "Use Goddard-type restrains? no(0:default), yes(1)\n",
        "0\n",
        "If flag above=1 then provide weight next\n",
        "0.00000\n",
        "Max distance cut-off scaling factor (1000.0 default)\n",
        "2.0\n"]

    filetemp = open('REPEAT_param.inp', 'w')
    filetemp.writelines(repeat_input)
    filetemp.close()


def mk_connectivity_ff(sym_tree):
    """Write connectivity.ff for input tree."""
    out_tree = ["%i\n" % len(sym_tree)]
    # Always sort to get consistent order
    for idx in sorted(sym_tree):
        out_tree.extend([
            "# %s tree\n" % str(idx),
            "%i\n" % len(sym_tree[idx]),
            " ".join(["%i" % (i+1) for i in sym_tree[idx]]),
            "\n"
        ])

    filetemp = open('connectivity.ff', 'w')
    filetemp.writelines(out_tree)
    filetemp.close()


def mk_incar(options, esp_grid=None):
    """Basic vasp INCAR; use defaults as much as possible."""
    # We need these options
    job_name = options.get('job_name')
    spin = options.getbool('spin')
    optim_h = options.getbool('optim_h')
    optim_all = options.getbool('optim_all')
    optim_cell = options.getbool('optim_cell')
    dispersion = options.getbool('dispersion')

    incar = ["SYSTEM  = %s\n" % job_name,
             "ALGO    = Fast\n",
             "EDIFF   = 1E-5\n",
             "EDIFFG  = -0.02\n",
             "POTIM   = 0.4\n",
             "LREAL   = Auto\n",
             "LVTOT   = .TRUE.\n",
             "LVHAR   = .TRUE.\n",
             "ISMEAR  = 0\n",
             "SIGMA   = 0.05\n",
             "NWRITE  = 0\n"]
    if optim_cell:
        # Positions will be fixed by selective dynamics
        info("Cell vectors will be optimized")
        incar.extend(["ENCUT   = 520\n",
                      "IBRION  = 2\n",
                      "NSW     = 800\n",
                      "ISIF    = 3\n"])
    elif optim_all or optim_h:
        # Just move positions
        incar.extend(["IBRION  = 2\n",
                      "NSW     = 600\n",
                      "ISIF    = 2\n"])
    else:
        # Single point energy
        info("Single point calculation")
        incar.extend(["IBRION  = 0\n",
                      "NSW     = 0\n",
                      "ISIF    = 0\n"])
    if spin:
        info("Spin polarised calculation")
        incar.append("ISPIN   = 2\n")

    if dispersion:
        info("Dispersion correction will be used")
        incar.append("LVDW    = .TRUE.\n")

    if esp_grid is not None:
        info("Changing FFT grid to %ix%ix%i" % esp_grid)
        incar.append("NGXF = %i ; NGYF = %i ; NGZF = %i\n" % esp_grid)

    # VASP recommends, for best performance:
    # NPAR = 4 - approx SQRT( number of cores)
    vasp_ncpu = options.getint('vasp_ncpu')
    npar = vasp_ncpu  # recommended up to 8 cpus
    if vasp_ncpu > 8:
        npar = 4*max(int((vasp_ncpu**0.5)/4.0), 1)
    incar.append("NPAR    = %i\n" % npar)

    return incar


def mk_kpoints(kpoints):
    """Defaults to gamma point only, or specified number."""
    if len(kpoints) != 3:
        error("kpoints specified incorectly; should be (i, i, i)")
    kpoints = [
        "Auto\n",
        "0\n",
        "Gamma\n",
        "%i %i %i\n" % tuple(kpoints),
        "0 0 0\n"]
    return kpoints


def mk_gcmc_control(temperature, pressures, options, guests, supercell=None):
    """Standard GCMC CONTROL file."""
    control = [
        "GCMC Run\n"
        "temperature  %f\n" % temperature,
        "steps    %i\n" % options.getint('mc_prod_steps'),
        "equilibration    %i\n" % options.getint('mc_eq_steps'),
        "cutoff          %f angstrom\n" % options.getfloat('mc_cutoff'),
        "delr            1.0 angstrom\n",
        "ewald precision  1d-6\n",
        "numguests %i\n" % options.getint('mc_numguests_freq'),
        "history %i\n" % options.getint('mc_history_freq')]
    if options.getbool('mc_jobcontrol'):
        control.append("jobcontrol\n")
    else:
        control.append("# jobcontrol\n")
    # Guest stuff
    if options.getbool('mc_probability_plot'):
        # We just comment out the probabilty plot lines if unneeded
        prob_on = ""
    else:
        prob_on = "# "
    if len(guests) > 1:
        gp_zip = zip(guests, pressures)
    else:
        gp_zip = zip(guests, [pressures])
    guest_count = 0
    for guest, press in gp_zip:
        guest_count += 1
        control.append("&guest %i\n" % guest_count)
        control.append("  pressure (bar) %f\n" % press)
        control.append("%s  probability %i\n" %
                       (prob_on, len(guest.probability)))
        for prob in guest.probability:
            control.append("%s  %i  " % (prob_on, len(prob)) +
                           "  ".join(["%i" % x for x in prob]) + "\n")
        control.append("&end\n")

    if options.getbool('fold') and supercell is not None:
        control.append('\ngrid factors %i %i %i\n' % supercell)
    elif supercell is not None:
        control.append('\n# grid factors %i %i %i\n' % supercell)

    resolution = options.getfloat('mc_probability_plot_spacing')
    if resolution != FASTMC_DEFAULT_GRID_SPACING:
        control.append('\ngrid spacing %f\n' % resolution)
    else:
        control.append('\n#grid spacing %f\n' % FASTMC_DEFAULT_GRID_SPACING)


    control.append("\nfinish\n")
    return control


def parse_qeq_params(param_tuple):
    """Convert an options tuple to a dict of values."""
    # group up into ((atom, electronegativity, 0.5*hardness), ... )
    param_tuple = subgroup(param_tuple, 3)
    param_dict = {}
    for param_set in param_tuple:
        try:
            atom_type = int(param_set[0])
            atom_type = ATOMIC_NUMBER[atom_type]
        except ValueError:
            # assume it is an element symbol instead
            atom_type = param_set[0]
        except IndexError:
            # Typed atom, ignore for now
            continue
        try:
            param_dict[atom_type] = (float(param_set[1]), float(param_set[2]))
        except IndexError:
            warning("Cannot read parameters for %s" % atom_type)

    return param_dict


def mk_egulp_params(param_tuple):
    """Convert an options tuple to an EGULP parameters file filling."""
    # group up into ((atom, electronegativity, 0.5*hardness), ... )
    param_tuple = subgroup(param_tuple, 3)
    # first line is the number of parametr sets
    #TODO(tdaff): try adding some error checking
    egulp_file = ["%s\n" % len(param_tuple)]
    for paramset in param_tuple:
        try:
            # catch if it is an atomic number or element symbol
            atomic_number = int(paramset[0])
        except ValueError:
            atomic_number = ATOMIC_NUMBER.index(paramset[0])
        egulp_file.append("%-4i %f %f\n" % (atomic_number, float(paramset[1]),
                                            float(paramset[2])))

    return egulp_file


def mk_egulp_ini(options):
    """Create a default ini file for egulp to do qeq calculation."""
    egulp_grid = int(options.getbool('egulp_grid'))
    egulp_potential = int(options.getbool('egulp_potential'))
    egulp_potential_difference = int(options.getbool('egulp_potential_difference'))

    egulp_grid_parameters = options.get('egulp_grid_parameters')

    egulp_ini = [
        "build_grid %i\n" % egulp_grid,
        "build_grid_from_scratch %s\n" % egulp_grid_parameters,
        "save_grid %i grid.cube\n" % egulp_grid,
        "calculate_pot_diff %i\n" % egulp_potential_difference,
        "calcaulte_pot %i repeat.cube\n" % egulp_potential,
        "skip_everything 0\n",
        "point_charges_present 0\n",
        "include_pceq 0\n",
        "imethod 0\n"]

    return egulp_ini


def validate_gulp_output(filename):
    """Check to see if gulp calculation has finished and return the energy."""

    final_energy = None
    final_gnorm = None

    try:
        gulp_output = open(filename)
    except IOError:
        warning("Gulp output not found; continuing anyway")
        return (False, final_energy, final_gnorm)

    finished_optimisation = False

    for line in gulp_output:
        if line.startswith("  Cycle:"):
            line = line.split()
            try:
                final_energy = float(line[3])
                final_gnorm = float(line[5])
            except ValueError:
                final_energy = 999999.9
                final_gnorm = 999999.9
        elif "Optimisation achieved" in line:
            # Great
            finished_optimisation = True
            break
        elif "Too many failed attempts to optimise" in line:
            warning("Gulp reports failed optimisation; check gnorm!")
            finished_optimisation = True
            break
        elif "no lower point can be found" in line:
            warning("Gulp reports failed optimisation; check gnorm!")
            finished_optimisation = True
            break

    if not finished_optimisation:
        warning("Gulp optimisation did not finish; check output!")

    if final_gnorm > 0.1:
        warning("Gulp optimisation has gnorm > 0.1; check output!")
        finished_optimisation = False

    debug("Gulp .out: %s" % [finished_optimisation, final_energy, final_gnorm])
    return (finished_optimisation, final_energy, final_gnorm)


def lorentz_berthelot(left, right):
    """Lorentz-Berthelot mixing rules for (sigma, epsilon) tuples."""
    sigma = (left[0] + right[0]) / 2.0
    epsilon = (left[1] * right[1])**0.5
    if epsilon == 0:
        sigma = 0
    return sigma, epsilon


def same_guests(base, other):
    """Test if the guests are the same index and order."""
    return [item.ident for item in base] == [item.ident for item in other]


def prod(seq):
    """Calculate the product of all members of a sequence."""
    # numpy.prod will silently overflow 32 bit integer values
    # so we can use python bignums natively
    product = 1
    for item in seq:
        product *= item
    return product


def dot3(vec1, vec2):
    """Calculate dot product for two 3d vectors."""
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]


def min_distance(first_atom, second_atom, cell=None):
    """Helper to find mimimum image criterion distance."""
    if cell is None:
        cell = first_atom._parent.cell.cell
    return min_dist(first_atom.pos,
                    first_atom.fractional,
                    second_atom.pos,
                    second_atom.fractional,
                    cell)


def min_dist(c_coa, f_coa, c_cob, f_cob_in, box):
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
        return vecdist3(c_coa, c_cob)
    else:
        new_b = [f_cob[0]*box[0][0] + f_cob[1]*box[1][0] + f_cob[2]*box[2][0],
                 f_cob[0]*box[0][1] + f_cob[1]*box[1][1] + f_cob[2]*box[2][1],
                 f_cob[0]*box[0][2] + f_cob[1]*box[1][2] + f_cob[2]*box[2][2]]
        return vecdist3(c_coa, new_b)


def vecdist3(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    vec = [coord2[0] - coord1[0],
           coord2[1] - coord1[1],
           coord2[2] - coord1[2]]

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5


def subgroup(iterable, width, itype=None):
    """Split an iterable into nested sub-itypes of width members."""
    # Return the same type as iterable
    if itype is None:
        if isinstance(iterable, list):
            itype = list
        else:
            itype = tuple
    # Will leave short groups if not enough members
    return itype([itype(iterable[x:x+width])
                  for x in range(0, len(iterable), width)])


def state_points(temperatures, pressures, individual, nguests):
    """Group temperatures and pressures and append given state point tuples."""
    # No error checking done, assume know what they are doing
    for index in range(0, len(individual), nguests+1):
        yield (individual[index], tuple(individual[index+1:index+nguests+1]))
    for temp in temperatures:
        for pressure in subgroup(pressures, nguests):
            yield (temp, pressure)


def name_from_types(sites, guest):
    """Generate a string that gives the atoms represented in sites."""
    stypes = []
    for site_ident in sites:
        if site_ident is 0:
            stypes.append('COM')
        else:
            stypes.append(guest.atoms[site_ident-1].element)
    stypes = unique(stypes)
    if len(stypes) > 1:
        site_name = "-".join(stypes)
    else:
        site_name = stypes[0]
    return site_name


def other_bond_index(bond, index):
    """Return the atom index for the other atom in a bond."""
    if bond[0] == index:
        return bond[1]
    elif bond[1] == index:
        return bond[0]
    else:
        raise ValueError("Index %s not found in bond %s" % (index, bond))


def welcome():
    """Print any important messages."""
    print(LOGO)
    print(("faps %s" % __version__).rjust(79))

