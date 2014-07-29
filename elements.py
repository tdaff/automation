"""
elements.py

Some standard element properties. Import as global constants.

"""

WEIGHT = {
    "H": 1.00794,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.96,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.9054519,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.500,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Th": 232.03806,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 262,
    "Rf": 265,
    "Db": 268,
    "Sg": 271,
    "Bh": 272,
    "Hs": 270,
    "Mt": 276,
    "Ds": 281,
    "Rg": 280,
    "Cn": 285,
    "Uut": 284,
    "Uuq": 289,
    "Uup": 288,
    "Uuh": 293,
    "Uuo": 294}

# If this is a nice list we can just .index or [slice] to get atomic numbers
ATOMIC_NUMBER = [
    "ZERO", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uuo"]

UFF = {
    "H": (2.5711, 0.0440),
    "He": (2.1043, 0.0560),
    "Li": (2.1836, 0.0250),
    "Be": (2.4455, 0.0850),
    "B": (3.6375, 0.1800),
    "C": (3.4309, 0.1050),
    "N": (3.2607, 0.0690),
    "O": (3.1181, 0.0600),
    "F": (2.9970, 0.0500),
    "Ne": (2.8892, 0.0420),
    "Na": (2.6576, 0.0300),
    "Mg": (2.6914, 0.1110),
    "Al": (4.0082, 0.5050),
    "Si": (3.8264, 0.4020),
    "P": (3.6946, 0.3050),
    "S": (3.5948, 0.2740),
    "Cl": (3.5164, 0.2270),
    "Ar": (3.4460, 0.1850),
    "K": (3.3961, 0.0350),
    "Ca": (3.0282, 0.2380),
    "Sc": (2.9355, 0.0190),
    "Ti": (2.8286, 0.0170),
    "V": (2.8010, 0.0160),
    "Cr": (2.6932, 0.0150),
    "Mn": (2.6380, 0.0130),
    "Fe": (2.5943, 0.0130),
    "Co": (2.5587, 0.0140),
    "Ni": (2.5248, 0.0150),
    "Cu": (3.1137, 0.0050),
    "Zn": (2.4616, 0.1240),
    "Ga": (3.9048, 0.4150),
    "Ge": (3.8130, 0.3790),
    "As": (3.7685, 0.3090),
    "Se": (3.7462, 0.2910),
    "Br": (3.7320, 0.2510),
    "Kr": (3.6892, 0.2200),
    "Rb": (3.6652, 0.0400),
    "Sr": (3.2438, 0.2350),
    "Y": (2.9801, 0.0720),
    "Zr": (2.7832, 0.0690),
    "Nb": (2.8197, 0.0590),
    "Mo": (2.7190, 0.0560),
    "Tc": (2.6709, 0.0480),
    "Ru": (2.6397, 0.0560),
    "Rh": (2.6094, 0.0530),
    "Pd": (2.5827, 0.0480),
    "Ag": (2.8045, 0.0360),
    "Cd": (2.5373, 0.2280),
    "In": (3.9761, 0.5990),
    "Sn": (3.9128, 0.5670),
    "Sb": (3.9378, 0.4490),
    "Te": (3.9823, 0.3980),
    "I": (4.0090, 0.3390),
    "Xe": (3.9235, 0.3320),
    "Cs": (4.0242, 0.0450),
    "Ba": (3.2990, 0.3640),
    "La": (3.1377, 0.0170),
    "Ce": (3.1680, 0.0130),
    "Pr": (3.2126, 0.0100),
    "Nd": (3.1850, 0.0100),
    "Pm": (3.1600, 0.0090),
    "Sm": (3.1360, 0.0080),
    "Eu": (3.1119, 0.0080),
    "Gd": (3.0005, 0.0090),
    "Tb": (3.0745, 0.0070),
    "Dy": (3.0540, 0.0070),
    "Ho": (3.0371, 0.0070),
    "Er": (3.0210, 0.0070),
    "Tm": (3.0059, 0.0060),
    "Yb": (2.9890, 0.2280),
    "Lu": (3.2429, 0.0410),
    "Hf": (2.7983, 0.0720),
    "Ta": (2.8241, 0.0810),
    "W": (2.7342, 0.0670),
    "Re": (2.6317, 0.0660),
    "Os": (2.7796, 0.0370),
    "Ir": (2.5302, 0.0730),
    "Pt": (2.4535, 0.0800),
    "Au": (2.9337, 0.0390),
    "Hg": (2.4099, 0.3850),
    "Tl": (3.8727, 0.6800),
    "Pb": (3.8282, 0.6630),
    "Bi": (3.8932, 0.5180),
    "Po": (4.1952, 0.3250),
    "At": (4.2318, 0.2840),
    "Rn": (4.2451, 0.2480),
    "Fr": (4.3654, 0.0500),
    "Ra": (3.2758, 0.4040),
    "Ac": (3.0985, 0.0330),
    "Th": (3.0255, 0.0260),
    "Pa": (3.0504, 0.0220),
    "U": (3.0246, 0.0220),
    "Np": (3.0504, 0.0190),
    "Pu": (3.0504, 0.0160),
    "Am": (3.0121, 0.0140),
    "Cm": (2.9631, 0.0130),
    "Bk": (2.9747, 0.0130),
    "Cf": (2.9515, 0.0130),
    "Es": (2.9391, 0.0120),
    "Fm": (2.9275, 0.0120),
    "Md": (2.9168, 0.0110),
    "No": (2.8936, 0.0110),
    "Lr": (2.8829, 0.0110)
}

VASP_PSEUDO_PREF = {
    "Li": "Li_sv",
    "Na": "Na_pv",
    "K": "K_sv",
    "Ca": "Ca_pv",
    "Sr": "Sr_sv",
    "Rb": "Rb_sv",
    "Cs": "Cs_sv",
    "Ba": "Ba_sv",
    "Ti": "Ti_pv",
    "V": "V_pv",
    "Cr": "Cr_pv",
    "Mn": "Mn_pv",
    "Sc": "Sc_sv",
    "Y": "Y_sv",
    "Zr": "Zr_sv",
    "Nb": "Nb_pv",
    "Mo": "Mo_pv",
    "Tc": "Tc_pv",
    "Hf": "Hf_pv",
    "Ta": "Ta_pv",
    "W": "W_pv",
    "Os": "Os_pv",
    "Ga": "Ga_d",
    "Ge": "Ge_d",
    "In": "In_d",
    "Sn": "Sn_d",
    "Tl": "Tl_d",
    "Pb": "Pb_d",
    "Bi": "Bi_d",
    "Po": "Po_d",
    "At": "At_d"
}

#TODO(jlo): make a dict with the pseudopotential valence electrons

QEQ_PARAMS = {
    # (chi, mu, rad)
    "H": (4.528, 6.9452, 0.371),
    "He": (9.66, 14.92, 1.3),
    "Li": (3.006, 2.386, 1.557),
    "Be": (4.877, 4.443, 1.24),
    "B": (5.11, 4.75, 0.822),
    "C": (5.343, 5.063, 0.759),
    "N": (6.899, 5.88, 0.715),
    "O": (8.741, 6.682, 0.669),
    "F": (10.874, 7.474, 0.706),
    "Ne": (11.04, 10.55, 1.768),
    "Na": (2.843, 2.296, 2.085),
    "Mg": (3.951, 3.693, 1.5),
    "Al": (4.06, 3.59, 1.201),
    "Si": (4.168, 3.487, 1.176),
    "P": (5.463, 4, 1.102),
    "S": (6.928, 4.486, 1.047),
    "Cl": (8.564, 4.946, 0.994),
    "Ar": (9.465, 6.355, 2.108),
    "K": (2.421, 1.92, 2.586),
    "Ca": (3.231, 2.88, 2),
    "Sc": (3.395, 3.08, 1.75),
    "Ti": (3.47, 3.38, 1.607),
    "V": (3.65, 3.41, 1.47),
    "Cr": (3.415, 3.865, 1.402),
    "Mn": (3.325, 4.105, 1.533),
    "Fe": (3.76, 4.14, 1.393),
    "Co": (4.105, 4.175, 1.406),
    "Ni": (4.465, 4.205, 1.398),
    "Cu": (3.729, 2.501, 1.434),
    "Zn": (5.106, 4.285, 1.4),
    "Ga": (3.641, 3.16, 1.211),
    "Ge": (4.051, 3.438, 1.189),
    "As": (5.188, 3.809, 1.204),
    "Se": (6.428, 4.131, 1.224),
    "Br": (7.79, 4.425, 1.141),
    "Kr": (8.505, 5.715, 2.27),
    "Rb": (2.331, 1.846, 2.77),
    "Sr": (3.024, 2.44, 2.415),
    "Y": (3.83, 2.81, 1.998),
    "Zr": (3.4, 3.55, 1.758),
    "Nb": (3.55, 3.38, 1.603),
    "Mo": (3.465, 3.755, 1.53),
    "Tc": (3.29, 3.99, 1.5),
    "Ru": (3.575, 4.015, 1.5),
    "Rh": (3.975, 4.005, 1.509),
    "Pd": (4.32, 4, 1.544),
    "Ag": (4.436, 3.134, 1.622),
    "Cd": (5.034, 3.957, 1.6),
    "In": (3.506, 2.896, 1.404),
    "Sn": (3.987, 3.124, 1.354),
    "Sb": (4.899, 3.342, 1.404),
    "Te": (5.816, 3.526, 1.38),
    "I": (6.822, 3.762, 1.333),
    "Xe": (7.595, 4.975, 2.459),
    "Cs": (2.183, 1.711, 2.984),
    "Ba": (2.814, 2.396, 2.442),
    "La": (2.8355, 2.7415, 2.071),
    "Ce": (2.744, 2.692, 1.925),
    "Pr": (2.858, 2.564, 2.007),
    "Nd": (2.8685, 2.6205, 2.007),
    "Pm": (2.881, 2.673, 2),
    "Sm": (2.9115, 2.7195, 1.978),
    "Eu": (2.8785, 2.7875, 2.227),
    "Gd": (3.1665, 2.9745, 1.968),
    "Tb": (3.018, 2.834, 1.954),
    "Dy": (3.0555, 2.8715, 1.934),
    "Ho": (3.127, 2.891, 1.925),
    "Er": (3.1865, 2.9145, 1.915),
    "Tm": (3.2514, 2.9329, 2),
    "Yb": (3.2889, 2.965, 2.158),
    "Lu": (2.9629, 2.4629, 1.896),
    "Hf": (3.7, 3.4, 1.759),
    "Ta": (5.1, 2.85, 1.605),
    "W": (4.63, 3.31, 1.538),
    "Re": (3.96, 3.92, 1.6),
    "Os": (5.14, 3.63, 1.7),
    "Ir": (5, 4, 1.866),
    "Pt": (4.79, 4.43, 1.557),
    "Au": (4.894, 2.586, 1.618),
    "Hg": (6.27, 4.16, 1.6),
    "Tl": (3.2, 2.9, 1.53),
    "Pb": (3.9, 3.53, 1.444),
    "Bi": (4.69, 3.74, 1.514),
    "Po": (4.21, 4.21, 1.48),
    "At": (4.75, 4.75, 1.47),
    "Rn": (5.37, 5.37, 2.2),
    "Fr": (2, 2, 2.3),
    "Ra": (2.843, 2.434, 2.2),
    "Ac": (2.835, 2.835, 2.108),
    "Th": (3.175, 2.905, 2.018),
    "Pa": (2.985, 2.905, 1.8),
    "U": (3.341, 2.853, 1.713),
    "Np": (3.549, 2.717, 1.8),
    "Pu": (3.243, 2.819, 1.84),
    "Am": (2.9895, 3.0035, 1.942),
    "Cm": (2.8315, 3.1895, 1.9),
    "Bk": (3.1935, 3.0355, 1.9),
    "Cf": (3.197, 3.101, 1.9),
    "Es": (3.333, 3.089, 1.9),
    "Fm": (3.4, 3.1, 1.9),
    "Md": (3.47, 3.11, 1.9),
    "No": (3.475, 3.175, 1.9),
    "Lr": (3.5, 3.2, 1.9)
}

CCDC_BOND_ORDERS = {
    # http://cima.chem.usyd.edu.au:8080/cif/skunkworks/html/ddl1/mif/bond.html
    'S': 1.0,  # single (two-electron) bond or sigma bond to metal
    'D': 2.0,  # double (four-electron) bond
    'T': 3.0,  # triple (six-electron) bond
    'Q': 4.0,  # quadruple (eight-electron, metal-metal) bond
    'A': 1.5,  # alternating normalized ring bond (aromatic)
    'C': 1.0,  # catena-forming bond in crystal structure
    'E': 1.5,  # equivalent (delocalized double) bond
    'P': 1.0,   # pi bond (metal-ligand pi interaction)
    'Am': 1.41,  # Amide bond (non standard)
    1.0: 'S',  # single (two-electron) bond or sigma bond to metal
    2.0: 'D',  # double (four-electron) bond
    3.0: 'T',  # triple (six-electron) bond
    4.0: 'Q',  # quadruple (eight-electron, metal-metal) bond
    1.5: 'A',  # alternating normalized ring bond (aromatic)
    1.41: 'Am'  # Amide bond (non standard)
}

GULP_BOND_ORDERS = {1: 'single', 1.41: 'amide', 1.5: 'resonant',
                    2: 'double', 3: 'triple', 4: 'quadruple'}

OB_BOND_ORDERS = {
    # openbabel needs integers and aromatics are 5?
    1.0: 1,  # single (two-electron) bond or sigma bond to metal
    2.0: 2,  # double (four-electron) bond
    3.0: 3,  # triple (six-electron) bond
    4.0: 4,  # quadruple (eight-electron, metal-metal) bond
    1.5: 5,  # alternating normalized ring bond (aromatic)
    1.41: 5  # Amide bond (non standard)
}

UFF_TYPES = {
    "H": "H_", #"H_b",
    "He": "He4+4",
    "Li": "Li",
    "Be": "Be3+2",
    "B": "B_3", #"B_2",
    "C": "C_3", #"C_R", "C_2", "C_1",
    "N": "N_3", #"N_R", "N_2", "N_1",
    "O": "O_3", #"O_3_z", "O_R", "O_2", "O_1",
    "F": "F_",
    "Ne": "Ne4+4",
    "Na": "Na",
    "Mg": "Mg3+2",
    "Al": "Al3",
    "Si": "Si3",
    "P": "P_3+3", #"P_3+5", "P_3+q",
    "S": "S_3+2", #"S_3+4", "S_3+6", "S_R", "S_2",
    "Cl": "Cl",
    "Ar": "Ar4+4",
    "K": "K",
    "Ca": "Ca6+2",
    "Sc": "Sc3+3",
    "Ti": "Ti3+4", #"Ti6+4",
    "V": "V_3+5",
    "Cr": "Cr6+3",
    "Mn": "Mn6+2",
    "Fe": "Fe3+2", #"Fe6+2",
    "Co": "Co6+2",
    "Ni": "Ni4+2",
    "Cu": "Cu3+1",
    "Zn": "Zn3+2",
    "Ga": "Ga3+3",
    "Ge": "Ge3",
    "As": "As3+3",
    "Se": "Se3+2",
    "Br": "Br",
    "Kr": "Kr4+4",
    "Rb": "Rb",
    "Sr": "Sr6+2",
    "Y": "Y_3+3",
    "Zr": "Zr3+4",
    "Nb": "Nb3+5",
    "Mo": "Mo6+6", #"Mo3+6",
    "Tc": "Tc6+5",
    "Ru": "Ru6+2", #"Rh6+3",
    "Pd": "Pd4+2",
    "Ag": "Ag1+1",
    "Cd": "Cd3+2",
    "In": "In3+3",
    "Sn": "Sn3",
    "Sb": "Sb3+3",
    "Te": "Te3+2",
    "I": "I_",
    "Xe": "Xe4+4",
    "Cs": "Cs",
    "Ba": "Ba6+2",
    "La": "La3+3",
    "Ce": "Ce6+3",
    "Pr": "Pr6+3",
    "Nd": "Nd6+3",
    "Pm": "Pm6+3",
    "Sm": "Sm6+3",
    "Eu": "Eu6+3",
    "Gd": "Gd6+3",
    "Tb": "Tb6+3",
    "Dy": "Dy6+3",
    "Ho": "Ho6+3",
    "Er": "Er6+3",
    "Tm": "Tm6+3",
    "Yb": "Yb6+3",
    "Lu": "Lu6+3",
    "Hf": "Hf3+4",
    "Ta": "Ta3+5",
    "W": "W_6+6", #"W_3+4", "W_3+6",
    "Re": "Re6+5", #"Re3+7",
    "Os": "Os6+6",
    "Ir": "Ir6+3",
    "Pt": "Pt4+2",
    "Au": "Au4+3",
    "Hg": "Hg1+2",
    "Tl": "Tl3+3",
    "Pb": "Pb3",
    "Bi": "Bi3+3",
    "Po": "Po3+2",
    "At": "At",
    "Rn": "Rn4+4",
    "Fr": "Fr",
    "Ra": "Ra6+2",
    "Ac": "Ac6+3",
    "Th": "Th6+4",
    "Pa": "Pa6+4",
    "U": "U_6+4",
    "Np": "Np6+4",
    "Pu": "Pu6+4",
    "Am": "Am6+4",
    "Cm": "Cm6+3",
    "Bk": "Bk6+3",
    "Cf": "Cf6+3",
    "Es": "Es6+3",
    "Fm": "Fm6+3",
    "Md": "Md6+3",
    "No": "No6+3",
    "Lr": "Lr6+3",
}

METALS = [3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
          37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57,
          58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
          75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93, 94,
          95, 96, 97, 98, 99, 100, 101, 102, 103]

COVALENT_RADII = {
    # Covalent radii revisited -- DOI:10.1039/B801115J
    "H": 0.31,
    "He": 0.28,
    "Li": 1.28,
    "Be": 0.96,
    "B": 0.84,
    "C": 0.76,  # for sp3; sp2 = 0.73; sp = 0.69
    "C_1": 0.69,
    "C_2": 0.73,
    "C_R": 0.73,
    "C_3": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Ne": 0.58,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ar": 1.06,
    "K": 2.03,
    "Ca": 1.76,
    "Sc": 1.7,
    "Ti": 1.6,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.61,  # low spin = 1.39
    "Fe": 1.52,  # low spin = 1.32
    "Co": 1.5,  # low spin = 1.26
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Ga": 1.22,
    "Ge": 1.2,
    "As": 1.19,
    "Se": 1.2,
    "Br": 1.2,
    "Kr": 1.16,
    "Rb": 2.2,
    "Sr": 1.95,
    "Y": 1.9,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Tc": 1.47,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "In": 1.42,
    "Sn": 1.39,
    "Sb": 1.39,
    "Te": 1.38,
    "I": 1.39,
    "Xe": 1.4,
    "Cs": 2.44,
    "Ba": 2.15,
    "La": 2.07,
    "Ce": 2.04,
    "Pr": 2.03,
    "Nd": 2.01,
    "Pm": 1.99,
    "Sm": 1.98,
    "Eu": 1.98,
    "Gd": 1.96,
    "Tb": 1.94,
    "Dy": 1.92,
    "Ho": 1.92,
    "Er": 1.89,
    "Tm": 1.9,
    "Yb": 1.87,
    "Lu": 1.87,
    "Hf": 1.75,
    "Ta": 1.7,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
    "Tl": 1.45,
    "Pb": 1.46,
    "Bi": 1.48,
    "Po": 1.4,
    "At": 1.5,
    "Rn": 1.5,
    "Fr": 2.6,
    "Ra": 2.21,
    "Ac": 2.15,
    "Th": 2.06,
    "Pa": 2,
    "U": 1.96,
    "Np": 1.9,
    "Pu": 1.87,
    "Am": 1.8,
    "Cm": 1.69,
}


UFF_FULL = {
    # Atom, r1, theta0, x1, D1, zeta, Z1, Vi, Uj, Xi, Hard, Radius
    "Du": (0.01, 180, 0.4, 5000, 12, 10, 0, 0, 9.66, 14.92, 0.7),
    "H_": (0.354, 180, 2.886, 0.044, 12, 0.712, 0, 0, 4.528, 6.9452, 0.371),
    "H_b": (0.46, 83.5, 2.886, 0.044, 12, 0.712, 0, 0, 4.528, 6.9452, 0.371),
    "He4+4": (0.849, 90, 2.362, 0.056, 15.24, 0.098, 0, 0, 9.66, 14.92, 1.3),
    "Li": (1.336, 180, 2.451, 0.025, 12, 1.026, 0, 2, 3.006, 2.386, 1.557),
    "Be3+2": (1.074, 109.47, 2.745, 0.085, 12, 1.565, 0, 2, 4.877, 4.443, 1.24),
    "B_3": (0.838, 109.47, 4.083, 0.18, 12.052, 1.755, 0, 2, 5.11, 4.75, 0.822),
    "B_2": (0.828, 120, 4.083, 0.18, 12.052, 1.755, 0, 2, 5.11, 4.75, 0.822),
    "C_3": (0.757, 109.47, 3.851, 0.105, 12.73, 1.912, 2.119, 2, 5.343, 5.063, 0.759),
    "C_R": (0.729, 120, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759),
    "C_2": (0.732, 120, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759),
    "C_1": (0.706, 180, 3.851, 0.105, 12.73, 1.912, 0, 2, 5.343, 5.063, 0.759),
    "N_3": (0.7, 106.7, 3.66, 0.069, 13.407, 2.544, 0.45, 2, 6.899, 5.88, 0.715),
    "N_R": (0.699, 120, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715),
    "N_2": (0.685, 111.2, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715),
    "N_1": (0.656, 180, 3.66, 0.069, 13.407, 2.544, 0, 2, 6.899, 5.88, 0.715),
    "O_3": (0.658, 104.51, 3.5, 0.06, 14.085, 2.3, 0.018, 2, 8.741, 6.682, 0.669),
    "O_3_z": (0.528, 146, 3.5, 0.06, 14.085, 2.3, 0.018, 2, 8.741, 6.682, 0.669),
    "O_3_M": (0.658, 109.47, 3.5, 0.06, 14.085, 2.3, 0.018, 2, 8.741, 6.682, 0.669),
    "O_R": (0.68, 110, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669),
    "O_2": (0.634, 120, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669),
    "O_1": (0.639, 180, 3.5, 0.06, 14.085, 2.3, 0, 2, 8.741, 6.682, 0.669),
    "F_": (0.668, 180, 3.364, 0.05, 14.762, 1.735, 0, 2, 10.874, 7.474, 0.706),
    "Ne4+4": (0.92, 90, 3.243, 0.042, 15.44, 0.194, 0, 2, 11.04, 10.55, 1.768),
    "Na": (1.539, 180, 2.983, 0.03, 12, 1.081, 0, 1.25, 2.843, 2.296, 2.085),
    "Mg3+2": (1.421, 109.47, 3.021, 0.111, 12, 1.787, 0, 1.25, 3.951, 3.693, 1.5),
    "Mg6": (1.421, 90, 3.021, 0.111, 12, 1.787, 0, 1.25, 3.951, 3.693, 1.5),
    "Al3": (1.244, 109.47, 4.499, 0.505, 11.278, 1.792, 0, 1.25, 4.06, 3.59, 1.201),
    "Si3": (1.117, 109.47, 4.295, 0.402, 12.175, 2.323, 1.225, 1.25, 4.168, 3.487, 1.176),
    "P_3+3": (1.101, 93.8, 4.147, 0.305, 13.072, 2.863, 2.4, 1.25, 5.463, 4, 1.102),
    "P_3+5": (1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.4, 1.25, 5.463, 4, 1.102),
    "P_3+q": (1.056, 109.47, 4.147, 0.305, 13.072, 2.863, 2.4, 1.25, 5.463, 4, 1.102),
    "S_3+2": (1.064, 92.1, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928, 4.486, 1.047),
    "S_3+4": (1.049, 103.2, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928, 4.486, 1.047),
    "S_3+6": (1.027, 109.47, 4.035, 0.274, 13.969, 2.703, 0.484, 1.25, 6.928, 4.486, 1.047),
    "S_R": (1.077, 92.2, 4.035, 0.274, 13.969, 2.703, 0, 1.25, 6.928, 4.486, 1.047),
    "S_2": (0.854, 120, 4.035, 0.274, 13.969, 2.703, 0, 1.25, 6.928, 4.486, 1.047),
    "Cl": (1.044, 180, 3.947, 0.227, 14.866, 2.348, 0, 1.25, 8.564, 4.946, 0.994),
    "Ar4+4": (1.032, 90, 3.868, 0.185, 15.763, 0.3, 0, 1.25, 9.465, 6.355, 2.108),
    "K_": (1.953, 180, 3.812, 0.035, 12, 1.165, 0, 0.7, 2.421, 1.92, 2.586),
    "Ca6+2": (1.761, 90, 3.399, 0.238, 12, 2.141, 0, 0.7, 3.231, 2.88, 2),
    "Sc3+3": (1.513, 109.47, 3.295, 0.019, 12, 2.592, 0, 0.7, 3.395, 3.08, 1.75),
    "Ti3+4": (1.412, 109.47, 3.175, 0.017, 12, 2.659, 0, 0.7, 3.47, 3.38, 1.607),
    "Ti6+4": (1.412, 90, 3.175, 0.017, 12, 2.659, 0, 0.7, 3.47, 3.38, 1.607),
    "V_3+5": (1.402, 109.47, 3.144, 0.016, 12, 2.679, 0, 0.7, 3.65, 3.41, 1.47),
    "Cr6+3": (1.345, 90, 3.023, 0.015, 12, 2.463, 0, 0.7, 3.415, 3.865, 1.402),
    "Mn6+2": (1.382, 90, 2.961, 0.013, 12, 2.43, 0, 0.7, 3.325, 4.105, 1.533),
    "Fe3+2": (1.27, 109.47, 2.912, 0.013, 12, 2.43, 0, 0.7, 3.76, 4.14, 1.393),
    "Fe6+2": (1.335, 90, 2.912, 0.013, 12, 2.43, 0, 0.7, 3.76, 4.14, 1.393),
    "Co6+3": (1.241, 90, 2.872, 0.014, 12, 2.43, 0, 0.7, 4.105, 4.175, 1.406),
    "Ni4+2": (1.164, 90, 2.834, 0.015, 12, 2.43, 0, 0.7, 4.465, 4.205, 1.398),
    "Cu3+1": (1.302, 109.47, 3.495, 0.005, 12, 1.756, 0, 0.7, 4.2, 4.22, 1.434),
    "Zn3+2": (1.193, 109.47, 2.763, 0.124, 12, 1.308, 0, 0.7, 5.106, 4.285, 1.4),
    "Ga3+3": (1.26, 109.47, 4.383, 0.415, 11, 1.821, 0, 0.7, 3.641, 3.16, 1.211),
    "Ge3": (1.197, 109.47, 4.28, 0.379, 12, 2.789, 0.701, 0.7, 4.051, 3.438, 1.189),
    "As3+3": (1.211, 92.1, 4.23, 0.309, 13, 2.864, 1.5, 0.7, 5.188, 3.809, 1.204),
    "Se3+2": (1.19, 90.6, 4.205, 0.291, 14, 2.764, 0.335, 0.7, 6.428, 4.131, 1.224),
    "Br": (1.192, 180, 4.189, 0.251, 15, 2.519, 0, 0.7, 7.79, 4.425, 1.141),
    "Kr4+4": (1.147, 90, 4.141, 0.22, 16, 0.452, 0, 0.7, 8.505, 5.715, 2.27),
    "Rb": (2.26, 180, 4.114, 0.04, 12, 1.592, 0, 0.2, 2.331, 1.846, 2.77),
    "Sr6+2": (2.052, 90, 3.641, 0.235, 12, 2.449, 0, 0.2, 3.024, 2.44, 2.415),
    "Y_3+3": (1.698, 109.47, 3.345, 0.072, 12, 3.257, 0, 0.2, 3.83, 2.81, 1.998),
    "Zr3+4": (1.564, 109.47, 3.124, 0.069, 12, 3.667, 0, 0.2, 3.4, 3.55, 1.758),
    "Nb3+5": (1.473, 109.47, 3.165, 0.059, 12, 3.618, 0, 0.2, 3.55, 3.38, 1.603),
    "Mo6+6": (1.467, 90, 3.052, 0.056, 12, 3.4, 0, 0.2, 3.465, 3.755, 1.53),
    "Mo3+6": (1.484, 109.47, 3.052, 0.056, 12, 3.4, 0, 0.2, 3.465, 3.755, 1.53),
    "Tc6+5": (1.322, 90, 2.998, 0.048, 12, 3.4, 0, 0.2, 3.29, 3.99, 1.5),
    "Ru6+2": (1.478, 90, 2.963, 0.056, 12, 3.4, 0, 0.2, 3.575, 4.015, 1.5),
    "Rh6+3": (1.332, 90, 2.929, 0.053, 12, 3.5, 0, 0.2, 3.975, 4.005, 1.509),
    "Pd4+2": (1.338, 90, 2.899, 0.048, 12, 3.21, 0, 0.2, 4.32, 4, 1.544),
    "Ag1+1": (1.386, 180, 3.148, 0.036, 12, 1.956, 0, 0.2, 4.436, 3.134, 1.622),
    "Cd3+2": (1.403, 109.47, 2.848, 0.228, 12, 1.65, 0, 0.2, 5.034, 3.957, 1.6),
    "In3+3": (1.459, 109.47, 4.463, 0.599, 11, 2.07, 0, 0.2, 3.506, 2.896, 1.404),
    "Sn3": (1.398, 109.47, 4.392, 0.567, 12, 2.961, 0.199, 0.2, 3.987, 3.124, 1.354),
    "Sb3+3": (1.407, 91.6, 4.42, 0.449, 13, 2.704, 1.1, 0.2, 4.899, 3.342, 1.404),
    "Te3+2": (1.386, 90.25, 4.47, 0.398, 14, 2.882, 0.3, 0.2, 5.816, 3.526, 1.38),
    "I_": (1.382, 180, 4.5, 0.339, 15, 2.65, 0, 0.2, 6.822, 3.762, 1.333),
    "Xe4+4": (1.267, 90, 4.404, 0.332, 12, 0.556, 0, 0.2, 7.595, 4.975, 2.459),
    "Cs": (2.57, 180, 4.517, 0.045, 12, 1.573, 0, 0.1, 2.183, 1.711, 2.984),
    "Ba6+2": (2.277, 90, 3.703, 0.364, 12, 2.727, 0, 0.1, 2.814, 2.396, 2.442),
    "La3+3": (1.943, 109.47, 3.522, 0.017, 12, 3.3, 0, 0.1, 2.8355, 2.7415, 2.071),
    "Ce6+3": (1.841, 90, 3.556, 0.013, 12, 3.3, 0, 0.1, 2.774, 2.692, 1.925),
    "Pr6+3": (1.823, 90, 3.606, 0.01, 12, 3.3, 0, 0.1, 2.858, 2.564, 2.007),
    "Nd6+3": (1.816, 90, 3.575, 0.01, 12, 3.3, 0, 0.1, 2.8685, 2.6205, 2.007),
    "Pm6+3": (1.801, 90, 3.547, 0.009, 12, 3.3, 0, 0.1, 2.881, 2.673, 2),
    "Sm6+3": (1.78, 90, 3.52, 0.008, 12, 3.3, 0, 0.1, 2.9115, 2.7195, 1.978),
    "Eu6+3": (1.771, 90, 3.493, 0.008, 12, 3.3, 0, 0.1, 2.8785, 2.7875, 2.227),
    "Gd6+3": (1.735, 90, 3.368, 0.009, 12, 3.3, 0, 0.1, 3.1665, 2.9745, 1.968),
    "Tb6+3": (1.732, 90, 3.451, 0.007, 12, 3.3, 0, 0.1, 3.018, 2.834, 1.954),
    "Dy6+3": (1.71, 90, 3.428, 0.007, 12, 3.3, 0, 0.1, 3.0555, 2.8715, 1.934),
    "Ho6+3": (1.696, 90, 3.409, 0.007, 12, 3.416, 0, 0.1, 3.127, 2.891, 1.925),
    "Er6+3": (1.673, 90, 3.391, 0.007, 12, 3.3, 0, 0.1, 3.1865, 2.9145, 1.915),
    "Tm6+3": (1.66, 90, 3.374, 0.006, 12, 3.3, 0, 0.1, 3.2514, 2.9329, 2),
    "Yb6+3": (1.637, 90, 3.355, 0.228, 12, 2.618, 0, 0.1, 3.2889, 2.965, 2.158),
    "Lu6+3": (1.671, 90, 3.64, 0.041, 12, 3.271, 0, 0.1, 2.9629, 2.4629, 1.896),
    "Hf3+4": (1.611, 109.47, 3.141, 0.072, 12, 3.921, 0, 0.1, 3.7, 3.4, 1.759),
    "Ta3+5": (1.511, 109.47, 3.17, 0.081, 12, 4.075, 0, 0.1, 5.1, 2.85, 1.605),
    "W_6+6": (1.392, 90, 3.069, 0.067, 12, 3.7, 0, 0.1, 4.63, 3.31, 1.538),
    "W_3+4": (1.526, 109.47, 3.069, 0.067, 12, 3.7, 0, 0.1, 4.63, 3.31, 1.538),
    "W_3+6": (1.38, 109.47, 3.069, 0.067, 12, 3.7, 0, 0.1, 4.63, 3.31, 1.538),
    "Re6+5": (1.372, 90, 2.954, 0.066, 12, 3.7, 0, 0.1, 3.96, 3.92, 1.6),
    "Re3+7": (1.314, 109.47, 2.954, 0.066, 12, 3.7, 0, 0.1, 3.96, 3.92, 1.6),
    "Os6+6": (1.372, 90, 3.12, 0.037, 12, 3.7, 0, 0.1, 5.14, 3.63, 1.7),
    "Ir6+3": (1.371, 90, 2.84, 0.073, 12, 3.731, 0, 0.1, 5, 4, 1.866),
    "Pt4+2": (1.364, 90, 2.754, 0.08, 12, 3.382, 0, 0.1, 4.79, 4.43, 1.557),
    "Au4+3": (1.262, 90, 3.293, 0.039, 12, 2.625, 0, 0.1, 4.894, 2.586, 1.618),
    "Hg1+2": (1.34, 180, 2.705, 0.385, 12, 1.75, 0, 0.1, 6.27, 4.16, 1.6),
    "Tl3+3": (1.518, 120, 4.347, 0.68, 11, 2.068, 0, 0.1, 3.2, 2.9, 1.53),
    "Pb3": (1.459, 109.47, 4.297, 0.663, 12, 2.846, 0.1, 0.1, 3.9, 3.53, 1.444),
    "Bi3+3": (1.512, 90, 4.37, 0.518, 13, 2.47, 1, 0.1, 4.69, 3.74, 1.514),
    "Po3+2": (1.5, 90, 4.709, 0.325, 14, 2.33, 0.3, 0.1, 4.21, 4.21, 1.48),
    "At": (1.545, 180, 4.75, 0.284, 15, 2.24, 0, 0.1, 4.75, 4.75, 1.47),
    "Rn4+4": (1.42, 90, 4.765, 0.248, 16, 0.583, 0, 0.1, 5.37, 5.37, 2.2),
    "Fr": (2.88, 180, 4.9, 0.05, 12, 1.847, 0, 0, 2, 2, 2.3),
    "Ra6+2": (2.512, 90, 3.677, 0.404, 12, 2.92, 0, 0, 2.843, 2.434, 2.2),
    "Ac6+3": (1.983, 90, 3.478, 0.033, 12, 3.9, 0, 0, 2.835, 2.835, 2.108),
    "Th6+4": (1.721, 90, 3.396, 0.026, 12, 4.202, 0, 0, 3.175, 2.905, 2.018),
    "Pa6+4": (1.711, 90, 3.424, 0.022, 12, 3.9, 0, 0, 2.985, 2.905, 1.8),
    "U_6+4": (1.684, 90, 3.395, 0.022, 12, 3.9, 0, 0, 3.341, 2.853, 1.713),
    "Np6+4": (1.666, 90, 3.424, 0.019, 12, 3.9, 0, 0, 3.549, 2.717, 1.8),
    "Pu6+4": (1.657, 90, 3.424, 0.016, 12, 3.9, 0, 0, 3.243, 2.819, 1.84),
    "Am6+4": (1.66, 90, 3.381, 0.014, 12, 3.9, 0, 0, 2.9895, 3.0035, 1.942),
    "Cm6+3": (1.801, 90, 3.326, 0.013, 12, 3.9, 0, 0, 2.8315, 3.1895, 1.9),
    "Bk6+3": (1.761, 90, 3.339, 0.013, 12, 3.9, 0, 0, 3.1935, 3.0355, 1.9),
    "Cf6+3": (1.75, 90, 3.313, 0.013, 12, 3.9, 0, 0, 3.197, 3.101, 1.9),
    "Es6+3": (1.724, 90, 3.299, 0.012, 12, 3.9, 0, 0, 3.333, 3.089, 1.9),
    "Fm6+3": (1.712, 90, 3.286, 0.012, 12, 3.9, 0, 0, 3.4, 3.1, 1.9),
    "Md6+3": (1.689, 90, 3.274, 0.011, 12, 3.9, 0, 0, 3.47, 3.11, 1.9),
    "No6+3": (1.679, 90, 3.248, 0.011, 12, 3.9, 0, 0, 3.475, 3.175, 1.9),
    "Lw6+3": (1.698, 90, 3.236, 0.011, 12, 3.9, 0, 0, 3.5, 3.2, 1.9)
}
