"""
Pen-Robinson equation of state for faps.

"""

from math import log, exp

import numpy as np


# Universal gas constant
R = 8.314

# From the excel implementation at http://www.egr.msu.edu/~lira/readcomp.htm
# Converted MPa -> to -> bar (multiply by 10)
properties = {
    '1,3-butadiene': {'t_crit': 425.4, 'p_crit': 43.3, 'omega': 0.193},
    '1-butanol': {'t_crit': 562.9, 'p_crit': 44.12, 'omega': 0.594},
    '1-butene': {'t_crit': 419.6, 'p_crit': 40.2, 'omega': 0.187},
    '1-methylnaphthalene': {'t_crit': 772, 'p_crit': 36.5, 'omega': 0.292},
    '1-pentene': {'t_crit': 464.8, 'p_crit': 35.29, 'omega': 0.233},
    'acetic acid': {'t_crit': 592.7, 'p_crit': 57.86, 'omega': 0.462},
    'acetone': {'t_crit': 508.2, 'p_crit': 47.01, 'omega': 0.306},
    'acetonitrile': {'t_crit': 545.5, 'p_crit': 48.33, 'omega': 0.353},
    'acetylene': {'t_crit': 308.3, 'p_crit': 61.39, 'omega': 0.187},
    'ammonia': {'t_crit': 406.6, 'p_crit': 112.7, 'omega': 0.252},
    'argon': {'t_crit': 150.9, 'p_crit': 48.98, 'omega': -0.004},
    'benzene': {'t_crit': 562.2, 'p_crit': 48.98, 'omega': 0.211},
    'biphenyl': {'t_crit': 789.3, 'p_crit': 38.47, 'omega': 0.366},
    'bromine': {'t_crit': 584.2, 'p_crit': 103.35, 'omega': 0.119},
    'carbon dioxide': {'t_crit': 304.2, 'p_crit': 73.82, 'omega': 0.228},
    'carbon disulfide': {'t_crit': 552, 'p_crit': 78, 'omega': 0.115},
    'carbon monoxide': {'t_crit': 132.9, 'p_crit': 34.99, 'omega': 0.066},
    'carbon tetrachloride': {'t_crit': 556.4, 'p_crit': 45, 'omega': 0.194},
    'chlorine': {'t_crit': 417.2, 'p_crit': 77.11, 'omega': 0.069},
    'chlorobenzene': {'t_crit': 632.4, 'p_crit': 44.6, 'omega': 0.249},
    'chloroform': {'t_crit': 536.4, 'p_crit': 54, 'omega': 0.216},
    'cumene': {'t_crit': 631.2, 'p_crit': 32.09, 'omega': 0.338},
    'cyclohexane': {'t_crit': 553.5, 'p_crit': 40.75, 'omega': 0.215},
    'cyclopentane': {'t_crit': 511.8, 'p_crit': 45.02, 'omega': 0.194},
    'diethyl ether': {'t_crit': 466.7, 'p_crit': 35.9, 'omega': 0.281},
    'diphenylmethane': {'t_crit': 768, 'p_crit': 29.2, 'omega': 0.461},
    'ethane': {'t_crit': 305.4, 'p_crit': 48.8, 'omega': 0.099},
    'ethanol': {'t_crit': 516.4, 'p_crit': 63.84, 'omega': 0.637},
    'ethylbenzene': {'t_crit': 617.2, 'p_crit': 36.09, 'omega': 0.304},
    'ethylene oxide': {'t_crit': 469, 'p_crit': 71, 'omega': 0.2},
    'ethylene': {'t_crit': 282.4, 'p_crit': 50.32, 'omega': 0.085},
    'freon-11': {'t_crit': 471.2, 'p_crit': 43.5, 'omega': 0.188},
    'freon-113': {'t_crit': 487.2, 'p_crit': 33.7, 'omega': 0.252},
    'freon-12': {'t_crit': 385, 'p_crit': 40.7, 'omega': 0.176},
    'freon-22': {'t_crit': 369.8, 'p_crit': 49.7, 'omega': 0.221},
    'helium-4': {'t_crit': 5.2, 'p_crit': 2.28, 'omega': 0},
    'hydrazine': {'t_crit': 653, 'p_crit': 145, 'omega': 0.328},
    'hydrogen chloride': {'t_crit': 324.6, 'p_crit': 82, 'omega': 0.12},
    'hydrogen cyanide': {'t_crit': 456.8, 'p_crit': 53.2, 'omega': 0.407},
    'hydrogen sulfide': {'t_crit': 373.5, 'p_crit': 89.37, 'omega': 0.081},
    'hydrogen': {'t_crit': 33.3, 'p_crit': 12.97, 'omega': -0.215},
    'isobutane': {'t_crit': 408.1, 'p_crit': 36.48, 'omega': 0.177},
    'isobutanol': {'t_crit': 547.7, 'p_crit': 42.95, 'omega': 0.589},
    'isobutene': {'t_crit': 417.9, 'p_crit': 39.99, 'omega': 0.189},
    'isopentane': {'t_crit': 460.4, 'p_crit': 33.81, 'omega': 0.228},
    'isoprene': {'t_crit': 484, 'p_crit': 38.5, 'omega': 0.158},
    'isopropanol': {'t_crit': 508.3, 'p_crit': 47.64, 'omega': 0.669},
    'krypton': {'t_crit': 209.4, 'p_crit': 55.02, 'omega': 0.001},
    'm-xylene': {'t_crit': 617.1, 'p_crit': 35.41, 'omega': 0.326},
    'methane': {'t_crit': 190.6, 'p_crit': 46.04, 'omega': 0.011},
    'methanol': {'t_crit': 512.6, 'p_crit': 80.96, 'omega': 0.566},
    'methyl chloride': {'t_crit': 416.3, 'p_crit': 65.9, 'omega': 0.156},
    'methyl ethyl ketone': {'t_crit': 535.6, 'p_crit': 41, 'omega': 0.329},
    'methylcyclohexane': {'t_crit': 572.2, 'p_crit': 34.71, 'omega': 0.235},
    'methylcyclopentane': {'t_crit': 532.8, 'p_crit': 37.85, 'omega': 0.23},
    'n-butane': {'t_crit': 425.2, 'p_crit': 37.97, 'omega': 0.193},
    'n-decane': {'t_crit': 618.5, 'p_crit': 21.23, 'omega': 0.484},
    'n-dodecane': {'t_crit': 658.2, 'p_crit': 18.24, 'omega': 0.575},
    'n-heptane': {'t_crit': 540.3, 'p_crit': 27.36, 'omega': 0.349},
    'n-hexadecane': {'t_crit': 720.6, 'p_crit': 14.19, 'omega': 0.747},
    'n-hexane': {'t_crit': 507.4, 'p_crit': 30.12, 'omega': 0.305},
    'n-nonane': {'t_crit': 595.7, 'p_crit': 23.06, 'omega': 0.437},
    'n-octane': {'t_crit': 568.8, 'p_crit': 24.86, 'omega': 0.396},
    'n-pentane': {'t_crit': 469.7, 'p_crit': 33.69, 'omega': 0.249},
    'n-tetradecane': {'t_crit': 696.9, 'p_crit': 14.38, 'omega': 0.57},
    'naphthalene': {'t_crit': 748.4, 'p_crit': 40.51, 'omega': 0.302},
    'neon': {'t_crit': 44.4, 'p_crit': 26.53, 'omega': -0.041},
    'neopentane': {'t_crit': 433.8, 'p_crit': 31.99, 'omega': 0.196},
    'nitric oxide': {'t_crit': 180.2, 'p_crit': 64.85, 'omega': 0.585},
    'nitrogen': {'t_crit': 126.1, 'p_crit': 33.94, 'omega': 0.04},
    'nitrous oxide': {'t_crit': 309.6, 'p_crit': 72.45, 'omega': 0.142},
    'o-xylene': {'t_crit': 630.4, 'p_crit': 37.34, 'omega': 0.313},
    'oxygen': {'t_crit': 154.6, 'p_crit': 50.43, 'omega': 0.022},
    'p-xylene': {'t_crit': 616.3, 'p_crit': 35.11, 'omega': 0.326},
    'propane': {'t_crit': 369.8, 'p_crit': 42.49, 'omega': 0.152},
    'propanol': {'t_crit': 536.7, 'p_crit': 51.7, 'omega': 0.628},
    'propylene': {'t_crit': 364.8, 'p_crit': 46.13, 'omega': 0.142},
    'sulfur dioxide': {'t_crit': 430.8, 'p_crit': 78.84, 'omega': 0.245},
    'sulfur trioxide': {'t_crit': 490.9, 'p_crit': 82.07, 'omega': 0.422},
    'tetralin': {'t_crit': 720.2, 'p_crit': 33, 'omega': 0.286},
    'thf': {'t_crit': 501.1, 'p_crit': 51.9, 'omega': 0.217},
    'toluene': {'t_crit': 591.8, 'p_crit': 41.09, 'omega': 0.264},
    'water': {'t_crit': 647.3, 'p_crit': 221.2, 'omega': 0.344},
    'xenon': {'t_crit': 289.7, 'p_crit': 58.4, 'omega': 0.012}
}

# Binary mixing rules taken from
# http://dx.doi.org/10.1016/j.jare.2012.03.004.
k12 = {'benzene': {'benzene': 0.0, 'carbon dioxide': 0.0774, 'ethane': 0.0322,
                   'heptane': 0.0011, 'hexane': 0.0093,
                   'hydrogen sulfide': 0.00293, 'methane': 0.0363,
                   'nitrogen': 0.1641},
       'butane': {'butane': 0.0, 'hydrogen sulfide': 0.11554, 'nitrogen': 0.08},
       'carbon dioxide': {'benzene': 0.0774, 'carbon dioxide': 0.0,
                          'decane': 0.1141, 'ethane': 0.1322, 'heptane': 0.1,
                          'i-butane': 0.12, 'i-pentane': 0.1219,
                          'm-xylene': 0.14339, 'methane': 0.0919,
                          'n-butane': 0.1333, 'n-hexane': 0.11,
                          'n-pentane': 0.1222, 'nitrogen': -0.017,
                          'octane': 0.13303, 'propane': 0.1241,
                          'toluene': 0.1056},
       'decane': {'carbon dioxide': 0.1141, 'decane': 0.0,
                  'hydrogen sulfide': 0.0333},
       'ethane': {'benzene': 0.0322, 'carbon dioxide': 0.1322, 'ethane': 0.0,
                  'heptane': 0.0067, 'hexane': -0.01,
                  'hydrogen sulfide': 0.0833, 'i-butane': -0.0067,
                  'methane': -0.0026, 'n-butane': 0.0096, 'nitrogen': 0.0515,
                  'octane': 0.0185, 'propane': 0.0011},
       'heptane': {'benzene': 0.0011, 'carbon dioxide': 0.1, 'ethane': 0.0067,
                   'heptane': 0.0, 'hydrogen sulfide': 0.06164,
                   'methane': 0.0352, 'nitrogen': 0.1441},
       'hexane': {'benzene': 0.0093, 'ethane': -0.01, 'hexane': 0.0,
                  'hydrogen sulfide': 0.05744, 'methane': 0.0422,
                  'nitrogen': 0.1496},
       'hydrogen sulfide': {'benzene': 0.00293, 'butane': 0.11554,
                            'decane': 0.0333, 'ethane': 0.0833,
                            'heptane': 0.06164, 'hexane': 0.05744,
                            'hydrogen sulfide': 0.0, 'i-butane': 0.0474,
                            'm-xylene': 0.0172, 'methane': 0.08857,
                            'nitrogen': 0.1767, 'pentane': 0.063,
                            'toluene': 0.00751},
       'i-butane': {'carbon dioxide': 0.12, 'ethane': -0.0067,
                    'hydrogen sulfide': 0.0474, 'i-butane': 0.0,
                    'methane': 0.0256, 'propane': -0.0078},
       'i-pentane': {'carbon dioxide': 0.1219, 'i-pentane': 0.0,
                     'propane': 0.0111},
       'm-xylene': {'carbon dioxide': 0.14339, 'hydrogen sulfide': 0.0172,
                    'm-xylene': 0.0, 'methane': 0.0844},
       'methane': {'benzene': 0.0363, 'carbon dioxide': 0.0919,
                   'ethane': -0.0026, 'heptane': 0.0352, 'hexane': 0.0422,
                   'hydrogen sulfide': 0.08857, 'i-butane': 0.0256,
                   'm-xylene': 0.0844, 'methane': 0.0, 'n-butane': 0.0133,
                   'n-decane': 0.0422, 'n-pentane': 0.023, 'nitrogen': 0.0311,
                   'nonane': 0.0474, 'propane': 0.014, 'toluene': 0.097},
       'n-butane': {'carbon dioxide': 0.1333, 'ethane': 0.0096,
                    'methane': 0.0133, 'n-butane': 0.0},
       'n-decane': {'methane': 0.0422, 'n-decane': 0.0},
       'n-hexane': {'carbon dioxide': 0.11, 'n-hexane': 0.0},
       'n-pentane': {'carbon dioxide': 0.1222, 'methane': 0.023,
                     'n-pentane': 0.0},
       'nitrogen': {'benzene': 0.1641, 'butane': 0.08, 'carbon dioxide': -0.017,
                    'ethane': 0.0515, 'heptane': 0.1441, 'hexane': 0.1496,
                    'hydrogen sulfide': 0.1767, 'methane': 0.0311,
                    'nitrogen': 0.0, 'octane': -0.41, 'pentane': 0.1,
                    'propane': 0.0852, 'toluene': 0.20142},
       'nonane': {'methane': 0.0474, 'nonane': 0.0},
       'octane': {'carbon dioxide': 0.13303, 'ethane': 0.0185,
                  'nitrogen': -0.41, 'octane': 0.0},
       'pentane': {'hydrogen sulfide': 0.063, 'nitrogen': 0.1, 'pentane': 0.0,
                   'toluene': 0.00845},
       'propane': {'carbon dioxide': 0.1241, 'ethane': 0.0011,
                   'i-butane': -0.0078, 'i-pentane': 0.0111, 'methane': 0.014,
                   'nitrogen': 0.0852, 'propane': 0.0},
       'toluene': {'carbon dioxide': 0.1056, 'hydrogen sulfide': 0.00751,
                   'methane': 0.097, 'nitrogen': 0.20142, 'pentane': 0.00845,
                   'toluene': 0.0}}


def peng_robinson(pressures, temperature=298.0):
    """
    Calculate the fugacity of gasses based on the Peng Robinson equation of
    state.

    :param: pressure is a dictionary of {"chemical species": pressure in bar}
    :param: temperature is in Kelvin

    :return: a dictionary of {"chemical species": fugacity in bar}

    :raises: KeyError if any species (or mixing) is unknown

    Based in part on https://www.e-education.psu.edu/png520/m11_p2.html
    With Tina Duren's implementation used for validation
    """

    kappa = {}
    alpha_t = {}
    agreek = {}
    bgreek = {}

    # Check if we have values for the parameters. Insert individual species
    # into the mixing rules if necessary with 0.0 for mixing.
    for species_i in pressures:
        if not species_i in properties:
            raise KeyError("Species %s unknown." % species_i)
        elif not species_i in k12:
            if len(pressures) == 1:
                # self mixing is all we need
                k12[species_i] = {species_i: 0.0}
        for species_j in pressures:
            if not species_j in k12[species_i]:
                raise KeyError("No mixing rule for %s and %s" %
                               (species_i, species_j))

    for species in pressures:
        kappa[species] = (0.37464 + 1.54226*properties[species]['omega'] -
                          0.26992*properties[species]['omega']**2)
        alpha_t[species] = (1 + kappa[species] *
                           (1 - (temperature/properties[species]['t_crit'])**0.5))**2
        agreek[species] = (0.45724 * (R * properties[species]['t_crit'])**2 /
                           properties[species]['p_crit'] * alpha_t[species])
        bgreek[species] = (0.0778 * R * properties[species]['t_crit'] /
                           properties[species]['p_crit'])

    # Calculate a and b factors for entire system
    total_pressure = sum(pressures.values())

    agreek_sys = 0.0
    bgreek_sys = 0.0

    if len(pressures) == 1:
        agreek_sys = agreek.values()[0]
        bgreek_sys = bgreek.values()[0]
    else:
        for species_i in pressures:
            for species_j in pressures:
                agreek_sys += (
                    (pressures[species_i]/total_pressure) *
                    (pressures[species_j]/total_pressure) *
                    ((agreek[species_i]*agreek[species_j])**0.5) *
                    (1 - k12[species_i][species_j]))
            bgreek_sys += ((pressures[species_i]/total_pressure) *
                           bgreek[species_i])

    # calculate A and B
    A = agreek_sys * total_pressure / (R * temperature) ** 2
    B = bgreek_sys * total_pressure / R / temperature

    # calculate alpha, beta, gamma
    alpha = B - 1
    beta_f = A - 3 * B * B - 2 * B
    gamma = -A * B + B * B + B * B * B

    Z1 = np.roots([1, alpha, beta_f, gamma])[0].real

    # calculate the sum(aik*zi)
    sum_za = dict((x, 0.0) for x in pressures)
    for species_i in pressures:
        for species_j in pressures:
            sum_za[species_i] += (
                (pressures[species_j]/total_pressure) *
                ((agreek[species_i]*agreek[species_j])**0.5) *
                (1 - k12[species_i][species_j]))

    # calculate fugacity coefficient
    fugacity_coefficient = {}
    fugacity = {}
    for species in pressures:
        fugacity_coefficient[species] = exp(
            bgreek[species] / bgreek_sys * (Z1 - 1) - log(Z1 - B) -
            agreek_sys / (2.0 * (2.0 ** 0.5) * bgreek_sys * R * temperature) *
            (2 * sum_za[species] / agreek_sys - bgreek[species] / bgreek_sys)
            * log((Z1 + (1 + (2.0 ** 0.5)) * B) /
                  (Z1 + (1 - (2.0 ** 0.5)) * B)))

        fugacity[species] = fugacity_coefficient[species] * pressures[species]

    return fugacity
