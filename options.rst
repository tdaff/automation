================
Customising faps
================

Faps aims to provide sensible defaults for most options, however it also aims
to offer options for customising calculations. Defaults for all options can be
found in the ``defaults.ini``; site-wide options such as the exe names and
directories should be set in the :ref:`site.ini <site-ini>`.

Config files use a simple ini format of ``options = value``.



# Method for calculating charge [repeat, ]
charge_method = repeat
# Number of CPUs to run CPMD on [int]
cpmd_ncpu = 8
# Command to convert [str]
cpmd2cube = cpmd2cube.x
# Method to use for dft/optimization [vasp, ]
dft_code = vasp
# Turn on dispersion corrections in dft? [bool]
dispersion = True
# Resolution, in Angstrom, of the esp grid [float]
esp_resolution = 0.1
# Source for ESP (same as dft, usually) [vasp, ]
esp_src = vasp
# Location of fastmc executable
fastmc_exe = fastmc
# Number of CPUs to run fastmc on [int]
fastmc_ncpu = 1
# Guest(s) to use in GCMC [string] [see guests.lib]
guests = CO2
# Try to read in data from a previous calculation [bool]
import = False
# Filetype for input file [str]
initial_structure_format = pdb
# Specify the input filename [str]
input_structure_file = None
# Enable interactive interface to the system [bool]
interactive = False
# Kpoint grid size for dft calculations [int, int, int]
kpoints = 1, 1, 1
# Method to use for Monte Carlo [fastmc, ]
mc_code = fastmc
# Potential cutoff in GCMC (also affects supercell) [float]
mc_cutoff = 12.5
# GCMC equilibration steps [int]
mc_eq_steps = 1000000
# How often to write the fastmc history [int]
mc_history_freq = 1000
# Add the 'jobcontrol' directive for manual stop of GCMC [bool]
mc_jobcontrol = False
# How often to write the fastmc numguests [int]
mc_numguests_freq = 1000
# GCMC pressure (bar) [float]
mc_pressure = 1.0
# Turn on probability plots in GCMC? [bool]
mc_probability_plot = True
# GCMC production steps [int]
mc_prod_steps = 10000000
# Supercell to use for GCMC (cutoff determines minimum) [(int, int, int)]
mc_supercell = (1, 1, 1)
# Temperature to use in GCMC (K) [float]
mc_temperature = 273
# How often to write the fastmc numguests/his (deprecated) [int]
mc_write_freq = 1000
# Skip the charge calculation step? [bool]
no_charges = False
# Skip the dft/optimization step? [bool]
no_dft = False
# Skip the gcmc step? [bool]
no_gcmc = False
# Do not submit jobs; just create input files
no_submit = False
# Optimize positions of hydrogens in dft/optimization step [bool]
optim_h = True
# Optimize all atom positions in dft/optimization step [bool]
optim_all = False
# Optimize cell vectors in dft/optimization step [bool]
optim_cell = False
# Location of VASP POTCARs; each element in a folder [str]
potcar_dir = /home/program/CHEMISTRY/VASP/potpaw_PBE_2011/
# Location of siesta psf pseudopotentials
psf_dir = .
# Queuing system to use [wooki, ]
queue = wooki
# Silence stdout [ignored; set on commandline]
quiet = False
# Location of REPEAT executable [str]
repeat_exe = repeat.x
# Cpus to use for REPEAT calculation [int]
repeat_ncpu = 1
# Run all the steps without stopping [bool]
run_all = True
# Maximum memory that can be used for serial calculations in GB [float]
serial_memory = 2.5
# General acucracy setting for siesta calcualtions [low, med, high]
siesta_accuracy = med
# Location of siesta executable
siesta_exe = siesta
# Number of CPUs to use for siesta
siesta_ncpu = 1
# Command to convert siesta ESP to .cube for repeat [str]
siesta_to_cube = siesta2repeat
# Spin polarized dft? [bool]
spin = False
# Codes that run with openmp threads, not mpi
threaded_codes = repeat
# Maximum memory to use for threaded calculations in GB [float]
threaded_memory = 12
# Re-read options on restart? [bool]
update_opts = True
# Name (location) of vasp executable [str]
vasp_exe = vasp
# Number of cpus to run vasp on [int]
vasp_ncpu = 8
# Command to convert LOCPOT to .cube for REPEAT [str]
vasp_to_cube = vasp2cube 2
# Print debugging information [ignored; set on commandline]
verbose = False
