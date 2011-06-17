#! /usr/bin/env python

import re
import os
import sys

EXIT_FAILURE = 1
EXIST_SUCCESS = 0

value = ""
regex = ""
VALUES = {}

DEBUGGING = False
#DEBUGGING = True


job_input_file = open("JOB_INPUT","r")
input_list=job_input_file.readlines()
for i in range(len(input_list)):
    input_list[i] = input_list[i][:input_list[i].find("#")].replace("\n","")

for i in range(len(input_list)):
    input_list[i] = input_list[i].strip().split()

while True:
 try:
   input_list.remove([])
 except:
   break


for i in range(len(input_list)):
    if input_list[i][0][:1] == "&":
        if len(input_list[i]) == 2:
            VALUES[input_list[i][0][1:]] = input_list[i][1]
        elif len(input_list[i]) == 1:
            VALUES[input_list[i][0][1:]] = ""
    elif input_list[i][0][:1] != "&":
        print "-"*120+"\nWARNING!: "+input_list[i][0]+" does not seem to be a properly formatted variable name, check JOB_INPUT file.\n"+"-"*120
         
#------------------

DEFAULTS = {
"CPMD_NUM_PROCESSORS"            :["Number of Processors to use in CPMD",8,"^[1-9]{1}$|^[0-9]+$",int],
"CPMD_CONVERGENCE_CRIT_ORBITALS" :["Orbital Convergence criteria in CPMD","1.0d-5","^[0-9]*\.?[0-9]+d[-+]?[0-9]",str],
"CPMD_MAX_INTERATIONS"           :["Maximum number of optimization steps",10000,"^[0-9]+$",int],
"CPMD_PSP_PATH"                  :["Absolute path for pseudopotential files","/home/program/CHEMISTRY/CPMD/Pseudo_potentials/","^((?:\/[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)*(?:\-[a-zA-Z0-9]+)*)+)\/?$",str],
"CPMD_DEFAULT_FUNCTIONAL"        :["Default functional type","PBE","^[A-Za-z]+",str],
"CPMD_DEFAULT_PSP_TYPE"          :["Default pseudopotential type","SG","^[A-Za-z0-9]+$",str],
"CPMD_DEFAULT_LMAX_FLAG"         :["Default highest orbital consideration","P","^[S,P,D,F]{1}",str],
"CPMD_CUTOFF"                    :["Cutoff for CPMD calculation in angstroms",70.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"CPMD_CHARGE"                    :["Charge of system for CPMD calculation",0,"^[+-]?[0-9]+$",int],
"CPMD_OPTIMIZE_HYDROGENS"        :["Flag for whether or not to optimize Hydrogens","False","^True$|^False$",bool],
"DLPOLY_GUEST_GRIDPOINTS"        :["number of gridlines per dimension for guest insertion search",10.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"DLPOLY_MIN_GUEST_DISTANCE"      :["Minimum distance to framework for guest insertion in Angstroms.",1.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"DLPOLY_MOL_TYPES"               :["Number of DL_POLY Molecular types",2,"^[1-9]{1}$|^[0-9]{1}[1-9]{1}$",int],
"DLPOLY_SUPERCELL_DIM_A"         :["Multiplicity of supercell in Direction A",2,"^[1-9]{1}$|^[0-9]{1}[1-9]{1}$",int],
"DLPOLY_SUPERCELL_DIM_B"         :["Multiplicity of supercell in Direction B",2,"^[1-9]{1}$|^[0-9]{1}[1-9]{1}$",int],
"DLPOLY_SUPERCELL_DIM_C"         :["Multiplicity of supercell in Direction C",2,"^[1-9]{1}$|^[0-9]{1}[1-9]{1}$",int],
"DLPOLY_LEVCFG"                  :["See DL_POLY manual",0,"^[0-2]$",int],
"DLPOLY_IMCON"                   :["See DL_POLY manual",3,"^[0-7]$",int],
"DLPOLY_VDW_FILE"                :["Filename for VDW params","VDW_DLP","^\.{0,2}((?:\/?[a-zA-Z0-9]+(?:_[a-zA-Z0-9\.]+)*(?:\-[a-zA-Z0-9\.]+)*)+)\/?$",str],
"DLPOLY_GUEST_FILE"              :["Filename for guest information","GUEST","^\.{0,2}((?:\/?[a-zA-Z0-9]+(?:_[a-zA-Z0-9\.]+)*(?:\-[a-zA-Z0-9\.]+)*)+)\/?$",str],
"DL_POLY_CUTOFF"                 :["Required forces cutoff (Angstroms)",8.5,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"DLPOLY_DELR"                    :["Verlet neighbour list shell width (Angstroms)",1.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"DL_POLY_EWALD_PRECISION"        :["Ewald sum for electrostatics, with automatic parameter optimisation (0 < f < 1.E . 4)","1.0d-6","^[0-9]*\.?[0-9]+d[-+]?[0-9]",str],
"REPEAT_ESP_CUBEFILE"            :["Input ESP file name in cube format","REPEAT_ESP.cube","^[A-Za-z]{1}[A-Za-z0-9_]*[\.]?[A-Za-z0-9\_]*$",str],
"REPEAT_PERIODICITY"             :["Fit molecular(0) or periodic(1:default) system?",1,"^[0-1]$",int],
"REPEAT_VDW_FACTOR"              :["van der Waals scaling factor (default = 1.0)",1.0,"^[0-9]+[\.]{1}[0-9]+$",float],
"REPEAT_RESP_PENALTIES"          :["Apply RESP penalties?, no(0:default), yes(1)",0,"^[0-1]$",int],
"REPEAT_READ_CUTOFF"             :["Read cutoff radius? no(0), yes(1:default)",1,"^[0-1]$",int],
"REPEAT_R_CUTOFF"                :["If flag above=1 provide R_cutoff next (in Bohrs)",20.0,"^[0-9]+[\.]{1}[0-9]+$",float],
"REPEAT_APPLY_SYMMETRY_RESTRAIN" :["Apply symmetry restrain? no(0:default), yes(1)",0,"^[0-1]$",int],
"REPEAT_USE_GODDARD"             :["Use Goddard-type restrains? no(0:default), yes(1)",0,"^[0-1]$",int],
"REPEAT_GODDARD_WEIGHT"          :["If flag above=1 then provide weight next",0.00000,"^[0-9]+[\.]{1}[0-9]+$",float],
"REPEAT_CHARGE"                  :["Enter total charge of the system",0.00000,"^[0-9]+[\.]{1}[0-9]+$",float],
"GCMC_NUMBER_TO_AVERAGE"         :["See GCMC README",1000000,"^[1-9]{1}$|^[0-9]+$",int ],
"GCMC_ACCEPTANCES"               :["See GCMC README",1000000,"^[1-9]{1}$|^[0-9]+$",int ],
"GCMC_PRODUCTION_AVERAGE"        :["See GCMC README",1000000,"^[1-9]{1}$|^[0-9]+$",int ],
"GCMC_IDEAL_GAS_PRESSURE"        :["See GCMC README",20000,"^[1-9]{1}$|^[0-9]+$",int ], 
"GCMC_PROBABILITY_PLOT"          :["See GCMC README",100,"^[1-9]{1}$|^[0-9]+$",int ],
"GCMC_WINDOW_SIZE"               :["Default window size for convergence criteria",100000,"^[1-9]{1}$|^[0-9]+$",int],
"GCMC_CONVERGENCE_CRITERIA"      :["GCMC Convergence criteria",0.001,"^[0-9]+[\.]{1}[0-9]+$",float],
"GCMC_GCMCPAR"                   :["GCMC_GCMCPAR","none","^[A-Za-z1-9]+$",str],
"GCMC_NUM_PROCESSORS"            :["Number of branches for GCMC run",8,"^[1-9]{1}$|^[0-9]+$",int],
"GCMC_WAIT_TIME"                 :["Wait time between convergence assessments in minutes.",20,"^[1-9]{1}$|^[0-9]+$",int ],
"GCMC_GUEST_NAME"                :["Name of Guest to pull from library","CO2","^[a-zA-Z0-9]*$",str ],
"AO_CONTROL_USE"                 :["Backup, Overwrite, or Use Current CONTROL file", "BACKUP","^BACKUP$|^OVERWRITE$|^USE$|^backup$|^overwrite$|^use$",str],
"DL_POLY_CONTROL_TEMP"           :["Temperature for CONTROL file",273,"^[1-9]{1}$|^[0-9]+$",int],
"MAIN_GCMC_PACKAGE"              :["GCMC Package to use, DLPOLY or FASTGCMC","DLPOLY","^DLPOLY$|^FASTGCMC$",str],
"FAST_PP_G1"                     :["Partial pressure of first guest in fastGCMC",1.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"FAST_PP_G2"                     :["Partial pressure of second guest in fastGCMC",1.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"FAST_PP_G3"                     :["Partial pressure of third guest in fastGCMC",1.0,"^([0-9]*|\d*\.\d{1}?\d*)$",float],
"FAST_NUM_GUESTS"                :["Number of guest types in fastGCMC",1,"^[1-9]{1}$|^[0-9]+$",int],
"FAST_SIM_LENGTH"                :["Simulation length for FASTGCMC (timesteps)",1000000,"^[1-9]{1}$|^[0-9]+$",int ],
"FAST_EQ_TIME"                   :["Simulation equilibration time for FASTGCMC (units)",100000,"^[1-9]{1}$|^[0-9]+$",int ],
"MAIN_ENABLE_SCHEDULE"           :["Enable of disable job scheduler","DISABLE","^ENABLE$|^DISABLE$",str]

}



def init_variable(variable_name):
    if DEFAULTS.__contains__(variable_name) and VALUES.__contains__(variable_name):
        if len(VALUES[variable_name]) != 0:
            if re.match(str(DEFAULTS[variable_name][2]),str(VALUES[variable_name])) != None:
                return VALUES[variable_name]
            elif(type(VALUES[variable_name]) != DEFAULTS[variable_name][3]):
                print "Unexpected value for " +str(variable_name)+". Expecting type "+str(DEFAULTS[variable_name][3])+" for "+variable_name+" . Using default. Adjust regexp or check type."
            else:
                print "Unexpected value for " +str(variable_name)+". Using default. Adjust regexp or check value."
        print "No value found for variable " + variable_name+" . Exiting"
        sys.exit(1)
    elif not(DEFAULTS.__contains__(variable_name)):
        print "Error in JOB_INPUT file: Variable has no profile in JOB_OPTIONS | "+variable_name
	sys.exit(EXIT_FAILURE)
    elif not(VALUES.__contains__(variable_name)):
        print "Error in JOB_INPUT file: Variable not in JOB_INPUT file | "+variable_name
	sys.exit(EXIT_FAILURE)

def test_all_vars():
    for i in DEFAULTS.keys():
        print "%(A)55s | %(B)55s" % {'A':"Submitting "+str(i)+" results in ",'B': str(init_variable(i))}



if DEBUGGING:
    test_all_vars()
