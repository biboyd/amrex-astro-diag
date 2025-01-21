PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 3

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE
USE_CXX_EOS = TRUE

MAX_ZONES := 16384

DEFINES += -DNPTS_MODEL=$(MAX_ZONES)

# programs to be compiled
EBASE := fneutrinos

# EOS and network
EOS_DIR := helmholtz

NETWORK_DIR := ignition_reaclib/URCA-simple
#NETWORK_INPUTS := triple_alpha_plus_o.net

Bpack := ./Make.package
Blocs := . ..

EXTERN_SEARCH += . ..

USE_AMR_CORE = TRUE

include $(MICROPHYSICS_HOME)/Make.Microphysics
