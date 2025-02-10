PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 3

COMP = g++

BL_NO_FORT = TRUE

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE
USE_CXX_EOS = TRUE

USE_ECAP = TRUE
USE_BETA = TRUE

DO_A21 = FALSE
DO_A23 = TRUE
DO_A25 = FALSE

ifeq ($(USE_ECAP), FALSE)
  ifeq ($(DO_A21), FALSE)
    DEFINES += -DSKIP_A21_ECAP
  endif
  ifeq ($(DO_A23), FALSE)
    DEFINES += -DSKIP_A23_ECAP
  endif
  ifeq ($(DO_A25), FALSE)
    DEFINES += -DSKIP_A25_ECAP
  endif
endif

ifeq ($(USE_BETA), FALSE)
  ifeq ($(DO_A21), FALSE)
    DEFINES += -DSKIP_A21_BETA
  endif
  ifeq ($(DO_A23), FALSE)
    DEFINES += -DSKIP_A23_BETA
  endif
  ifeq ($(DO_A25), FALSE)
    DEFINES += -DSKIP_A25_BETA
  endif
endif

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
