PRECISION  = DOUBLE
DIM        = 3
COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE

EBASE = main

BL_NO_FORT = TRUE

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  ?= ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR = custom
NETWORK_CUSTOM_PATH := $(shell pwd)

CUSTOM_PARAMETERS := _standalone_parameters

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := . $(MICROPHYSICS_HOME)/unit_test

CEXE_sources += main.cpp

include $(MICROPHYSICS_HOME)/Make.Microphysics


