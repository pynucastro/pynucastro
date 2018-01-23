# A set of useful macros.

# include the main Makefile stuff
include $(PYNUCASTRO_HOME)/util/F_mk/GMakedefs.mak

PYNUCASTRO_CORE += $(PYNUCASTRO_HOME)/pynucastro/screening
PYNUCASTRO_CORE += $(PYNUCASTRO_HOME)/util/BLAS
PYNUCASTRO_CORE += $(PYNUCASTRO_HOME)/util/LINPACK
PYNUCASTRO_CORE += $(PYNUCASTRO_HOME)/util/VODE
PYNUCASTRO_CORE += $(PYNUCASTRO_HOME)/util/VODE/vode_source

#-----------------------------------------------------------------------------
# core pynucastro directories
Fmpack := $(foreach dir, $(PYNUCASTRO_CORE), $(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(PYNUCASTRO_CORE), $(dir))
Fmincs :=

# auxillary directories
Fmpack += $(foreach dir, $(PROGRAM_DIR), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(PROGRAM_DIR), $(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# we need a probin.f90, since the various microphysics routines can
# have runtime parameters
f90sources += probin.f90

PROBIN_TEMPLATE := $(PYNUCASTRO_HOME)/util/F_mk/dummy.probin.template
PROBIN_PARAMETER_DIRS += $(PYNUCASTRO_HOME)/pynucastro/templates/basefort
EXTERN_PARAMETER_DIRS += $(PYNUCASTRO_CORE)


PROBIN_PARAMETERS := $(shell $(PYNUCASTRO_HOME)/util/F_mk/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(PYNUCASTRO_HOME)/util/F_mk/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(PYNUCASTRO_HOME)/util/F_mk/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "


# vpath defines the directories to search for the source files

#  VPATH_LOCATIONS to first search in the problem directory
#  Note: GMakerules.mak will include '.' at the start of the
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(PYNUCASTRO_HOME)/util/F_mk/makebuildinfo.py \
           --modules "$(Fmdirs) $(PYNUCASTRO_CORE) $(UNIT_DIR)" \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --amrex_home "$(PYNUCASTRO_HOME)" \
           --source_home "$(PYNUCASTRO_HOME)" \
           --network "$(PYNUCASTRO_HOME)" \
           --integrator "$(PYNUCASTRO_HOME)" \
           --eos "$(PYNUCASTRO_HOME)"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean ::

realclean ::
