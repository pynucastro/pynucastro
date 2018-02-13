# source files
$(warning Using SRCDIR: $(SRCDIR))
f90sources := $(foreach dir, $(SRCDIR), $(notdir $(wildcard $(dir)/*.f90)))
fsources := $(foreach dir, $(SRCDIR), $(notdir $(wildcard $(dir)/*.f)))
f90sources_dir += $(foreach dir, $(SRCDIR), $(wildcard $(dir)/*.f90))

$(warning Using f90sources: $(f90sources))
$(warning Using f90sources_dir: $(f90sources_dir))
$(warning Using fsources: $(fsources))
$(warning Using VPATH: $(VPATH))
$(warning Using program: $(program))

# dependencies
deps: $(f90sources_dir)
	../../../util/dep.py --program $(program) --vpath $(VPATH) --files $(f90sources_dir) --output deps 

include deps

$(warning Using f90sources: $(f90sources))

# Libraries to link
LINKLIBS :=

# set the compiler flags for those compilers we know about
ifeq ($(FC),gfortran)
  #FFLAGS := -c -O2 -g -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
  FFLAGS := -c -O2 -g -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero,overflow,underflow -finit-real=snan -ffree-line-length-none

  ifdef ACC
    FFLAGS += -fopenacc
  endif

  ifdef OMP
    FFLAGS += -fopenmp
  endif

  wrapper := $(FC)
  link := $(wrapper)

  # SUNDIALS libraries
  SUNLIBDIR := /home/eugene/local/sundials/instdir/lib
  LINKLIBS += -L${SUNLIBDIR} -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial

  # LAPACK
  LAPACKDIR := /usr/lib64
  LINKLIBS += -L${LAPACKDIR} -llapack

  # BLAS	 
  BLASDIR := /usr/lib64
  LINKLIBS += -L${BLASDIR} -lblas

else
  $(error ERROR: compiler $(FC) invalid)
endif

# default rule for building the object files
%.o:: %.f90
	$(wrapper) $(FFLAGS) -o $@ $<

%.o:: %.f
	$(wrapper) $(FFLAGS) -o $@ $<

# create the list of dependencies for the final build (all the .o files)
OBJECTS := $(f90sources:.f90=.o)
OBJECTS += $(fsources:.f=.o)

print-%: ; @echo $* is $($*)
