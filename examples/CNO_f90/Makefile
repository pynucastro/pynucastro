# Make the Integrator

SUNINCDIR = /home/eugene/local/sundials/instdir/include
SUNLIBDIR = /home/eugene/local/sundials/instdir/lib

F90 	= gfortran
FFLAGS 	= -O3 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -g
SUNINCLUDE = -I ${SUNINCDIR}
INCLUDE = 
F90LINK	= gfortran
CC	= gfortran

LLIBSS	= -L${SUNLIBDIR} -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial
LLIBSM	= -L/usr/lib64 -llapack -lblas

SRCS 	= network.f90 net_rates.f90 integrator.f90 physical_constants.f90 data_wrangler.f90 parameters.f90 cvode_parameters.f90 table_rates.f90
OBJS	= network.o net_rates.o cvode_parameters.o data_wrangler.o integrator.o parameters.o physical_constants.o table_rates.o

MAIN	= intnet.exe

all:	${MAIN}
	@echo	Compilation complete!

%.mod : %.o
	@echo Compiled module.

call_net : net_rates.mod physical_constants.mod
	${F90} ${FFLAGS} call_net.f90 net_rates.o physical_constants.o -o call_net.exe

call_table : table_rates.mod
	${F90} ${FFLAGS} call_table.f90 table_rates.o -o call_table.exe

output_table : table_rates.mod
	${F90} ${FFLAGS} output_table.f90 table_rates.o -o output_table.exe

table_rates.o : table_rates.f90
	${F90} ${FFLAGS} -c $< -o $@

integrator.o : integrator.f90 network.mod net_rates.mod data_wrangler.mod parameters.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

cvode_parameters.o : cvode_parameters.f90
	${F90} ${FFLAGS} -c $< -o $@

net_rates.o : net_rates.f90 physical_constants.mod table_rates.mod
	${F90} ${FFLAGS} -c $< -o $@

network.o : network.f90 net_rates.mod physical_constants.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

data_wrangler.o : data_wrangler.f90
	${F90} ${FFLAGS} -c $< -o $@

parameters.o : parameters.f90 cvode_parameters.mod net_rates.mod
	${F90} ${FFLAGS} -c $< -o $@

physical_constants.o : physical_constants.f90
	${F90} ${FFLAGS} -c $< -o $@

$(MAIN): ${OBJS}
	@echo '${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}'
	${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}

clean:
	rm *.o *.mod *.exe *~ ${MAIN}

depend: ${SRCS}
	makedepend ${INCLUDE} $^
