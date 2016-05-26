# Fortran: C12-burning

The C_f90 example integrates a reaction network which implements
simple C12-burning until C-12 is depleted. The integration final time
and abundance output save cadence is set by the timing in net.par
where DT_SAVE is the time interval between file writes and NDT_SAVE is
the number of file writes to perform. The simulation end time is
NDT_SAVE*DT_SAVE unless the stopping-condition root solver in
subroutine FCVROOTFN (network.f90) finds a root.

```
  cv_pars%DT_SAVE  = 2.0 
  cv_pars%NDT_SAVE = 10000
```

To run:

1) Generate the reaction network using c.py (pyreaclib directory must
be in your PYTHONPATH):

```
$ python c.py
```

2) Edit rhs.f90 and revise the stopping root solver condition in
subroutine FCVROOTFN if desired. For example, to stop when C12 is
depleted, use the following:

```
    G(1) = Y(jc12)
```

3) Add your SUNDIALS, LAPACK, and BLAS libraries to ../GMake.common, eg.

```
  # SUNDIALS libraries
  SUNLIBDIR := /home/eugene/local/sundials/instdir/lib
  LINKLIBS += -L${SUNLIBDIR} -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial

  # LAPACK
  LAPACKDIR := /usr/lib64
  LINKLIBS += -L${LAPACKDIR} -llapack

  # BLAS	 
  BLASDIR := /usr/lib64
  LINKLIBS += -L${BLASDIR} -lblas
```

4) Make and run

```
$ make
$ ./integrator.gfortran.exe
```
