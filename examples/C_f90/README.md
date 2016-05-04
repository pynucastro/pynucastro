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

2) Edit network.f90 and fill in the binding energies (MeV) in the
init_net_info subroutine: (relative to proton BE)

```
    ebind_per_nucleon(jn)   = 0.782347d0
    ebind_per_nucleon(jp)   = 0.0d0
    ebind_per_nucleon(jhe4)   = 7.073915d0
    ebind_per_nucleon(jc12)   = 7.680144d0
    ebind_per_nucleon(jne20)   = 8.032240d0
    ebind_per_nucleon(jna23)   = 8.111493d0
    ebind_per_nucleon(jmg23)   = 7.901104d0
```

3) Edit rhs.f90 and revise the stopping root solver condition in
subroutine FCVROOTFN if desired. For example, to stop when C12 is
depleted, use the following:

```
    G(1) = Y(jc12)
```

4) Add your SUNDIALS, LAPACK, and BLAS libraries to ../GMake.common, eg.

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

5) Make and run

```
$ make
$ ./integrator.gfortran.exe
```
