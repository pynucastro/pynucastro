# Fortran: CNO cycle

The CNO_f90 example uses the CVODE implementation of SUNDIALS to
integrate the CNO reaction network corresponding to the python-based
CNO example. The integration final time and abundance output save
cadence is set by the timing in net.par where DT_SAVE is the time
interval between file writes and NDT_SAVE is the number of file writes
to perform. The simulation end time is NDT_SAVE*DT_SAVE unless the
stopping-condition root solver in subroutine FCVROOTFN (network.f90)
finds a root.

```
  cv_pars%DT_SAVE  = 2.0 
  cv_pars%NDT_SAVE = 10000
```

To run:

1) Generate the reaction network using cno.py (pyreaclib directory
must be in your PYTHONPATH):

```
$ python cno.py
```

2) Edit net_rates.f90 and fill in the binding energies (MeV) in the
init_net_info subroutine:

```
    self%ebind_per_nucleon(self%ip) = 0.0d0
    self%ebind_per_nucleon(self%ihe4) = 7.073915d0
    self%ebind_per_nucleon(self%ic12) = 7.680144d0
    self%ebind_per_nucleon(self%ic13) = 7.469849d0
    self%ebind_per_nucleon(self%in13) = 7.238863d0
    self%ebind_per_nucleon(self%in14) = 7.475614d0
    self%ebind_per_nucleon(self%in15) = 7.699460d0
    self%ebind_per_nucleon(self%io14) = 7.052301d0
    self%ebind_per_nucleon(self%io15) = 7.46369d0
```

3) Edit network.f90 and revise the stopping root solver condition in
subroutine FCVROOTFN if desired. For example, to stop when H is
depleted, use the following:

```
    G(1) = Y(net_meta%ip)
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