# Fortran: triple-alpha

Fortran network for triple-alpha fusion to C-12 and the reverse rate.

Now with the correct prefactor for C-12 photodisintegration.

To run:

1) Generate the reaction network using triple-alpha.py (pynucastro directory
must be in your PYTHONPATH):

```
$ python triple-alpha.py
```

2) Edit rhs.f90 and revise the stopping root solver condition in
subroutine FCVROOTFN if desired. If you don't want a stopping
condition, start with the recommended timing DT_SAVE and NDT_SAVE
above. For example, to stop when either He-4 or C-12 are depleted,
use the following:

```
    G(1) = min(Y(jhe4),Y(jc12))
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
