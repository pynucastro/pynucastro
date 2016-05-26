# Fortran: C12-burning with A=23 URCA reactions

The URCA_f90 example integrates a reaction network which implements
simple C12-burning along with the A=23 URCA reactions until either
23Ne or 23Na are depleted.

The integration final time and abundance output save
cadence is set by the timing in net.par where DT_SAVE is the time
interval between file writes and NDT_SAVE is the number of file writes
to perform. The simulation end time is NDT_SAVE*DT_SAVE unless the
stopping-condition root solver in subroutine FCVROOTFN (network.f90)
finds a root.

If you don't add a stopping condition (step 3 is recommended), the
following timing is recommended to use in net.par:

```
  cv_pars%DT_SAVE  = 1.0d-7
  cv_pars%NDT_SAVE = 10000
```

To run:

1) Generate the reaction network using urca.py (pyreaclib directory
must be in your PYTHONPATH):

```
$ python urca.py
```

2) Edit rhs.f90 and revise the stopping root solver condition in
subroutine FCVROOTFN if desired. If you don't want a stopping
condition, start with the recommended timing DT_SAVE and NDT_SAVE
above. For example, to stop when either Ne-23 or Na-23 are depleted,
use the following:

```
    G(1) = min(Y(jne23),Y(jna23))
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

--------------------------------------------------------------------------------

The URCA example uses A=23 URCA rates from Suzuki et al., 2016, ApJ 817:163

They can be downloaded from the web at:

http://w3p.phys.chs.nihon-u.ac.jp/~suzuki/data2/link.html

(as of March 16, 2016)
