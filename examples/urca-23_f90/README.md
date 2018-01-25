# Fortran: C12-burning with A=23 URCA reactions

The URCA_f90 example integrates a reaction network which implements
simple C12-burning along with the A=23 URCA reactions until either
23Ne or 23Na are depleted.

To run:

1) Generate the reaction network using urca.py (pynucastro directory
must be in your PYTHONPATH):

```
$ python urca.py
```

2) Edit `inputs_integration` to set the integration inital conditions
and timestepping.

3) Make and run

```
$ make
$ ./main.Linux.gfortran.exe inputs_integration
```

--------------------------------------------------------------------------------

The URCA example uses A=23 URCA rates from Suzuki et al., 2016, ApJ 817:163

They can be downloaded from the web at:

http://w3p.phys.chs.nihon-u.ac.jp/~suzuki/data2/link.html

(as of March 16, 2016)
