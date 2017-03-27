# Fortran: C12-burning with A=23 URCA reactions

The URCA-boxlib example creates a reaction network which implements
simple C12-burning along with the A=23 URCA reactions until either
23Ne or 23Na are depleted. It is intended to generate a network for
BoxLib/Microphysics for use in MAESTRO or CASTRO.

To use:

1) Generate the reaction network using urca.py (pyreaclib directory
must be in your PYTHONPATH):

```
$ python urca.py
```

2) Edit actual_network.f90 and fill in the binding energies (MeV) in the
actual_network_init subroutine: (relative to proton BE)

```
    ebind_per_nucleon(jn)   = 0.782347d0
    ebind_per_nucleon(jp)   = 0.0d0
    ebind_per_nucleon(jhe4)   = 7.073915d0
    ebind_per_nucleon(jc12)   = 7.680144d0
    ebind_per_nucleon(jne20)   = 8.032240d0
    ebind_per_nucleon(jne23)   = 7.955255d0
    ebind_per_nucleon(jna23)   = 8.111493d0
    ebind_per_nucleon(jmg23)   = 7.901104d0
```

--------------------------------------------------------------------------------

The URCA example uses A=23 URCA rates from Suzuki et al., 2016, ApJ 817:163

They can be downloaded from the web at:

http://w3p.phys.chs.nihon-u.ac.jp/~suzuki/data2/link.html

(as of March 16, 2016)
