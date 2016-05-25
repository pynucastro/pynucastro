# Fortran: C12-burning

The URCA-boxlib example creates a reaction network which implements
simple C12-burning. It is intended to generate a network for
BoxLib/Microphysics for use in MAESTRO or CASTRO.

To use:

1) Generate the reaction network using c.py (pyreaclib directory
must be in your PYTHONPATH):

```
$ python c.py
```

2) Edit actual_network.f90 and fill in the binding energies (MeV) in the
actual_network_init subroutine: (relative to proton BE)

```
    ebind_per_nucleon(jn)   = 0.782347d0
    ebind_per_nucleon(jp)   = 0.0d0
    ebind_per_nucleon(jhe4)   = 7.073915d0
    ebind_per_nucleon(jc12)   = 7.680144d0
    ebind_per_nucleon(jne20)   = 8.032240d0
    ebind_per_nucleon(jna23)   = 8.111493d0
    ebind_per_nucleon(jmg23)   = 7.901104d0
```

