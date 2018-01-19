# Fortran: C12-burning with A=23 URCA reactions

The URCA-starkiller example creates a reaction network which implements
simple C12-burning along with the A=23 URCA reactions until either
23Ne or 23Na are depleted. It is intended to generate a network for
StarKiller/Microphysics for use in MAESTRO or CASTRO.

To use, generate the reaction network using urca.py (pynucastro directory
must be in your PYTHONPATH):

```
$ python urca.py
```

--------------------------------------------------------------------------------

The URCA example uses A=23 URCA rates from Suzuki et al., 2016, ApJ 817:163

They can be downloaded from the web at:

http://w3p.phys.chs.nihon-u.ac.jp/~suzuki/data2/link.html

(as of March 16, 2016)
