The CNO directory builds a python network. C_f90, CNO_f90, and
URCA_f90 build Fortran 90 networks.

To make your own Fortran network, copy
pyreaclib/templates/sundials-cvode to your own directory. Obtain
Reaclib v.2 rate files for your problem and write your own version of
C_f90/c.py, URCA_f90/urca.py, or CNO_f90/cno.py, replacing the list of
'files' with your own. Follow the steps in the URCA_f90 readme,
filling in the binding energy per nucleon in MeV.

Finally, by default the network will output data based on the
cv_pars%DT_SAVE (timestep) and cv_pars%NDT_SAVE (number of time steps
to take) parameters in the net.par file. If you would like to add
custom stopping criteria, edit the FCVROOTFN subroutine in network.f90
so G(1) is the root to solve for.

If you would like to include tabular rates, for now they must be in
the form of, eg. 23Na-23Ne_electroncapture.dat in the URCA_f90
example, indexed by density*ye and temperature, with the same number
and order of variables. Prepare a rate file to describe how to read
the table in the form of, eg. URCA_f90/na23--ne23-toki:

```
t
       [parent nuclide]  [daughter nuclide]
[rate table file name]
[number of header lines to eat before the first line of data]
[number of density*ye values]
[number of temperature values]
```

