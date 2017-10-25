The CNO directory builds a python network.

C_f90, CNO_f90, and URCA_f90 build Fortran 90 networks.

To make your own Fortran network, it will be useful to look at the
readme files for a problem in examples/ similar to your own, but here
are the main steps:

* Create a directory under examples/ (e.g. examples/working)

* Obtain Reaclib v.2 rate files for your problem and either place them
  in examples/working or pyreaclib/reaclib-rates.

* Write your own version of C_f90/c.py, URCA_f90/urca.py, or
CNO_f90/cno.py, replacing the list of 'files' with your own (e.g.
working.py)

* Choose either the 'sundials', 'python', or 'boxlib' types of output networks.

* Run your python script (pyreaclib must be in your PYTHONPATH).

```
$ python working.py
```

* (SUNDIALS) By default the sundials networks will output data based on the
cv_pars%DT_SAVE (timestep) and cv_pars%NDT_SAVE (number of time steps
to take) parameters in the net.par file. If you would like to add
custom stopping criteria, edit the FCVROOTFN subroutine in network.f90
so G(1) is the root to solve for.

* If you would like to include tabular rates, for now they must be in
the form of, eg. 23Na-23Ne_electroncapture.dat in examples/URCA_f90,
indexed by density*ye and temperature, with the same number and order
of variables. Prepare a rate file to describe how to read the table in
the form of, eg. URCA_f90/na23--ne23-toki:

```
t
       [parent nuclide]  [daughter nuclide]
[rate table file name]
[number of header lines to eat before the first line of data]
[number of density*ye values]
[number of temperature values]
```

* (SUNDIALS) Make

```
$ make
```

* (SUNDIALS) Specify initial abundances, density, temperature, ye in net.par

* (SUNDIALS) Run 

```
$ ./integrator.gfortran.exe
```